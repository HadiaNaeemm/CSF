clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Simulation                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%value of constants
a1=0.1;a2=0.2;
omega1=5;omega2=4;
G=1;C12=0.01;C21=0.02;
dt=0.01;   %step size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1(1)=0.5;
y1(1)=0.5;
x2(1)=1;
y2(1)=1;
for i=2:1000
    x1(i)=x1(i-1)+((a1-x1(i-1)^2-y1(i-1)^2)*x1(i-1)-omega1*y1(i-1)+G*C12*(x2(i-1)-x1(i-1)))*dt;
    y1(i)=y1(i-1)+((a1-x1(i-1)^2-y1(i-1)^2)*y1(i-1)+omega1*x1(i-1)+G*C12*(y2(i-1)-y1(i-1)))*dt;
    x2(i)=x2(i-1)+((a2-x2(i-1)^2-y2(i-1)^2)*x2(i-1)-omega2*y2(i-1)+G*C21*(x1(i-1)-x2(i-1)))*dt;
    y2(i)=y2(i-1)+((a2-x2(i-1)^2-y2(i-1)^2)*y2(i-1)+omega2*x2(i-1)+G*C21*(y1(i-1)-y2(i-1)))*dt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Ofservation                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ndatamin=1;Ndatamax=20;
Data=Ndatamin:Ndatamax;
  for N_m=Ndatamin:Ndatamax
    N_measurements=N_m;
    N_basis=11;
    index=randi([100,499],1,N_measurements);
    Xdot1=zeros([1,N_measurements]);
    Ydot1=zeros([1,N_measurements]);
    Xdot2=zeros([1,N_measurements]);
    Ydot2=zeros([1,N_measurements]);
    for ni=1:N_measurements
        Xdot1(ni)=(x1(index(ni)+1)-x1(index(ni)))/dt;
        Ydot1(ni)=(y1(index(ni)+1)-y1(index(ni)))/dt;
        Xdot2(ni)=(x2(index(ni)+1)-x2(index(ni)))/dt;
        Ydot2(ni)=(y2(index(ni)+1)-y2(index(ni)))/dt;
    end
    M=zeros([N_measurements,N_basis]);
    for i=1:N_measurements
        for j=1:N_basis 
            if j==1
               M(i,j)=1;
            elseif j==2
               M(i,j)=x1(index(i));
            elseif j==3
               M(i,j)=y1(index(i));
            elseif j==4
               M(i,j)=x2(index(i));
            elseif j==5
               M(i,j)=x1(index(i))*y1(index(i));
            elseif j==6
               M(i,j)=x1(index(i))^2;
            elseif j==7
               M(i,j)=y1(index(i))^2;
            elseif j==8
               M(i,j)=x1(index(i))^3;
            elseif j==9
               M(i,j)=y1(index(i))^3;
            elseif j==10
               M(i,j)=x1(index(i))^2*y1(index(i));
            else
               M(i,j)=x1(index(i))*y1(index(i))^2;
            
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       RIP                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Norms=zeros([1,N_basis]);
    for j=1:N_basis
        Norms(j)=norm(M(:,j));
    end
    for i=1:N_measurements
        for j=1:N_basis
            M(i,j)=M(i,j)/Norms(j);
        end
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda=0.001;
    Xdot1=Xdot1';
    cd 'D:\MATLAB\cvx-w64\cvx'
    cvx_setup
    cvx_begin 
      variable ksaix1(11); 
      minimize(norm(M*ksaix1-Xdot1)+ lambda*norm(ksaix1,1)); 
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Xdot2=Xdot2';
    cvx_setup
    cvx_begin  
      variable ksaix2(11); 
      minimize(norm(M*ksaix2-Xdot2)+ lambda*norm(ksaix2,1)); 
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       RIP                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:N_basis
        ksaix1(j)=ksaix1(j)/Norms(j);
    end
    ksaix1;
  end



