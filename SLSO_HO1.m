close all; clear all; clc;
OMEGA = 3;
AMPLITUDEavg=0.01;
        dt=0.01;
        x(1)=0.5;
        y(1)=0.5;
        %value of constants
        omega= OMEGA ;
        a=AMPLITUDEavg;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:10000
    x(i)=x(i-1)+((a-x(i-1)^2-y(i-1)^2)*x(i-1)-omega*y(i-1))*dt;
    y(i)=y(i-1)+((a-x(i-1)^2-y(i-1)^2)*y(i-1)+omega*x(i-1))*dt;
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Trim the initial transient part
     x = x(1000:10000);
     y = y(1000:10000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Observation                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ndatamin=1;Ndatamax=20;
Data=Ndatamin:Ndatamax;
 for N_m=Ndatamin:Ndatamax
    N_measurements=N_m;
    N_basis=15;
    index=randi([100,499],1,N_m);
    Xdot=zeros([1,N_measurements]);
    Ydot=zeros([1,N_measurements]);
for ni=1:N_measurements
        Xdot(ni)=(x(index(ni)+1)-x(index(ni)))/dt;
        Ydot(ni)=(y(index(ni)+1)-y(index(ni)))/dt;
end
    M=zeros([N_measurements,N_basis]);
    for i=1:N_measurements
        for j=1:N_basis
            if j==1
               M(i,j)=1;
            elseif j==2
               M(i,j)=x(index(i));
            elseif j==3
               M(i,j)=y(index(i));
            elseif j==4
               M(i,j)=x(index(i))*y(index(i));
            elseif j==5
               M(i,j)=x(index(i))^2;
            elseif j==6
               M(i,j)=y(index(i))^2;
            elseif j==7
               M(i,j)=x(index(i))^3;
            elseif j==8
               M(i,j)=y(index(i))^3;
            elseif j==9
               M(i,j)=x(index(i))^2*y(index(i));
            elseif j==10
               M(i,j)=x(index(i))*y(index(i))^2;
            elseif j==11
               M(i,j)=x(index(i))^4;
            elseif j==12
               M(i,j)=y(index(i))^4;
             elseif j==13
               M(i,j)=x(index(i))^3*y(index(i));
             elseif j==14
               M(i,j)=x(index(i))*y(index(i))^3;
            else
               M(i,j)=x(index(i))^2*y(index(i))^2;
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       Reconstruction                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda=0.001;
    Xdot=Xdot';
    cd 'D:\MATLAB\cvx-w64\cvx'
    cvx_setup
    cvx_begin 
      variable ksaix(15); 
      minimize(norm(M*ksaix-Xdot)+ lambda*norm(ksaix,1)); 
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       RIP                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:N_basis
        ksaix(j)=ksaix(j)/Norms(j);
    end
    ksaix;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end 