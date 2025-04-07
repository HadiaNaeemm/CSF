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
Ndatamin=1;Ndatamax=30;
Data=Ndatamin:Ndatamax;
Realization=25;
for r=1:Realization

Error_nonzero_x1=zeros([1,Ndatamax-Ndatamin+1]);%initialization
Error_zero_x1=zeros([1,Ndatamax-Ndatamin+1]);%initialization
Error_nonzero_x2=zeros([1,Ndatamax-Ndatamin+1]);%initialization
Error_zero_x2=zeros([1,Ndatamax-Ndatamin+1]);%initialization
count=0;
  for N_m=Ndatamin:Ndatamax
    count=count+1;
    N_measurements=N_m;
    N_basis=20;
    index=randi([100,999],1,N_measurements);
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
               M(i,j)=y2(index(i));
            elseif j==6
               M(i,j)=x1(index(i))*y1(index(i))^2;
            elseif j==7
                M(i,j)=x2(index(i))*y2(index(i))^2;
            elseif j==8
                M(i,j)=x1(index(i))^3;
            elseif j==9
                M(i,j)=x2(index(i))^3;
            elseif j==10
                M(i,j)=y1(index(i))^2;
            elseif j==11
                M(i,j)=y2(index(i))^2;
            elseif j==12
                M(i,j)=y1(index(i))^3;  
            elseif j==13
                M(i,j)=y2(index(i))^3;
            elseif j==14
                M(i,j)=y2(index(i))^2*y1(index(i));
            elseif j==15
                M(i,j)=y1(index(i))^2*x2(index(i));
            elseif j==16
                M(i,j)=y1(index(i))^2*y2(index(i));
            elseif j==17
                M(i,j)=x1(index(i))^2*y1(index(i));
            elseif j==18
                M(i,j)=y2(index(i))^2*x1(index(i)); 
            elseif j==19
                M(i,j)=x2(index(i))^2*x1(index(i));
            else
                M(i,j)=x2(index(i))^2*y1(index(i));
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ksaix1_ref=[0 (a1-G*C12) -omega1 G*C12 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0];
    ksaix2_ref=[0 G*C21 0 (a2-G*C21) -omega2 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0];
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
    cd '/home/hadia/cvx/CLSON/L1'
    cvx_setup
    cvx_begin 
      variable ksaix1(20); 
      minimize(norm(M*ksaix1-Xdot1)+ lambda*norm(ksaix1,1)); 
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Xdot2=Xdot2';
    cvx_setup
    cvx_begin  
      variable ksaix2(20); 
      minimize(norm(M*ksaix2-Xdot2)+ lambda*norm(ksaix2,1)); 
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       RIP                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:N_basis
        ksaix1(j)=ksaix1(j)/Norms(j);
        ksaix2(j)=ksaix2(j)/Norms(j);
    end
    ksaix1;
    ksaix2;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Error_nonzero_x1(count)=0;%initialization
    Error_zero_x1(count)=0;%initialization
    Num_nz=0;%initialization
    Num_z=0;%initialization
    for n=1:20
         if  n==2||n==3||n==4||n==6||n==8
             Error_nonzero_x1(count)=Error_nonzero_x1(count)+abs((ksaix1(n)-ksaix1_ref(n))/abs(ksaix1_ref(n)));
             Num_nz=Num_nz+1;
         else
           Error_zero_x1(count)=Error_zero_x1(count)+abs(ksaix1(n));
            Num_z=Num_z+1;
         end  
     end
    % Store errors for the current AMPLITUDEavg(num)
        Error_nonzero_x1(count) = Error_nonzero_x1(count)/Num_nz;
        Error_zero_x1 (count)= Error_zero_x1(count)/Num_z;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Error_nonzero_x2(count)=0;%initialization
    Error_zero_x2(count)=0;%initialization
    Num_nz=0;%initialization
    Num_z=0;%initialization
    for n=1:20
         if  n==2||n==4||n==5||n==7||n==9
             Error_nonzero_x2(count)=Error_nonzero_x2(count)+abs((ksaix2(n)-ksaix2_ref(n))/abs(ksaix2_ref(n)));
             Num_nz=Num_nz+1;
         else
           Error_zero_x2(count)=Error_zero_x2(count)+abs(ksaix2(n));
            Num_z=Num_z+1;
         end  
     end
    % Store errors for the current AMPLITUDEavg(num)
        Error_nonzero_x2(count) = Error_nonzero_x2(count)/Num_nz;
        Error_zero_x2 (count)= Error_zero_x2(count)/Num_z;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end 
  % Define the directory where you want to save the files
save_dir = '/home/hadia/cvx/CLSON/L1';

% Ensure the directory exists (create it if it doesnâ€™t)
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% Save results after each realization
save(fullfile(save_dir, sprintf('Error_nonzero_x1_realization_%d.mat', r)), 'Error_nonzero_x1');
save(fullfile(save_dir, sprintf('Error_zero_x1_realization_%d.mat', r)), 'Error_zero_x1');
save(fullfile(save_dir, sprintf('Error_nonzero_x2_realization_%d.mat', r)), 'Error_nonzero_x2');
save(fullfile(save_dir, sprintf('Error_zero_x2_realization_%d.mat', r)), 'Error_zero_x2');
end 



