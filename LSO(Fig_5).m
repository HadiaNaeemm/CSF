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

Realization=100;
for r=1:Realization

Error_nonzero_x=zeros([1,Ndatamax-Ndatamin+1]);%initialization
Error_zero_x=zeros([1,Ndatamax-Ndatamin+1]);%initialization
Error_nonzero_y=zeros([1,Ndatamax-Ndatamin+1]);%initialization
Error_zero_y=zeros([1,Ndatamax-Ndatamin+1]);%initialization

count=0;
  for N_m=Ndatamin:Ndatamax
    count=count+1;
    N_measurements=N_m;
    N_basis=16;
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
                M(i,j)=x(index(i))^2*y(index(i))^2;
            elseif j==8
                M(i,j)=x(index(i))^2*y(index(i));
            elseif j==9
                M(i,j)=x(index(i))*y(index(i))^2;
            elseif j==10
                M(i,j)=x(index(i))^3;
            elseif j==11
                M(i,j)=y(index(i))^3;  
            elseif j==12
                M(i,j)=x(index(i))^3*y(index(i));  
            elseif j==13
                M(i,j)=x(index(i))*y(index(i))^3;
            elseif j==14
                M(i,j)=x(index(i))^3*y(index(i))^2;  
            elseif j==15
                M(i,j)=x(index(i))^2*y(index(i))^3;  
            else
                M(i,j)=x(index(i))^3*y(index(i))^3;  
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
    ksaix_ref=[0 a -omega 0 0 0 0 0 -1 -1 0 0 0 0 0 0];%nonzero:a -omega -1 -1
    ksaiy_ref=[0 omega a  0 0 0 0 -1 0 0 -1 0 0 0 0 0];%nonzero:omega a  -1 -1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       Reconstruction                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda=0.001;
    Xdot=Xdot';
    cd '/home/hadia/cvx/LSO_1/L1'
    cvx_setup
    cvx_begin 
      variable ksaix(16); 
      minimize(norm(M*ksaix-Xdot)+ lambda*norm(ksaix,1)); 
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ydot=Ydot';
    cvx_setup
    cvx_begin 
      variable ksaiy(16); 
      minimize(norm(M*ksaiy-Ydot)+ lambda*norm(ksaiy,1)); 
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       RIP                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:N_basis
        ksaix(j)=ksaix(j)/Norms(j);
        ksaiy(j)=ksaiy(j)/Norms(j);
    end
    ksaix;
    ksaiy;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Error_nonzero_x(count)=0;%initialization
    Error_zero_x(count)=0;%initialization
    Num_nz=0;%initialization
    Num_z=0;%initialization
    for n=1:16
         if  n==2||n==3||n==9||n==10
             Error_nonzero_x(count)=Error_nonzero_x(count)+abs((ksaix(n)-ksaix_ref(n))/abs(ksaix_ref(n)));
             Num_nz=Num_nz+1;
         else
           Error_zero_x(count)=Error_zero_x(count)+abs(ksaix(n));
            Num_z=Num_z+1;
         end  
     end
    % Store errors for the current AMPLITUDEavg(num)
        Error_nonzero_x(count) = Error_nonzero_x(count)/Num_nz
        Error_zero_x(count) = Error_zero_x(count)/Num_z
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Error_nonzero_y(count)=0;%initialization
    Error_zero_y(count)=0;%initialization
    Num_nz=0;%initialization
    Num_z=0;%initialization
    for n=1:16
         if  n==2||n==3||n==8||n==11
             Error_nonzero_y(count)=Error_nonzero_y(count)+abs((ksaiy(n)-ksaiy_ref(n))/abs(ksaiy_ref(n)));
             Num_nz=Num_nz+1;
         else
           Error_zero_y(count)=Error_zero_y(count)+abs(ksaiy(n));
            Num_z=Num_z+1;
         end  
     end
        Error_nonzero_y(count) = Error_nonzero_y(count)/Num_nz
        Error_zero_y(count) = Error_zero_y(count)/Num_z

  end
% Define the directory where you want to save the files
save_dir = '/home/hadia/cvx/LSO_1/L1';

% Ensure the directory exists (create it if it doesn’t)
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% Save results after each realization
save(fullfile(save_dir, sprintf('Error_nonzero_x_realization_%d.mat', r)), 'Error_nonzero_x');
save(fullfile(save_dir, sprintf('Error_zero_x_realization_%d.mat', r)), 'Error_zero_x');
save(fullfile(save_dir, sprintf('Error_nonzero_y_realization_%d.mat', r)), 'Error_nonzero_y');
save(fullfile(save_dir, sprintf('Error_zero_y_realization_%d.mat', r)), 'Error_zero_y');

end 

