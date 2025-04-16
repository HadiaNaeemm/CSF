% Clear workspace, close figures, and clear command window
close all; clear all; clc;

% Parameters of the Lorenz system
sigma = 10;     % Prandtl number
rho = 28;       % Rayleigh number
beta = 8/3;     % Geometric factor

% Time step and total iterations
dt = 0.01;      % Time step
nSteps = 10000; % Number of time steps

% Initial conditions
x(1) = 0.5; % Initial x value
y(1) = 0.5; % Initial y value
z(1) = 0.5; % Initial z value

% Time evolution using an explicit iterative approach
for i = 2:nSteps
    x(i) = x(i-1) + sigma * (y(i-1) - x(i-1)) * dt;
    y(i) = y(i-1) + (x(i-1) * (rho - z(i-1)) - y(i-1)) * dt;
    z(i) = z(i-1) + (x(i-1) * y(i-1) - beta * z(i-1))  * dt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Observation                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ndatamin=1;Ndatamax=20;
Data=Ndatamin:Ndatamax;

Realization=100;
for r=1:Realization

Error_nonzero_y=zeros([1,Ndatamax-Ndatamin+1]);%initialization
Error_zero_y=zeros([1,Ndatamax-Ndatamin+1]);%initialization

count=0;
  for N_m=Ndatamin:Ndatamax
    count=count+1;
    N_measurements=N_m;
    N_basis=12;
    index=randi([100,499],1,N_m);
    Xdot=zeros([1,N_measurements]);
    Ydot=zeros([1,N_measurements]);
    Zdot=zeros([1,N_measurements]);
for ni=1:N_measurements
        Xdot(ni)=(x(index(ni)+1)-x(index(ni)))/dt;
        Ydot(ni)=(y(index(ni)+1)-y(index(ni)))/dt;
        Zdot(ni)=(z(index(ni)+1)-z(index(ni)))/dt;
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
               M(i,j)=z(index(i)); 
            elseif j==5
               M(i,j)=x(index(i))*y(index(i));
            elseif j==6
               M(i,j)=x(index(i))*z(index(i));
            elseif j==7
                M(i,j)=y(index(i))*z(index(i));
           elseif j==8
              M(i,j)=x(index(i))*y(index(i))*z(index(i));
            elseif j==9
                M(i,j)=x(index(i))^2;
            elseif j==10
                M(i,j)=y(index(i))^2;
            elseif j==11
                M(i,j)=z(index(i))^2;
            else
                M(i,j)=x(index(i))^2*y(index(i));  
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
    ksaix_ref=[0 -sigma sigma 0 0 0 0 0 0 0 0 0];
    ksaiy_ref=[0  rho  -1  0  0  -1 0 0 0  0 0 0];
    ksaiz_ref=[0  0  0  -beta  1  0 0 0 0  0 0 0];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       Reconstruction                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda=0.001;
    Xdot=Xdot';
    cd '/home/hadia/cvx/LOY_1/L1'
    cvx_setup
    cvx_begin 
      variable ksaix(12); 
      minimize(norm(M*ksaix-Xdot)+ lambda*norm(ksaix,1)); 
    cvx_end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Ydot=Ydot';
    cvx_setup
    cvx_begin 
      variable ksaiy(12); 
      minimize(norm(M*ksaiy-Ydot)+ lambda*norm(ksaiy,1)); 
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       RIP                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:N_basis
        ksaiy(j)=ksaiy(j)/Norms(j);
    end
    ksaiy;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Error_nonzero_y(count)=0;%initialization
    Error_zero_y(count)=0;%initialization
    Num_nz=0;%initialization
    Num_z=0;%initialization
    for n=1:12
         if n==2||n==3||n==6
             Error_nonzero_y(count)=Error_nonzero_y(count)+abs((ksaiy(n)-ksaiy_ref(n))/abs(ksaiy_ref(n)));
             Num_nz=Num_nz+1;
         else
           Error_zero_y(count)=Error_zero_y(count)+abs(ksaiy(n));
            Num_z=Num_z+1;
         end  
     end
    % Store errors for the current AMPLITUDEavg(num)
        Error_nonzero_y(count) = Error_nonzero_y(count)/Num_nz;
        Error_zero_y(count) = Error_zero_y(count)/Num_z;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
% Define the directory where you want to save the files
save_dir = '/home/hadia/cvx/LOY_1/L1';

% Ensure the directory exists (create it if it doesn’t)
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% Save results after each realization
save(fullfile(save_dir, sprintf('Error_nonzero_y_realization_%d.mat', r)), 'Error_nonzero_y');
save(fullfile(save_dir, sprintf('Error_zero_y_realization_%d.mat', r)), 'Error_zero_y');

end 




