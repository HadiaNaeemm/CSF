% Clear workspace, close figures, and clear command window
close all; clear; clc;
% Parameters of the Lorenz system
sigma = 10;     
rho = 28;       
beta = 8/3;     

% Coupling strength
d1 = 0.1;    
d2 = 0.2;
% Time step and total iterations
dt = 0.01;      
nSteps = 10000; 

% Initial conditions for both oscillators
x1(1) = 0.5; y1(1) = 0.5; z1(1) = 0.5;
x2(1) = -0.5; y2(1) = -0.5; z2(1) = -0.5;

% Time evolution using an explicit iterative approach
for i = 2:nSteps
    x1(i) = x1(i-1) + (sigma * (y1(i-1) - x1(i-1)) + d1 * (x2(i-1) - x1(i-1))) * dt;
    y1(i) = y1(i-1) + (x1(i-1) * (rho - z1(i-1)) - y1(i-1)) * dt;
    z1(i) = z1(i-1) + (x1(i-1) * y1(i-1) - beta * z1(i-1)) * dt;

    x2(i) = x2(i-1) + (sigma * (y2(i-1) - x2(i-1)) + d2 * (x1(i-1) - x2(i-1))) * dt;
    y2(i) = y2(i-1) + (x2(i-1) * (rho - z2(i-1)) - y2(i-1)) * dt;
    z2(i) = z2(i-1) + (x2(i-1) * y2(i-1) - beta * z2(i-1)) * dt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Observation                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    N_basis=12;
    index=randi([100,499],1,N_m);
    X1dot=zeros([1,N_measurements]);
    Y1dot=zeros([1,N_measurements]);
    Z1dot=zeros([1,N_measurements]);
    X2dot=zeros([1,N_measurements]);
    Y2dot=zeros([1,N_measurements]);
    Z2dot=zeros([1,N_measurements]);
for ni=1:N_measurements
        X1dot(ni)=(x1(index(ni)+1)-x1(index(ni)))/dt;
        Y1dot(ni)=(y1(index(ni)+1)-y1(index(ni)))/dt;
        Z1dot(ni)=(z1(index(ni)+1)-z1(index(ni)))/dt;
        X2dot(ni)=(x2(index(ni)+1)-x2(index(ni)))/dt;
        Y2dot(ni)=(y2(index(ni)+1)-y2(index(ni)))/dt;
        Z2dot(ni)=(z2(index(ni)+1)-z2(index(ni)))/dt;
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
                M(i,j)=z1(index(i));
            elseif j==7
                M(i,j)=z2(index(i));
            elseif j==8
                M(i,j)=x1(index(i))^2;
            elseif j==9
               M(i,j)=x1(index(i))*y1(index(i));
            elseif j==10
               M(i,j)=x1(index(i))*z1(index(i));
            elseif j==11
                M(i,j)=y1(index(i))*z1(index(i));
            else
              M(i,j)=x1(index(i))*y1(index(i))*z1(index(i));    
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
    ksaix1_ref=[0 -(sigma+d1) sigma d1 0 0 0 0 0 0 0 0];
    ksaix2_ref=[0 d2 0 -(sigma+d2) sigma 0 0 0 0 0 0 0];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       Reconstruction                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda=0.001;
    X1dot=X1dot';
    cd '/home/hadia/cvx/CLON/L1'
    cvx_setup
    cvx_begin 
      variable ksaix1(12); 
      minimize(norm(M*ksaix1-X1dot)+ lambda*norm(ksaix1,1)); 
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda=0.001;
    X2dot=X2dot';
    cd '/home/hadia/cvx/CLON/L1'
    cvx_setup
    cvx_begin 
      variable ksaix2(12); 
      minimize(norm(M*ksaix2-X2dot)+ lambda*norm(ksaix2,1)); 
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    for n=1:12
         if  n==2||n==3||n==4
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
    for n=1:12
         if  n==2||n==4||n==5
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
save_dir = '/home/hadia/cvx/CLON/L1';

% Ensure the directory exists (create it if it doesn’t)
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% Save results after each realization
save(fullfile(save_dir, sprintf('Error_nonzero_x1_realization_%d.mat', r)), 'Error_nonzero_x1');
save(fullfile(save_dir, sprintf('Error_zero_x1_realization_%d.mat', r)), 'Error_zero_x1');
save(fullfile(save_dir, sprintf('Error_nonzero_x2_realization_%d.mat', r)), 'Error_nonzero_x2');
save(fullfile(save_dir, sprintf('Error_zero_x2_realization_%d.mat', r)), 'Error_zero_x2');
end 

