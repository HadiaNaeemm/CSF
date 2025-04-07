close all; clear all; clc;
avgnum=100;
OMEGA = 1:1:10; %column
AMPLITUDEavg=0.01:0.01:0.1; %row
% Initialize arrays to accumulate errors over 'avgnum' iterations
Error_nonzero_x_avg = zeros(length(OMEGA), length(AMPLITUDEavg)); % Accumulate nonzero error
Error_zero_x_avg = zeros(length(OMEGA), length(AMPLITUDEavg));    % Accumulate zero error

Error_nonzero_y_avg = zeros(length(OMEGA), length(AMPLITUDEavg)); % Accumulate nonzero error
Error_zero_y_avg = zeros(length(OMEGA), length(AMPLITUDEavg));    % Accumulate zero error

for avg=1:avgnum
    % Initialize arrays to store errors for each 'AMPLITUDEavg(num)'
    Error_nonzero_x_all = zeros(length(OMEGA), length(AMPLITUDEavg)); % Store nonzero error for each num
    Error_zero_x_all = zeros(length(OMEGA), length(AMPLITUDEavg));    % Store zero error for each num

    Error_nonzero_y_all = zeros(length(OMEGA), length(AMPLITUDEavg)); % Store nonzero error for each num
    Error_zero_y_all = zeros(length(OMEGA), length(AMPLITUDEavg));    % Store zero error for each num

    for numa=1:length(AMPLITUDEavg)
         for numw=1:length(OMEGA)
        %initializing
        dt=0.01;
        x(1)=0.5;
        y(1)=0.5;
        %value of constants
        omega=OMEGA(numw);
        a=AMPLITUDEavg(numa);
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
    N_m=2;
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
    cd '/home/hadia/cvx/SLSO/R1'
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
    ksaix
    ksaiy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Error_nonzero_x=0;%initialization
    Error_zero_x=0;%initialization
    Num_nz=0;%initialization
    Num_z=0;%initialization
    for n=1:16
         if  n==2||n==3||n==9||n==10
             Error_nonzero_x=Error_nonzero_x+abs((ksaix(n)-ksaix_ref(n))/abs(ksaix_ref(n)));
             Num_nz=Num_nz+1;
         else
           Error_zero_x=Error_zero_x+abs(ksaix(n));
            Num_z=Num_z+1;
         end  
     end
    % Store errors for the current AMPLITUDEavg(num)
        Error_nonzero_x_all(numw,numa) = Error_nonzero_x/Num_nz;
        Error_zero_x_all(numw,numa) = Error_zero_x/Num_z;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Error_nonzero_y=0;%initialization
    Error_zero_y=0;%initialization
    Num_nz=0;%initialization
    Num_z=0;%initialization
    for n=1:16
         if  n==2||n==3||n==8||n==11
             Error_nonzero_y=Error_nonzero_y+abs((ksaiy(n)-ksaiy_ref(n))/abs(ksaiy_ref(n)));
             Num_nz=Num_nz+1;
         else
           Error_zero_y=Error_zero_y+abs(ksaiy(n));
            Num_z=Num_z+1;
         end  
     end
    % Store errors for the current AMPLITUDEavg(num)
        Error_nonzero_y_all(numw,numa) = Error_nonzero_y/Num_nz;
        Error_zero_y_all(numw,numa) = Error_zero_y/Num_z;
    end
    end 
    % Accumulate errors over 'avgnum' iterations
    Error_nonzero_x_avg = Error_nonzero_x_avg + Error_nonzero_x_all;
    Error_zero_x_avg = Error_zero_x_avg + Error_zero_x_all;
    
    Error_nonzero_y_avg = Error_nonzero_y_avg + Error_nonzero_y_all;
    Error_zero_y_avg = Error_zero_y_avg + Error_zero_y_all;
    
end
% Calculate the average errors by dividing by 'avgnum'
Error_nonzero_x_avg = Error_nonzero_x_avg / avgnum;
Error_zero_x_avg = Error_zero_x_avg / avgnum;
Error_nonzero_y_avg = Error_nonzero_y_avg / avgnum;
Error_zero_y_avg = Error_zero_y_avg / avgnum;


saveDir = '/home/hadia/cvx/SLSO/R1';  
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% Save each error matrix as a separate .mat file
% File paths for each error type
file_nonzero_x = fullfile(saveDir, 'Error_nonzero_x_avg.mat');
file_zero_x = fullfile(saveDir, 'Error_zero_x_avg.mat');
file_nonzero_y = fullfile(saveDir, 'Error_nonzero_y_avg.mat');
file_zero_y = fullfile(saveDir, 'Error_zero_y_avg.mat');

% Save matrices
save(file_nonzero_x, 'Error_nonzero_x_avg');
save(file_zero_x, 'Error_zero_x_avg');
save(file_nonzero_y, 'Error_nonzero_y_avg');
save(file_zero_y, 'Error_zero_y_avg');
