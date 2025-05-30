close all;
clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
avgnum=100;
%value of constants
G=1;C12=0.02;C21=0.05;
dt=0.01;   %step size
OMEGA_bar = 1:1:10;
delta_omega_range = -0.5:0.02:0.5;
AMPLITUDE=0.001:0.001:0.009;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize arrays to accumulate errors over 'avgnum' iterations
% Error_nonzero_x11_avg = zeros(length(OMEGA_bar), length(AMPLITUDE), length(delta_omega_range)); % Accumulate nonzero error
% Error_zero_x11_avg = zeros(length(OMEGA_bar), length(AMPLITUDE), length(delta_omega_range));    % Accumulate zero error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for avg=1:avgnum
    % Initialize arrays to store errors for each 'AMPLITUDEavg(num)'
    Error_nonzero_x11_all = zeros(length(OMEGA_bar), length(AMPLITUDE), length(delta_omega_range)); % Store nonzero error for each num
    Error_zero_x11_all = zeros(length(OMEGA_bar), length(AMPLITUDE), length(delta_omega_range));    % Store zero error for each num
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for numa=1:length(AMPLITUDE)
% Iterate over each omega_bar value
   a=AMPLITUDE(numa);
for numw = 1:length(OMEGA_bar)
    omega_bar = OMEGA_bar(numw);
    % Iterate over delta omega range
    for numd = 1:length(delta_omega_range)
        delta_omega = delta_omega_range(numd);
        omega1 = (2 * omega_bar - delta_omega) / 2;
        omega2 = (2 * omega_bar + delta_omega) / 2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x11(1)=0.5; y11(1)=0.5; x21(1)=1; y21(1)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:10000
    %%Pair 1
    x11(i)=x11(i-1)+((a-x11(i-1)^2-y11(i-1)^2)*x11(i-1)-omega1*y11(i-1)+G*C12*(x21(i-1)-x11(i-1)))*dt; %osc1
    y11(i)=y11(i-1)+((a-x11(i-1)^2-y11(i-1)^2)*y11(i-1)+omega1*x11(i-1)+G*C12*(y21(i-1)-y11(i-1)))*dt; %osc1
    
    x21(i)=x21(i-1)+((a-x21(i-1)^2-y21(i-1)^2)*x21(i-1)-omega2*y21(i-1)+G*C21*(x11(i-1)-x21(i-1)))*dt; %osc2
    y21(i)=y21(i-1)+((a-x21(i-1)^2-y21(i-1)^2)*y21(i-1)+omega2*x21(i-1)+G*C21*(y11(i-1)-y21(i-1)))*dt; %osc2 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   x11 = x11(1000:10000); y11 = y11(1000:10000); x21 = x21(1000:10000); y21 = y21(1000:10000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Observation                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N_m=17;
    N_measurements=N_m;
    N_basis=17;
    index=randi([100,499],1,N_m);
    X11dot1=zeros([1,N_measurements]);
    Y11dot1=zeros([1,N_measurements]);
    X21dot2=zeros([1,N_measurements]);
    Y21dot2=zeros([1,N_measurements]);
    
    for ni=1:N_measurements
        X11dot1(ni)=(x11(index(ni)+1)-x11(index(ni)))/dt;
        Y11dot1(ni)=(y11(index(ni)+1)-y11(index(ni)))/dt;
        X21dot2(ni)=(x21(index(ni)+1)-x21(index(ni)))/dt;
        Y21dot2(ni)=(y21(index(ni)+1)-y21(index(ni)))/dt;
    end
    MX1=zeros([N_measurements,N_basis]);
    for i=1:N_measurements
        for j=1:N_basis
            if j==1
               MX1(i,j)=1;
            elseif j==2
               MX1(i,j)=x11(index(i));
            elseif j==3
               MX1(i,j)=y11(index(i));
            elseif j==4
               MX1(i,j)=x11(index(i))*y11(index(i));
            elseif j==5
               MX1(i,j)=x11(index(i))^2;
            elseif j==6
               MX1(i,j)=y11(index(i))^2;
               elseif j==7
                MX1(i,j)=x11(index(i))^2*y11(index(i))^2;
            elseif j==8
                MX1(i,j)=x11(index(i))^2*y11(index(i));
            elseif j==9
                MX1(i,j)=x11(index(i))*y11(index(i))^2;
            elseif j==10
                MX1(i,j)=x11(index(i))^3;
            elseif j==11
                MX1(i,j)=y11(index(i))^3;  
            elseif j==12
                MX1(i,j)=x11(index(i))^3*y11(index(i));  
            elseif j==13
                MX1(i,j)=x11(index(i))*y11(index(i))^3;
            elseif j==14
                MX1(i,j)=x11(index(i))^3*y11(index(i))^2;  
            elseif j==15
                MX1(i,j)=x11(index(i))^2*y11(index(i))^3;  
            elseif j==16
                MX1(i,j)=x11(index(i))^3*y11(index(i))^3;  
            else 
                MX1(i,j)=x21(index(i));
            end    
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       RIP                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Norms=zeros([1,N_basis]);
    for j=1:N_basis
        Norms(j)=norm(MX1(:,j));
    end
    for i=1:N_measurements
        for j=1:N_basis
            MX1(i,j)=MX1(i,j)/Norms(j);
        end
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ksaix11_ref=[0 (a-G*C12) -omega1  0 0 0 0 0 -1 -1 0 0 0 0 0 0 G*C12];%nonzero:(a1-G*C12) -omega1 G*C12 -1 -1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                       Reconstruction                      %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda=0.001;
    X11dot1=X11dot1';
    cd '/home/hadia/cvx'
    cvx_setup
    cvx_begin 
      variable ksaix11(17); 
      minimize(norm(MX1*ksaix11-X11dot1)+ lambda*norm(ksaix11,1)); 
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       RIP                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:N_basis
        ksaix11(j)=ksaix11(j)/Norms(j);
    end 
     ksaix11;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %                Reconstruction Errors                      %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data(count)=N_measurements/N_basis;
    Error_nonzero_x11=0;%initialization
    Error_zero_x11=0;%initialization
    Num_nz=0;%initialization
    Num_z=0;%initialization
    for n=1:17
        if n==2||n==3||n==9||n==10||n==17
             Error_nonzero_x11=Error_nonzero_x11+abs((ksaix11(n)-ksaix11_ref(n))/abs(ksaix11_ref(n)));
             Num_nz=Num_nz+1;
         else
           Error_zero_x11=Error_zero_x11+abs(ksaix11(n));
            Num_z=Num_z+1;
         end  
     end
        Error_nonzero_x11_all(numw,numa,numd) = Error_nonzero_x11/Num_nz;
        Error_zero_x11_all(numw,numa,numd) = Error_zero_x11/Num_z;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
    
    % Save results for each delta_omega
    saveDir = '/home/hadia/cvx/del1';  % Adjust to your desired directory
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    
    for numd = 1:length(delta_omega_range)
        delta_omega = delta_omega_range(numd);
        error_slice_nonzero_x11 = Error_nonzero_x11_all(:,:,numd);
        error_slice_zero_x11 = Error_zero_x11_all(:,:,numd);
        filename = sprintf('Error_slices_delta_omega_%.4f_avg%d.mat', delta_omega, avg);
        
        fullFilePath = fullfile(saveDir, filename);
        save(fullFilePath, 'error_slice_nonzero_x11', 'error_slice_zero_x11');
        
        fprintf('Saved error slices for delta_omega = %.4f avg%d as %s\n', delta_omega, avg, fullFilePath);
    end
end