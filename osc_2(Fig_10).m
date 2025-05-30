close all;
clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
avgnum=100;
%value of constants
G=1;C12=-0.5;C21=0.2;
dt=0.01;   %step size
ALPHA_bar = 0.03:0.01:0.09;
delta_alpha_range = -0.05:0.002:0.05;
OMEGA=1:1:10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Error_nonzero_x21_avg = zeros(length(ALPHA_bar), length(OMEGA), length(delta_alpha_range)); 
Error_zero_x21_avg = zeros(length(ALPHA_bar), length(OMEGA), length(delta_alpha_range));    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for avg=1:avgnum
    Error_nonzero_x21_all = zeros(length(ALPHA_bar), length(OMEGA), length(delta_alpha_range)); 
    Error_zero_x21_all = zeros(length(ALPHA_bar), length(OMEGA), length(delta_alpha_range));    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for numw=1:length(OMEGA)
    w=OMEGA(numw);
for numa = 1:length(ALPHA_bar)
    alpha_bar = ALPHA_bar(numa);
    for numd = 1:length(delta_alpha_range)
        delta_alpha = delta_alpha_range(numd);
        a1 = (2 * alpha_bar - delta_alpha) / 2;
        a2 = (2 * alpha_bar + delta_alpha) / 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x11(1)=0.5; y11(1)=0.5; x21(1)=1; y21(1)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:10000
    %%Pair 1
    x11(i)=x11(i-1)+((a1-x11(i-1)^2-y11(i-1)^2)*x11(i-1)-w*y11(i-1)+G*C12*(x21(i-1)-x11(i-1)))*dt; %osc1
    y11(i)=y11(i-1)+((a1-x11(i-1)^2-y11(i-1)^2)*y11(i-1)+w*x11(i-1)+G*C12*(y21(i-1)-y11(i-1)))*dt; %osc1
    
    x21(i)=x21(i-1)+((a2-x21(i-1)^2-y21(i-1)^2)*x21(i-1)-w*y21(i-1)+G*C21*(x11(i-1)-x21(i-1)))*dt; %osc2
    y21(i)=y21(i-1)+((a2-x21(i-1)^2-y21(i-1)^2)*y21(i-1)+w*x21(i-1)+G*C21*(y11(i-1)-y21(i-1)))*dt; %osc2 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  x11 = x11(1000:10000); y11 = y11(1000:10000); x21 = x21(1000:10000); y21 = y21(1000:10000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N_m=100;
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
    MX2=zeros([N_measurements,N_basis]);
    for i=1:N_measurements
        for j=1:N_basis
            if j==1
               MX2(i,j)=1;
            elseif j==2
               MX2(i,j)=x21(index(i));
            elseif j==3
               MX2(i,j)=y21(index(i));
            elseif j==4
               MX2(i,j)=x21(index(i))*y21(index(i));
            elseif j==5
               MX2(i,j)=x21(index(i))^2;
            elseif j==6
               MX2(i,j)=y21(index(i))^2;
               elseif j==7
                MX2(i,j)=x21(index(i))^2*y21(index(i))^2;
            elseif j==8
                MX2(i,j)=x21(index(i))^2*y21(index(i));
            elseif j==9
                MX2(i,j)=x21(index(i))*y21(index(i))^2;
            elseif j==10
                MX2(i,j)=x21(index(i))^3;
            elseif j==11
                MX2(i,j)=y21(index(i))^3;  
            elseif j==12
                MX2(i,j)=x21(index(i))^3*y21(index(i));  
            elseif j==13
                MX2(i,j)=x21(index(i))*y21(index(i))^3;
            elseif j==14
                MX2(i,j)=x21(index(i))^3*y21(index(i))^2;  
            elseif j==15
                MX2(i,j)=x21(index(i))^2*y21(index(i))^3;  
            elseif j==16
                MX2(i,j)=x21(index(i))^3*y21(index(i))^3;  
            else 
                MX2(i,j)=x11(index(i));
            end    
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       RIP                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Norms=zeros([1,N_basis]);
    for j=1:N_basis
        Norms(j)=norm(MX2(:,j));
    end
    for i=1:N_measurements
        for j=1:N_basis
            MX2(i,j)=MX2(i,j)/Norms(j);
        end
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ksaix21_ref=[0 (a2-G*C21) -w  0 0 0 0 0 -1 -1 0 0 0 0 0 0 G*C21];%nonzero:(a2-G*C21) -w G*C21 -1 -1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                       Reconstruction                      %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  lambda=0.001;
  X21dot2=X21dot2';
  cd '/home/hadia/cvx'
  cvx_setup
  cvx_begin 
  variable ksaix21(17); 
  minimize(norm(MX2*ksaix21-X21dot2)+ lambda*norm(ksaix21,1)); 
  cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       RIP                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:N_basis
        ksaix21(j)=ksaix21(j)/Norms(j);
    end 
     ksaix21;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %                Reconstruction Errors                      %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data(count)=N_measurements/N_basis;
    Error_nonzero_x21=0;%initialization
    Error_zero_x21=0;%initialization
    Num_nz=0;%initialization
    Num_z=0;%initialization
    for n=1:17
        if n==2||n==3||n==9||n==10||n==17
             Error_nonzero_x21=Error_nonzero_x21+abs((ksaix21(n)-ksaix21_ref(n))/abs(ksaix21_ref(n)));
             Num_nz=Num_nz+1;
         else
           Error_zero_x21=Error_zero_x21+abs(ksaix21(n));
            Num_z=Num_z+1;
         end  
     end
        Error_nonzero_x21_all(numa,numw,numd) = Error_nonzero_x21/Num_nz;
        Error_zero_x21_all(numa,numw,numd) = Error_zero_x21/Num_z;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end 
end 
 saveDir = '/home/hadia/cvx/delamp11';
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    for numd = 1:length(delta_alpha_range)
        delta_alpha = delta_alpha_range(numd);  % Fix variable name
        error_slice_nonzero_x21 = Error_nonzero_x21_all(:,:,numd);
        error_slice_zero_x21 = Error_zero_x21_all(:,:,numd);

        filename = sprintf('Error_slices_delta_alpha_%d_avg%d.mat', round(delta_alpha * 10000), avg); % Fix filename
        fullFilePath = fullfile(saveDir, filename);

        save(fullFilePath, 'error_slice_nonzero_x21', 'error_slice_zero_x21');
        fprintf('Saved error slices for delta_alpha = %.4f avg%d as %s\n', delta_alpha, avg, fullFilePath);
    end
end