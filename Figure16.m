close all;
clear;
clc;

% Parameters for Code 1
G1 = 1; 
C12_1 = 0.2; 
C21_1 = -0.5;
ALPHA_bar_1 = 0.03; % Fixed base value for alpha_bar
delta_alpha_1 = -0.05; % Adjust as needed

% Parameters for Code 2
G2 = 1; 
C12_2 = -0.5; 
C21_2 = 0.2;
ALPHA_bar_2 = 0.03; % Fixed base value for alpha_bar
delta_alpha_2 = -0.05; % Adjust as needed

% Common parameters
alpha1_values = linspace(0.01, 0.2, 600); % Range of alpha1 to sweep
alpha2_values = linspace(0.01, 0.2, 600); % Range of alpha2 to sweep
iterations = 10000; % Total number of iterations for each alpha1
transient_cutoff = 5000; % Ignore initial transient part

% Initialize storage for steady-state values
rho1_steady_1 = [];
rho2_steady_1 = [];
rho1_steady_2 = [];
rho2_steady_2 = [];

% Code 1: Bifurcation Diagram for rho1 vs alpha1
alpha2_1 = (2 * ALPHA_bar_1 + delta_alpha_1) / 2; % Fixed alpha2
for alpha1 = alpha1_values
    % Initial conditions
    rho1 = sqrt(alpha1); % \rho_1^(0)
    rho2 = sqrt(alpha2_1); % \rho_2^(0)
    
    % Iterative computation
    rho1_series = zeros(1, iterations); % Pre-allocate
    for n = 1:iterations
        rho1_next = sqrt(alpha1 + G1 * C12_1 * (rho2 / rho1 - 1));
        rho2_next = sqrt(alpha2_1 + G1 * C21_1 * (rho1 / rho2 - 1));
        rho1 = rho1_next;
        rho2 = rho2_next;
        rho1_series(n) = rho1; % Store rho1 values
    end
    
    % Remove transient part and store the steady-state values
    rho1_steady_1 = [rho1_steady_1; alpha1 * ones(length(rho1_series(transient_cutoff:end)), 1), ...
                     rho1_series(transient_cutoff:end)'];
end

% Code 1: Bifurcation Diagram for rho2 vs alpha2
alpha1_1 = (2 * ALPHA_bar_1 - delta_alpha_1) / 2; % Fixed alpha1
for alpha2 = alpha2_values
    % Initial conditions
    rho1 = sqrt(alpha1_1); % \rho_1^(0)
    rho2 = sqrt(alpha2); % \rho_2^(0)
    
    % Iterative computation
    rho2_series = zeros(1, iterations); % Pre-allocate
    for n = 1:iterations
        rho1_next = sqrt(alpha1_1 + G1 * C12_1 * (rho2 / rho1 - 1));
        rho2_next = sqrt(alpha2 + G1 * C21_1 * (rho1 / rho2 - 1));
        rho1 = rho1_next;
        rho2 = rho2_next;
        rho2_series(n) = rho2; % Store rho2 values
    end
    
    % Remove transient part and store the steady-state values
    rho2_steady_1 = [rho2_steady_1; alpha2 * ones(length(rho2_series(transient_cutoff:end)), 1), ...
                     rho2_series(transient_cutoff:end)'];
end

% Code 2: Bifurcation Diagram for rho1 vs alpha1
alpha2_2 = (2 * ALPHA_bar_2 + delta_alpha_2) / 2; % Fixed alpha2
for alpha1 = alpha1_values
    % Initial conditions
    rho1 = sqrt(alpha1); % \rho_1^(0)
    rho2 = sqrt(alpha2_2); % \rho_2^(0)
    
    % Iterative computation
    rho1_series = zeros(1, iterations); % Pre-allocate
    for n = 1:iterations
        rho1_next = sqrt(alpha1 + G2 * C12_2 * (rho2 / rho1 - 1));
        rho2_next = sqrt(alpha2_2 + G2 * C21_2 * (rho1 / rho2 - 1));
        rho1 = rho1_next;
        rho2 = rho2_next;
        rho1_series(n) = rho1; % Store rho1 values
    end
    
    % Remove transient part and store the steady-state values
    rho1_steady_2 = [rho1_steady_2; alpha1 * ones(length(rho1_series(transient_cutoff:end)), 1), ...
                     rho1_series(transient_cutoff:end)'];
end

% Code 2: Bifurcation Diagram for rho2 vs alpha2
alpha1_2 = (2 * ALPHA_bar_2 - delta_alpha_2) / 2; % Fixed alpha1
for alpha2 = alpha2_values
    % Initial conditions
    rho1 = sqrt(alpha1_2); % \rho_1^(0)
    rho2 = sqrt(alpha2); % \rho_2^(0)
    
    % Iterative computation
    rho2_series = zeros(1, iterations); % Pre-allocate
    for n = 1:iterations
        rho1_next = sqrt(alpha1_2 + G2 * C12_2 * (rho2 / rho1 - 1));
        rho2_next = sqrt(alpha2 + G2 * C21_2 * (rho1 / rho2 - 1));
        rho1 = rho1_next;
        rho2 = rho2_next;
        rho2_series(n) = rho2; % Store rho2 values
    end
    
    % Remove transient part and store the steady-state values
    rho2_steady_2 = [rho2_steady_2; alpha2 * ones(length(rho2_series(transient_cutoff:end)), 1), ...
                     rho2_series(transient_cutoff:end)'];
end

% Plotting all diagrams as subplots
figure;

% Subplot 1: Bifurcation Diagram for rho1 vs alpha1 (Code 1)
subplot(2, 2, 1);
plot(rho1_steady_1(:, 1), rho1_steady_1(:, 2), 'k.', 'MarkerSize', 2);
xlabel('$\alpha_1$', 'Interpreter', 'latex');
ylabel('$\rho_1$', 'Interpreter', 'latex');
title('Bifurcation Diagram: $\rho_1$ vs $\alpha_1$ (Code 1)', 'Interpreter', 'latex');
grid on;

% Subplot 2: Bifurcation Diagram for rho2 vs alpha2 (Code 1)
subplot(2, 2, 2);
plot(rho2_steady_1(:, 1), rho2_steady_1(:, 2), 'k.', 'MarkerSize', 2);
xlabel('$\alpha_2$', 'Interpreter', 'latex');
ylabel('$\rho_2$', 'Interpreter', 'latex');
title('Bifurcation Diagram: $\rho_2$ vs $\alpha_2$ (Code 1)', 'Interpreter', 'latex');
grid on;

% Subplot 3: Bifurcation Diagram for rho1 vs alpha1 (Code 2)
subplot(2, 2, 3);
plot(rho1_steady_2(:, 1), rho1_steady_2(:, 2), 'k.', 'MarkerSize', 2);
xlabel('$\alpha_1$', 'Interpreter', 'latex');
ylabel('$\rho_1$', 'Interpreter', 'latex');
title('Bifurcation Diagram: $\rho_1$ vs $\alpha_1$ (Code 2)', 'Interpreter', 'latex');
grid on;

% Subplot 4: Bifurcation Diagram for rho2 vs alpha2 (Code 2)
subplot(2, 2, 4);
plot(rho2_steady_2(:, 1), rho2_steady_2(:, 2), 'k.', 'MarkerSize', 2);
xlabel('$\alpha_2$', 'Interpreter', 'latex');
ylabel('$\rho_2$', 'Interpreter', 'latex');
title('Bifurcation Diagram: $\rho_2$ vs $\alpha_2$ (Code 2)', 'Interpreter', 'latex');
grid on;