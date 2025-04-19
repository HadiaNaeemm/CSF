% Clear workspace
close all;
clear;
clc;

% MATLAB script to plot bifurcation diagrams for rho1 and rho2

% Parameters
G = 1; 
C12_1 = 0.2;   % For rho1 plot
C21_1 = -0.5;

C12_2 = -0.5;  % For rho2 plot
C21_2 = 0.2;

dt = 0.01;
ALPHA_bar = 0.03;
alpha_bar = ALPHA_bar;

% Initial conditions
rho1_0 = sqrt(alpha_bar);
rho2_0 = sqrt(alpha_bar);

% Iteration setup
iterations = 5000;
initial_iterations = 1000;

% Delta alpha range
delta_alpha_range =  -0.05:0.000001:0.05; % Slightly coarser for faster plot

% Preallocate
bifurcation_rho1 = zeros(size(delta_alpha_range));
bifurcation_rho2 = zeros(size(delta_alpha_range));

% Loop for rho1 vs delta alpha (using first coupling set)
for i = 1:length(delta_alpha_range)
    delta_alpha = delta_alpha_range(i);
    bifurcation_rho1(i) = iterate_rho_get(delta_alpha, C12_1, C21_1, alpha_bar, G, iterations, initial_iterations, rho1_0, rho2_0, 1);
end

% Loop for rho2 vs delta alpha (using second coupling set)
for i = 1:length(delta_alpha_range)
    delta_alpha = delta_alpha_range(i);
    bifurcation_rho2(i) = iterate_rho_get(delta_alpha, C12_2, C21_2, alpha_bar, G, iterations, initial_iterations, rho1_0, rho2_0, 2);
end

% Plotting
figure;
plot(delta_alpha_range, bifurcation_rho1, 'b.', 'MarkerSize', 5); hold on;
plot(delta_alpha_range, bifurcation_rho2, 'r.', 'MarkerSize', 5);
xlabel('$\Delta_\alpha$', 'Interpreter', 'latex');
ylabel('$\rho$', 'Interpreter', 'latex');
legend({'$\rho_1$ for $C_{12}=0.2$, $C_{21}=-0.5$', ...
        '$\rho_2$ for $C_{12}=-0.5$, $C_{21}=0.2$'}, ...
        'Interpreter', 'latex', 'Location', 'best');
title('Bifurcation Diagram of $\rho_1$ and $\rho_2$', 'Interpreter', 'latex');
% grid on;
ylim([0.1 0.4]);  % Adjust based on your values

% Function definition
function rho_out = iterate_rho_get(delta_alpha, C12, C21, alpha_bar, G, iterations, initial_iterations, rho1_0, rho2_0, target)
    alpha1 = (2 * alpha_bar - delta_alpha) / 2;
    alpha2 = (2 * alpha_bar + delta_alpha) / 2;
    rho1 = zeros(1, iterations + 1);
    rho2 = zeros(1, iterations + 1);
    rho1(1) = rho1_0;
    rho2(1) = rho2_0;
    for n = 1:iterations
        rho1(n+1) = sqrt(alpha1 + G * C12 * (rho2(n)/rho1(n) - 1));
        rho2(n+1) = sqrt(alpha2 + G * C21 * (rho1(n)/rho2(n) - 1));
    end
    if target == 1
        rho_out = mean(rho1(initial_iterations+1:end));
    else
        rho_out = mean(rho2(initial_iterations+1:end));
    end
end
