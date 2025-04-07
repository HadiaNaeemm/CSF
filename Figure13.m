% Clear workspace
close all;
clear;
clc;

% Define parameters
G = 1; 
dt = 0.01; % Step size
ALPHA_bar = 0.03;
delta_alpha = -0.05; % Adjust as needed
OMEGA = 6; % Omega value
iterations = 80; % Total number of iterations
theta = linspace(0, 2 * pi, iterations + 1); % Angle values spread over 0 to 2π

% Condition 1
C12_1 = 0.2; 
C21_1 = -0.5;
alpha1_1 = (2 * ALPHA_bar - delta_alpha) / 2;
alpha2_1 = (2 * ALPHA_bar + delta_alpha) / 2;
rho1_0_1 = sqrt(alpha1_1); % ρ₁₀
rho2_0_1 = sqrt(alpha2_1); % ρ₂₀

% Condition 2
C12_2 = -0.5; 
C21_2 = 0.2;
alpha1_2 = (2 * ALPHA_bar - delta_alpha) / 2;
alpha2_2 = (2 * ALPHA_bar + delta_alpha) / 2;
rho1_0_2 = sqrt(alpha1_2); % ρ₁₀
rho2_0_2 = sqrt(alpha2_2); % ρ₂₀
% set(gcf, 'Position', [100, 100, 800, 1000]);

% Create a figure with 2x2 subplots
figure;

% Subplot 1: Linear plot for Condition 1
subplot(3, 2, 1);
rho1_1 = zeros(1, iterations + 1); % Pre-allocate for speed
rho2_1 = zeros(1, iterations + 1); % Pre-allocate for speed
rho1_1(1) = rho1_0_1; % Initial condition
rho2_1(1) = rho2_0_1; % Initial condition
for n = 1:iterations
    rho1_1(n + 1) = sqrt(alpha1_1 + G * C12_1 * (rho2_1(n) / rho1_1(n) - 1));
    rho2_1(n + 1) = sqrt(alpha2_1 + G * C21_1 * (rho1_1(n) / rho2_1(n) - 1));
end
plot(0:iterations, rho1_1, 'b-', 'LineWidth', 2); hold on;
plot(0:iterations, rho2_1, 'r--', 'LineWidth', 2);
xlabel('Iteration Number');
ylabel('$\rho$', 'Interpreter', 'latex');
% title('Condition 1: Convergence of ρ_1 and ρ_2');
% title('Chaotic Behavior');
legend({'$\rho_1$', '$\rho_2$'}, 'Interpreter', 'latex');

% Subplot 2: Polar plot for Condition 1
subplot(3, 2, 2);
polarplot(theta, rho1_1, 'b-', 'LineWidth', 2); hold on;
polarplot(theta, rho2_1, 'r--', 'LineWidth', 2);
% title('Condition 1: Circular Representation of ρ_1 and ρ_2');
% xlabel('x');
% ylabel('y')
legend({'$\rho_1$', '$\rho_2$'}, 'Interpreter', 'latex', 'Location', 'best');
% Access the polar axes and remove grid and ticks
ax = gca; % Get the current axis

% Remove angular labels and grid
ax.ThetaTickLabel = []; % Remove angular tick labels
ax.RLim = [0 1]; % Set radial limits from 0 to 1 (or another range as needed)
ax.RTickLabel = []; % Remove radial ticks

% Turn off the grid lines
ax.ThetaGrid = 'off'; % Remove angular grid
ax.RGrid = 'off'; % Remove radial grid

% Remove the circular border line (black circle)
ax.Box = 'off'; % This removes the box around the polar plot


% Subplot 3: Linear plot for Condition 2
subplot(3, 2, 3);
rho1_2 = zeros(1, iterations + 1); % Pre-allocate for speed
rho2_2 = zeros(1, iterations + 1); % Pre-allocate for speed
rho1_2(1) = rho1_0_2; % Initial condition
rho2_2(1) = rho2_0_2; % Initial condition
for n = 1:iterations
    rho1_2(n + 1) = sqrt(alpha1_2 + G * C12_2 * (rho2_2(n) / rho1_2(n) - 1));
    rho2_2(n + 1) = sqrt(alpha2_2 + G * C21_2 * (rho1_2(n) / rho2_2(n) - 1));
end
plot(0:iterations, rho1_2, 'b-', 'LineWidth', 2); hold on;
plot(0:iterations, rho2_2, 'r--', 'LineWidth', 2);
xlabel('Iteration Number');
ylabel('$\rho$', 'Interpreter', 'latex');
% title('Condition 2: Convergence of ρ_1 and ρ_2');
% title('Stable Behavior');
legend({'$\rho_1$', '$\rho_2$'}, 'Interpreter', 'latex');
% grid on;

% Subplot 4: Polar plot for Condition 2
subplot(3, 2, 4);
polarplot(theta, rho1_2, 'b-', 'LineWidth', 2); hold on;
polarplot(theta, rho2_2, 'r--', 'LineWidth', 2);
% title('Condition 2: Circular Representation of ρ_1 and ρ_2');
legend({'$\rho_1$', '$\rho_2$'}, 'Interpreter', 'latex', 'Location', 'best');
% Access the polar axes and remove grid and ticks
ax = gca; % Get the current axis

% Remove angular labels and grid
ax.ThetaTickLabel = []; % Remove angular tick labels
ax.RLim = [0 1]; % Set radial limits from 0 to 1 (or another range as needed)
ax.RTickLabel = []; % Remove radial ticks

% Turn off the grid lines
ax.ThetaGrid = 'off'; % Remove angular grid
ax.RGrid = 'off'; % Remove radial grid

% Remove the circular border line (black circle)
ax.Box = 'off'; % This removes the box around the polar plot

% Parameters for the first set (initial)
G1 = 1; 
C12_1 = 0.2; 
C21_1 = -0.5;
ALPHA_bar1 = 0.03;
delta_alpha1 = -0.05; % Adjust as needed

% Derived parameters for the first set
alpha1_1 = (2 * ALPHA_bar1 - delta_alpha1) / 2;
alpha2_1 = (2 * ALPHA_bar1 + delta_alpha1) / 2;

% Initial conditions for both sets
rho1_0_1 = sqrt(alpha1_1); % Initial value for ρ_1 for first set
rho2_0_1 = sqrt(alpha2_1); % Initial value for ρ_2 for first set


% Iteration setup
iterations = 100000; % Total number of iterations
transient_cutoff = 1000; % Number of transient iterations to discard

% Pre-allocate for speed
rho1_1 = zeros(1, iterations + 1);
rho2_1 = zeros(1, iterations + 1);

% Initial conditions
rho1_1(1) = rho1_0_1;
rho2_1(1) = rho2_0_1;


% Iterative computation for the first set of parameters
for n = 1:iterations
    rho1_1(n + 1) = sqrt(alpha1_1 + G1 * C12_1 * (rho2_1(n) / rho1_1(n) - 1));
    rho2_1(n + 1) = sqrt(alpha2_1 + G1 * C21_1 * (rho1_1(n) / rho2_1(n) - 1));
end


% Remove transient part for phase diagrams
rho1_1_phase = rho1_1(transient_cutoff:end);
rho2_1_phase = rho2_1(transient_cutoff:end);


% Create the figure with 4 subplots
% figure;

% Phase diagram for ρ_1 (first set)
subplot(3, 2, 5);
plot(rho1_1_phase(1:end-1), rho1_1_phase(2:end), 'b.', 'MarkerSize', 3);
xlabel('$\rho_1(t)$', 'Interpreter', 'latex');
ylabel('$\rho_1(t+1)$', 'Interpreter', 'latex');
title('$\rho_1$', 'Interpreter', 'latex');
% grid on;

% Phase diagram for ρ_2 (first set)
subplot(3, 2, 6);
plot(rho2_1_phase(1:end-1), rho2_1_phase(2:end), 'r.', 'MarkerSize', 3);
xlabel('$\rho_2(t)$', 'Interpreter', 'latex');
ylabel('$\rho_2(t+1)$', 'Interpreter', 'latex');
title('$\rho_2$', 'Interpreter', 'latex');
% grid on;


