% Clear workspace
close all;
clear;
clc;

% Define parameters
G = 1; 
dt = 0.01; % Step size
ALPHA_bar = 0.03;
delta_alpha = 0; % Adjust as needed
OMEGA = 0; % Omega value
iterations = 100; % Total number of iterations
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

% Create a figure with 2x2 subplots
figure;

% Subplot 1: Linear plot for Condition 1
subplot(2, 2, 1);
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
legend({'$\rho_1$', '$\rho_2$'}, 'Interpreter', 'latex');

% Subplot 2: Linear plot for Condition 2
subplot(2, 2, 3);
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
legend({'$\rho_1$', '$\rho_2$'}, 'Interpreter', 'latex');

% Subplot 3: Polar plot for Condition 1
subplot(2, 2, 2);
rho1_1_scaled = rho1_1 * 2;  % Scale the data to enlarge it
rho2_1_scaled = rho2_1 * 2;  % Scale the data to enlarge it

polarplot(theta, rho1_1_scaled, 'b-', 'LineWidth', 2); hold on;
polarplot(theta, rho2_1_scaled, 'r--', 'LineWidth', 2);

legend({'$\rho_1$', '$\rho_2$'}, 'Interpreter', 'latex', 'Location', 'best');
ax = gca;
ax.ThetaTickLabel = []; % Remove angular tick labels
ax.RLim = [0 max([rho1_1_scaled, rho2_1_scaled])]; % Set radial limits to the maximum scaled value
ax.RTickLabel = []; % Remove radial ticks
ax.ThetaGrid = 'off'; % Remove angular grid
ax.RGrid = 'off'; % Remove radial grid
ax.Box = 'off'; % Remove the circular border

% Subplot 4: Polar plot for Condition 2
subplot(2, 2, 4);
rho1_2_scaled = rho1_2 * 2;  % Scale the data to enlarge it
rho2_2_scaled = rho2_2 * 2;  % Scale the data to enlarge it

polarplot(theta, rho1_2_scaled, 'b-', 'LineWidth', 2); hold on;
polarplot(theta, rho2_2_scaled, 'r--', 'LineWidth', 2);

legend({'$\rho_1$', '$\rho_2$'}, 'Interpreter', 'latex', 'Location', 'best');
ax = gca;
ax.ThetaTickLabel = []; % Remove angular tick labels
ax.RLim = [0 max([rho1_2_scaled, rho2_2_scaled])]; % Set radial limits to the maximum scaled value
ax.RTickLabel = []; % Remove radial ticks
ax.ThetaGrid = 'off'; % Remove angular grid
ax.RGrid = 'off'; % Remove radial grid
ax.Box = 'off'; % Remove the circular border
