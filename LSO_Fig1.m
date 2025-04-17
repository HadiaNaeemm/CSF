close all; clear all; clc;

% Parameters
a_values = [-0.09, 0, 0.009]; % Values of 'a' for each regime
omega = 1; % Angular frequency of oscillation
time = 0:0.01:300; % Time vector for trajectory plots
time_phase = [0 100]; % Time range for phase portraits
[x, y] = meshgrid(-2:0.2:2, -2:0.2:2); % Grid for phase portraits

% Initialize arrays to store results
x_results = cell(1,3);
y_results = cell(1,3);

% Colors for trajectories and phase portraits
trajectory_colors = {[1, 0.4, 0.6], [0, 0.2, 0.9], [0, 0.7, 0.3]};  % Magenta, Royal Blue, Emerald Green
phase_colors = {[0.6350, 0.0780, 0.1840], [0.6 0.3 0.9], [1, 0.5, 0.2]};        % Deep Blue, Violet, Burnt Orange
arrow_color = [0.2 0.2 0.2];   % Dark gray for vector field

% Create figure
 % figure('Position', [100, 100, 1000, 600]); % Large canvas for clarity

% --- Trajectory Plots (Top Row) ---
for i = 1:3
    a = a_values(i);
    x0 = 0.1;
    y0 = 0.1;
    z0 = x0 + 1i*y0;
    
    % Stuart-Landau Equation
    f = @(t, z) (a + 1i * omega) * z - abs(z)^2 * z;
    [t, z] = ode45(@(t, z) f(t, z), time, z0);
    
    % Extract real and imaginary parts
    x_results{i} = real(z);
    y_results{i} = imag(z);
    
    % Create axes (no subplot warning)
    ax1 = axes('Position', [0.1 + (i-1)*0.3, 0.55, 0.25, 0.35]);
    plot(x_results{i}, y_results{i}, 'Color', trajectory_colors{i}, 'LineWidth', 2.5);
    xlabel('Re(z)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Im(z)', 'FontSize', 10, 'FontWeight', 'bold');
    axis equal;
    set(ax1, 'LineWidth', 1.2);
end

