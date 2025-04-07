% MATLAB script to plot bifurcation diagrams for rho1 and rho2 in the same panel

% Parameters
G = 1; 
C12_1 = 0.2; 
C21_1 = -0.5;
C12_2 = -0.5; 
C21_2 = 0.2;
dt = 0.01; % Step size
ALPHA_bar = 0.03;

% Derived parameters
alpha_bar = ALPHA_bar;

% Initial conditions
rho1_0 = sqrt(alpha_bar); % ρ_1^(0)
rho2_0 = sqrt(alpha_bar); % ρ_2^(0)

% Iteration setup
iterations = 5000; % Total number of iterations
initial_iterations = 1000; % Initial iterations to discard transient behavior

% Delta Alpha range
delta_alpha_range = -0.05:0.000001:0.05;

% Function to compute bifurcation values
compute_bifurcation = @(C12, C21) arrayfun(@(delta_alpha) ...
    mean(iterate_rho(delta_alpha, C12, C21, alpha_bar, G, iterations, initial_iterations, rho1_0, rho2_0)), ...
    delta_alpha_range);

% Compute bifurcation diagrams for both cases
bifurcation_rho1 = compute_bifurcation(C12_1, C21_1);
bifurcation_rho2 = compute_bifurcation(C12_2, C21_2);

% Plot both on the same figure
% figure;
% hold on;
% plot(delta_alpha_range, bifurcation_rho1, 'b.', 'MarkerSize', 2);
% plot(delta_alpha_range, bifurcation_rho2, 'r.', 'MarkerSize', 2);
% xlabel('Delta Alpha');
% ylabel('Rho');
% title('Bifurcation Diagram for \rho_1 and \rho_2');
% legend({'$\rho_1$', '$\rho_2$'}, 'Interpreter', 'latex');
% grid on;
% hold off;

% Plotting the bifurcation diagrams on the same plot
figure;
plot(delta_alpha_range, bifurcation_rho1, 'b.', 'MarkerSize', 3); % Bifurcation for rho1
hold on;
plot(delta_alpha_range, bifurcation_rho2, 'r.', 'MarkerSize', 3); % Bifurcation for rho2
xlabel('$\Delta_{\alpha}$', 'Interpreter', 'latex');
ylabel('$\rho$', 'Interpreter', 'latex');
title('Combined Bifurcation Diagrams for \rho_1 and \rho_2');
legend({'$\rho_1$', '$\rho_2$'}, 'Interpreter', 'latex');
grid on;
% Adjust vertical range (ylim) to make bifurcations more visible
ylim([0.2 0.5]); % You can adjust this range based on your data

% Function to iterate rho1 and rho2
function rho_values = iterate_rho(delta_alpha, C12, C21, alpha_bar, G, iterations, initial_iterations, rho1_0, rho2_0)
    alpha1 = (2 * alpha_bar - delta_alpha) / 2;
    alpha2 = (2 * alpha_bar + delta_alpha) / 2;
    rho1 = zeros(1, iterations + 1);
    rho2 = zeros(1, iterations + 1);
    rho1(1) = rho1_0;
    rho2(1) = rho2_0;
    for n = 1:iterations
        rho1(n + 1) = sqrt(alpha1 + G * C12 * (rho2(n) / rho1(n) - 1));
        rho2(n + 1) = sqrt(alpha2 + G * C21 * (rho1(n) / rho2(n) - 1));
    end
    rho_values = mean(rho1(initial_iterations+1:end));
end
