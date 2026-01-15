clear; clc; close all;

% Fixed Parameters
R   = 8.314;          % J/(mol·K)
T   = 298.15;         % Temperature
m   = -0.001e6;       % soil water potential
L   = 1.0;            % Phloem (sieve tube) Length
a   = 1.0e-5;         % Phloem (sieve tube) radius
rho = 1000;           % Density of water in kg/m^3
g   = 9.81;           % Acceleration due to gravity in m/s^2

% Leaf water potential 
u = m - rho * g * L;

% Concentration range 
c_min     = 0;
c_max     = 2000;
num_points = 400;
c = linspace(c_min, c_max, num_points)';   % mol/m^3

% Range of k values
k_values = linspace(1.0e-15, 1.0e-12, 10);

% Dynamic viscosity
nu = Dynamic_Viscosity_T_Sucrose_CW(T, c);
nu = nu(:);

% Prepare plot 
figure('Color','w');
hold on;
grid on; box on;
xlabel('Concentration, c (mol/m^3)', 'Color', 'k');
ylabel('Flow, J(c)', 'Color', 'k');
title('J vs c for different k values', 'Color', 'k');

ax = gca;
ax.XColor      = 'k';
ax.YColor      = 'k';
ax.ZColor      = 'k';
ax.Title.Color = 'k';
ax.GridColor   = [0.8 0.8 0.8];
ax.GridAlpha   = 0.2;
set(ax,'Layer','top');

colors = lines(length(k_values));    

fprintf('-----------------------------------------------\n');
fprintf('   k-value         c_opt (mol/m^3)      J_max\n');
fprintf('-----------------------------------------------\n');

% Loop over k values
for i = 1:length(k_values)
    k = k_values(i);

    % Compute J(c)
    N = -m .* k .* L .* a.^2 .* R .* T .* (c.^2) ...
        - m .* u .* a.^2 .* k .* L .* c;
    D = R .* T .* c .* (a.^3 + 4 .* k .* nu .* L.^2) ...
        + m .* (a.^3 - 4 .* k .* nu .* L.^2);
    J = N ./ D;

    % Find maximum J and c_opt
    [J_max, idx] = max(J);
    c_opt = c(idx);

    % Plot J(c) (no legend)
    plot(c, J, 'Color', colors(i,:), 'LineWidth', 1.6);

    % Mark the optimum point
    plot(c_opt, J_max, 'o', 'Color', colors(i,:), ...
         'MarkerFaceColor', colors(i,:), 'MarkerSize', 5);

    % Print results 
    fprintf('  %.1e\t\t%8.2f\t\t%.4e\n', k, c_opt, J_max);
end

fprintf('-----------------------------------------------\n');

% White arrow showing decreasing k direction
annotation('textarrow', ...
    [0.75 0.60], ...   % x: start → end
    [0.80 0.65], ...   % y: start → end
    'String', 'Decreasing k \rightarrow', ...
    'Color', 'w', ...          % WHITE arrow
    'FontSize', 12, ...
    'LineWidth', 1.5);
