%%
clear; clc; close all;

% L values (phloem length)
L_values = logspace(-2, 2, 20);   % 20 values from 0.01 m to 100 m

% Fixed parameters 
R = 8.314;          % J/(mol·K)
T = 298.15;         % Temperature
rho = 1000;         % Density of water in kg/m^3
g = 9.81;           % Acceleration due to gravity in m/s^2
m = -0.001e6;       % soil water potential
a = 1.0e-5;         % Phloem (sieve tube) radius
k = 1.0e-12;        % permeability 

% Concentration range 
c_min = 0;
c_max = 2000;
num_points = 400;
c = linspace(c_min, c_max, num_points)';   % mol/m^3

% Dynamic viscosity (depends only on T and c)
nu = Dynamic_Viscosity_T_Sucrose_CW(T, c);
nu = nu(:);   % ensure column vector

% Plot setup 
figure('Color','w');
hold on; grid on; box on;
xlabel('Concentration, c (mol/m^3)');
ylabel('Flow, J(c)');
title('J vs c for different phloem lengths L (T = 298.15 K, k,a fixed)');
grid on; box on;
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
ax.ZColor = 'k';
ax.Title.Color = 'k';

ax.GridColor = [0.8 0.8 0.8];
ax.GridAlpha = 0.2;
set(ax,'Layer','top');
colors = lines(length(L_values));

fprintf('-----------------------------------------------\n');
fprintf('   L (m)          c_opt (mol/m^3)       J_max\n');
fprintf('-----------------------------------------------\n');

% Loop over L values
for i = 1:length(L_values)
    L = L_values(i);

    % leaf water potential
    u = m - rho*g*L;

    % Compute J(c)
    N = -m .* k .* L .* a.^2 .* R .* T .* (c.^2) - m .* u .* a.^2 .* k .* L .* c;
    D = R .* T .* c .* (a.^3 + 4 .* k .* nu .* L.^2) + m .* (a.^3 - 4 .* k .* nu .* L.^2);
    J = N ./ D;

    % Find maximum point
    [J_max, idx] = max(J);
    c_opt = c(idx);

    % Plot curve (no legend)
    plot(c, J, 'Color', colors(i,:), 'LineWidth', 1.6);

    % Mark maximum
    plot(c_opt, J_max, 'o', 'Color', colors(i,:), ...
         'MarkerFaceColor', colors(i,:), 'MarkerSize', 5);

    % Print values
    fprintf('  %6.2f\t\t%10.2f\t\t%.4e\n', L, c_opt, J_max);
end

fprintf('-----------------------------------------------\n');

% White arrow showing decreasing L direction
annotation('textarrow', ...
    [0.75 0.6], [0.8 0.65], ...   % arrow position on figure
    'String', 'Decreasing L \rightarrow', ...
    'Color', 'w', 'FontSize', 12, 'LineWidth', 1.5);
