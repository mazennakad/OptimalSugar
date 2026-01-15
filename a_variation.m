clear; clc; close all;

% Parameters
R = 8.314;          % J/(mol·K)
T = 298.15;         % Temperature
rho = 1000;         % Density of water in kg/m^3
g = 9.81;           % Acceleration due to gravity in m/s^2
m = -0.001e6;       % soil water potential
L = 1.0;            % Phloem (sieve tube) Length
u = m - rho*g*L;    % leaf water potential
k = 1.0e-12;        % permeability 

% Concentration range 
c_min = 0;
c_max = 2000;
num_points = 400;
c = linspace(c_min, c_max, num_points)';   % mol/m^3

% Range of a values (radius)
a_values = linspace(1e-6, 1e-4, 10);   % 1 µm to 100 µm

% Dynamic viscosity 
nu = Dynamic_Viscosity_T_Sucrose_CW(T, c);
nu = nu(:);

% Plot setup 
figure('Color','w');
hold on;
grid on; box on;
xlabel('Concentration, c (mol/m^3)', 'Color', 'k');
ylabel('Flow, J(c)', 'Color', 'k');
title('Effect of phloem radius (a) on J vs c at T = 298.15 K', 'Color', 'k');

ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
ax.ZColor = 'k';
ax.Title.Color = 'k';
ax.GridColor = [0.8 0.8 0.8];
ax.GridAlpha  = 0.2;
set(ax,'Layer','top');

colors = lines(length(a_values));

fprintf('-------------------------------------------------------------\n');
fprintf('   a (m)             c_opt (mol/m^3)         J_max\n');
fprintf('-------------------------------------------------------------\n');

% Loop over a values 
for i = 1:length(a_values)
    
    a = a_values(i);

    % J(c)
    N = -m .* k .* L .* a.^2 .* R .* T .* (c.^2) ...
        - m .* u .* a.^2 .* k .* L .* c;

    D = R .* T .* c .* (a.^3 + 4 .* k .* nu .* L.^2) ...
        + m .* (a.^3 - 4 .* k .* nu .* L.^2);

    J = N ./ D;

    % Find maximum J and c_opt
    [J_max, idx] = max(J);
    c_opt = c(idx);

    % Plot J(c)
    plot(c, J, 'Color', colors(i,:), 'LineWidth', 1.6);

    % Mark optimum
    plot(c_opt, J_max, 'o', 'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), 'MarkerSize', 5);

    % Print results
    fprintf('  %.2e\t\t%8.2f\t\t%.4e\n', a, c_opt, J_max);
end

fprintf('-------------------------------------------------------------\n');

% Remove legend
% (none added now)

%% White arrow showing decreasing a-direction
annotation('textarrow', ...
   [0.72 0.58], ...   % x coordinates (normalized figure)
   [0.75 0.55], ...   % y coordinates
   'String', 'Decreasing a  →', ...
   'Color', 'w', ...
   'FontSize', 12, ...
   'LineWidth', 1.5);

