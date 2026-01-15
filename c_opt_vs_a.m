clear; clc; close all;

% Parameters
R   = 8.314;          % J/(mol·K)
T   = 298.15;         % Temperature
rho = 1000;           % Density of water in kg/m^3
g   = 9.81;           % Acceleration due to gravity in m/s^2
m   = -0.001e6;       % soil water potential
L   = 1.0;            % Phloem (sieve tube) Length
u   = m - rho*g*L;    % leaf water potential
k   = 1.0e-12;        % permeability (fixed)

% Concentration range 
c_min      = 0;
c_max      = 2000;
num_points = 400;
c          = linspace(c_min, c_max, num_points)';   % mol/m^3

% Range of a values (same as your J(c) plot)
a_values = linspace(1e-6, 1e-4, 100);   % 1 µm to 100 µm

% Dynamic viscosity 
nu = Dynamic_Viscosity_T_Sucrose_CW(T, c);
nu = nu(:);

% Storage for c_opt(a)
c_opt_values = zeros(size(a_values));

% Storage for c_opt(a)
J_opt_values = zeros(size(a_values));

fprintf('-------------------------------------------------------------\n');
fprintf('   a (m)             c_opt (mol/m^3)\n');
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

    % Find c_opt
    [~, idx] = max(J);
    c_opt = c(idx);

    c_opt_values(i) = c_opt;
    J_opt_values(i) = max(J);

    % Print numeric results
    fprintf('  %.2e\t\t%8.2f\n', a, c_opt);

    % figure
    % plot(c,J)

end

fprintf('-------------------------------------------------------------\n');

%% Plot c_opt vs a

figure('Color','w');
plot(a_values, c_opt_values, 'o-', 'LineWidth', 1.8, 'MarkerSize', 6);

figure('Color','w');
plot(a_values, J_opt_values, 'o-', 'LineWidth', 1.8, 'MarkerSize', 6);

grid on; box on;
xlabel('Phloem radius, a (m)', 'Color', 'k');
ylabel('Optimal concentration, c_{opt} (mol/m^3)', 'Color', 'k');
title('c_{opt} as a function of phloem radius a', 'Color', 'k');

ax = gca;
ax.XColor      = 'k';
ax.YColor      = 'k';
ax.Title.Color = 'k';
ax.GridColor   = [0.7 0.7 0.7];
ax.GridAlpha   = 0.3;
set(ax,'Layer','top');
