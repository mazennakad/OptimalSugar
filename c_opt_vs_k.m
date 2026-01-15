% clear; clc; close all;

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

% Range of k values (same as your original code)
k_values = linspace(1.0e-15, 1.0e-12, 100);

% Dynamic viscosity
nu = Dynamic_Viscosity_T_Sucrose_CW(T, c);
nu = nu(:);

% Storage for c_opt results
c_opt_values = zeros(size(k_values));

% Storage for L_opt(L)
J_opt_values = zeros(size(k_values));

fprintf('-----------------------------------------------\n');
fprintf('   k-value         c_opt (mol/m^3)\n');
fprintf('-----------------------------------------------\n');

% Loop over k values to compute c_opt
for i = 1:length(k_values)
    
    k = k_values(i);

    % Compute J(c)
    N = -m .* k .* L .* a.^2 .* R .* T .* (c.^2) ...
        - m .* u .* a.^2 .* k .* L .* c;

    D = R .* T .* c .* (a.^3 + 4 .* k .* nu .* L.^2) ...
        + m .* (a.^3 - 4 .* k .* nu .* L.^2);

    J = N ./ D;

    % Find maximum J and corresponding c_opt
    [~, idx] = max(J);
    c_opt = c(idx);

    c_opt_values(i) = c_opt;
    J_opt_values(i) = max(J);

    fprintf('  %.1e\t\t%8.2f\n', k, c_opt);

end

fprintf('-----------------------------------------------\n');


%% Plot c_opt vs k (normal plot)

figure('Color','w');
plot(k_values, c_opt_values, 'o-', 'LineWidth', 1.8, 'MarkerSize', 6);

figure('Color','w');
plot(k_values, J_opt_values, 'o-', 'LineWidth', 1.8, 'MarkerSize', 6);

grid on; box on;
xlabel('Permeability, k (m^2)', 'Color', 'k');
ylabel('Optimal concentration, c_{opt} (mol/m^3)', 'Color', 'k');
title('c_{opt} as a function of k', 'Color', 'k');

ax = gca;
ax.XColor      = 'k';
ax.YColor      = 'k';
ax.Title.Color = 'k';
ax.GridColor   = [0.7 0.7 0.7];
ax.GridAlpha   = 0.3;
set(ax,'Layer','top');
