%%
clear; clc; close all;

% L values (phloem length)
L_values = logspace(-2, 2, 100);   % 20 values from 0.01 m to 100 m

% Fixed parameters 
R   = 8.314;          % J/(mol·K)
T   = 298.15;         % Temperature
rho = 1000;           % Density of water in kg/m^3
g   = 9.81;           % Acceleration due to gravity in m/s^2
m   = -0.001e6;       % soil water potential
a   = 1.0e-5;         % Phloem radius
k   = 1.0e-12;        % permeability

% Concentration range
c_min = 0;
c_max = 2000;
num_points = 400;
c = linspace(c_min, c_max, num_points)';   % mol/m^3

% Dynamic viscosity
nu = Dynamic_Viscosity_T_Sucrose_CW(T, c);
nu = nu(:);

% Storage for c_opt(L)
c_opt_values = zeros(size(L_values));

fprintf('-----------------------------------------------\n');
fprintf('   L (m)              c_opt (mol/m^3)\n');
fprintf('-----------------------------------------------\n');

% Compute c_opt for each L
for i = 1:length(L_values)
    
    L = L_values(i);

    % Leaf water potential
    u = m - rho*g*L;

    % Compute J(c)
    N = -m .* k .* L .* a.^2 .* R .* T .* (c.^2) ...
        - m .* u .* a.^2 .* k .* L .* c;

    D = R .* T .* c .* (a.^3 + 4 .* k .* nu .* L.^2) ...
        + m .* (a.^3 - 4 .* k .* nu .* L.^2);

    J = N ./ D;

    % Find optimal concentration c_opt
    [~, idx] = max(J);
    c_opt = c(idx);

    c_opt_values(i) = c_opt;

    fprintf('   %.3f\t\t%10.2f\n', L, c_opt);

end

fprintf('-----------------------------------------------\n');

%% Plot c_opt vs L
figure('Color','w');
plot(L_values, c_opt_values, 'o-', 'LineWidth', 1.8, 'MarkerSize', 6);

grid on; box on;
xlabel('Phloem length, L (m)', 'Color', 'k');
ylabel('Optimal concentration, c_{opt} (mol/m^3)', 'Color', 'k');
title('c_{opt} as a function of phloem length L', 'Color', 'k');

ax = gca;
ax.XColor      = 'k';
ax.YColor      = 'k';
ax.Title.Color = 'k';
ax.GridColor   = [0.7 0.7 0.7];
ax.GridAlpha   = 0.3;
set(ax,'Layer','top');
