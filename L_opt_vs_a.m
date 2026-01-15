%%
clear; clc; close all;

% Fixed parameters 
R   = 8.314;          % J/(mol·K)
T   = 298.15;         % K
rho = 1000;           % kg/m^3 (density of water)
g   = 9.81;           % m/s^2 (gravity)
m   = -0.001e6;       % Pa (soil water potential)
k   = 1.0e-12;        % permeability (fixed)
                      
% Concentration range 
c_min      = 0;
c_max      = 2000;
num_points = 400;
c          = linspace(c_min, c_max, num_points)';   % mol/m^3

% Dynamic viscosity 
nu = Dynamic_Viscosity_T_Sucrose_CW(T, c);
nu = nu(:);

% Radius values (a) for which we compute L_opt
a_values = linspace(1e-6, 1e-4, 20);   % you can change number of points

% Range of L values used to search for L_opt at each a
L_values = logspace(-2, 2, 100);       % from 0.01 m to 100 m

% Storage for L_opt and corresponding J_max
L_opt_values   = zeros(size(a_values));
J_opt_values   = zeros(size(a_values));

fprintf('-------------------------------------------------------------\n');
fprintf('   a (m)             L_opt (m)             J_max(L_opt)\n');
fprintf('-------------------------------------------------------------\n');

% Loop over radius values a
for ia = 1:length(a_values)
    
    a = a_values(ia);
    
    % For this a, we will find J_max as a function of L
    J_max_over_L = zeros(size(L_values));
    
    for iL = 1:length(L_values)
        
        L = L_values(iL);
        
        % Leaf water potential depends on L
        u = m - rho*g*L;    % Pa
        
        % J(c)
        N = -m .* k .* L .* a.^2 .* R .* T .* (c.^2) ...
            - m .* u .* a.^2 .* k .* L .* c;
        
        D = R .* T .* c .* (a.^3 + 4 .* k .* nu .* L.^2) ...
            + m .* (a.^3 - 4 .* k .* nu .* L.^2);
        
        J = N ./ D;
        
        % Max J for this (a, L)
        J_max_over_L(iL) = max(J);
        
    end
    
    % For this a, pick the L that gives the largest J_max
    [J_opt, idx_best_L] = max(J_max_over_L);
    L_opt = L_values(idx_best_L);
    
    L_opt_values(ia) = L_opt;
    J_opt_values(ia) = J_opt;
    
    fprintf('  %.3e\t\t%10.3f\t\t% .4e\n', a, L_opt, J_opt);
end

fprintf('-------------------------------------------------------------\n');

%% Plot L_opt as a function of a

figure('Color','w');
plot(a_values, L_opt_values, 'o-', 'LineWidth', 1.8, 'MarkerSize', 6);
grid on; box on;

xlabel('Phloem radius, a (m)', 'Color', 'k');
ylabel('Optimal length, L_{opt} (m)', 'Color', 'k');
title('L_{opt} as a function of phloem radius a', 'Color', 'k');

ax = gca;
ax.XColor      = 'k';
ax.YColor      = 'k';
ax.Title.Color = 'k';
ax.GridColor   = [0.7 0.7 0.7];
ax.GridAlpha   = 0.3;
set(ax,'Layer','top');
