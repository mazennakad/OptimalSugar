%%
clear; clc; close all;

%% Ranges for L and a

% L values (phloem length) used to search for L_threshold
L_values = logspace(-2, 2, 100);      % 100 values from 0.01 m to 100 m

% Radius values (a) for which we compute L_threshold
a_values = linspace(1e-6, 1e-4, 50);  % from 1 µm to 100 µm

%% Fixed parameters 
R   = 8.314;          % J/(mol·K)
T   = 298.15;         % Temperature
rho = 1000;           % Density of water in kg/m^3
g   = 9.81;           % Acceleration due to gravity in m/s^2
m   = -0.001e6;       % Soil water potential
k   = 1.0e-12;        % Permeability (fixed here)

%% Concentration range 
c_min      = 0;
c_max      = 2000;
num_points = 400;
c          = linspace(c_min, c_max, num_points)';   % mol/m^3 (column vector)

% Dynamic viscosity (depends only on T and c)
nu = Dynamic_Viscosity_T_Sucrose_CW(T, c);
nu = nu(:);    % ensure column vector

% Concentration threshold that defines L_threshold
c_threshold = 1500;    % mol/m^3

%% Storage for L_threshold(a)
L_threshold = NaN(size(a_values));

fprintf('-----------------------------------------------------------\n');
fprintf('         a (m)            L_threshold (m)\n');
fprintf('-----------------------------------------------------------\n');

%% Loop over radius values a
for ia = 1:length(a_values)
    
    a = a_values(ia);
    
    % For this a, compute c_opt for each L
    c_opt_vec = NaN(size(L_values));
    
    for iL = 1:length(L_values)
        
        L = L_values(iL);
        
        % Leaf water potential (depends on L)
        u = m - rho * g * L;      % leaf water potential
        
        % Compute J(c)
        N = -m .* k .* L .* a.^2 .* R .* T .* (c.^2) ...
            - m .* u .* a.^2 .* k .* L .* c;
        
        D = R .* T .* c .* (a.^3 + 4 .* k .* nu .* L.^2) ...
            + m .* (a.^3 - 4 .* k .* nu .* L.^2);
        
        J = N ./ D;
        
        % Find optimal concentration c_opt for this (a, L)
        [~, idx_max] = max(J);
        c_opt_vec(iL) = c(idx_max);
        
    end
    
    % Indices where c_opt is below the threshold
    below_idx = find(c_opt_vec < c_threshold);
    
    if ~isempty(below_idx)
        % Choose the L whose c_opt is closest to c_threshold from below
        diff_to_thr = c_threshold - c_opt_vec(below_idx);
        [~, local_idx] = min(diff_to_thr);
        best_iL = below_idx(local_idx);
        L_threshold(ia) = L_values(best_iL);
    else
        % No c_opt below the threshold for this a
        L_threshold(ia) = NaN;
    end
    
    fprintf('  %12.4e\t    %12.4e\n', a, L_threshold(ia));
end

fprintf('-----------------------------------------------------------\n');

%% Plot L_threshold as a function of a

figure('Color','w');
loglog(a_values, L_threshold, 'o-', 'LineWidth', 1.6, 'MarkerSize', 6);
grid on; box on;

xlabel('Radius, a (m)', 'Color', 'k');
ylabel('L_{threshold} (m)', 'Color', 'k');
title(sprintf('L_{threshold} vs a (c_{threshold} = %.0f mol/m^3)', c_threshold), ...
      'Color', 'k');

ax = gca;
ax.XColor      = 'k';
ax.YColor      = 'k';
ax.Title.Color = 'k';
ax.GridColor   = [0.7 0.7 0.7];
ax.GridAlpha   = 0.4;
set(ax, 'Layer', 'top');
