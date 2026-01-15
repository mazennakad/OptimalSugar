%%
clear; clc; close all;

% Fixed constants
R   = 8.314;          % J/(mol·K)
T   = 298.15;         % K
m   = -0.001e6;       % Pa (soil water potential)
k   = 1.0e-15;        % fixed permeability (m^2)
rho = 1000;           % kg/m^3 (density of water)
g   = 9.81;           % m/s^2 (gravity)

% c-grid 
c  = linspace(0, 2000, 400)';     % mol/m^3
nu = Dynamic_Viscosity_T_Sucrose_CW(T, c); 
nu = nu(:);                       % ensure column

% Parameter grids (5 by 5)
a_values = linspace(5e-6, 3e-5, 5);   % m
L_values = linspace(0.1, 20, 5);      % m

% Storage: [L, a, c_max, J_max] 
results = zeros(numel(L_values)*numel(a_values), 4);
row = 0;

% Plot setup
figure('Color','w'); 
hold on; grid on; box on;
xlabel('Concentration, c (mol/m^3)', 'Color', 'k');
ylabel('J(c)', 'Color', 'k');
title(sprintf('J vs c for 5×5 grid of (L, a) with k = %.0e', k), 'Color', 'k');

colors = lines(numel(L_values)*numel(a_values));
curve_idx = 0;

ax = gca;
ax.XColor       = 'k';
ax.YColor       = 'k';
ax.Title.Color  = 'k';
ax.XLabel.Color = 'k';
ax.YLabel.Color = 'k';
ax.GridColor    = [1 1 1];   % white grid lines
ax.GridAlpha    = 0.2;
set(ax,'Layer','top');       % draw grid on top of curves

% Sweep L and a 
for ia = 1:numel(a_values)
    a = a_values(ia);
    for iL = 1:numel(L_values)
        L = L_values(iL);
        curve_idx = curve_idx + 1;

        % Leaf water potential depends on L
        u = m - rho*g*L;    % Pa

        % J(c)
        N = -m .* k .* L .* a.^2 .* R .* T .* (c.^2) ...
            - m .* u .* a.^2 .* k .* L .* c;
        D = R .* T .* c .* (a.^3 + 4 .* k .* nu .* L.^2) ...
            + m .* (a.^3 - 4 .* k .* nu .* L.^2);
        J = N ./ D;

        % Find J_max & c_max via for+if
        J_max = J(1);  
        c_max = c(1);
        for j = 2:numel(J)
            if J(j) > J_max
                J_max = J(j);
                c_max = c(j);
            end
        end

        % Store and plot
        row = row + 1;
        results(row, :) = [L, a, c_max, J_max];

        plot(c, J, 'Color', colors(curve_idx,:), 'LineWidth', 1.3);
        plot(c_max, J_max, 'o', 'MarkerSize', 5, ...
             'MarkerFaceColor', colors(curve_idx,:), ...
             'MarkerEdgeColor', colors(curve_idx,:));
    end
end

% Print all combinations
fprintf('All (L, a) combinations with their maxima:\n');
fprintf('---------------------------------------------------------------\n');
fprintf('   L (m)        a (m)       c_max (mol/m^3)       J_max\n');
fprintf('---------------------------------------------------------------\n');
for r = 1:size(results,1)
    fprintf('%9.3f   %9.3e      %10.3f         % .4e\n', ...
        results(r,1), results(r,2), results(r,3), results(r,4));
end
fprintf('---------------------------------------------------------------\n');

% Top 5 by J_max
[~, ord] = sort(results(:,4), 'descend');
top5 = results(ord(1:5), :);

fprintf('\nTop 5 combinations by J_max (L, a):\n');
fprintf('---------------------------------------------------------------\n');
fprintf('   L (m)        a (m)       c_max (mol/m^3)       J_max\n');
fprintf('---------------------------------------------------------------\n');
for r = 1:5
    fprintf('%9.3f   %9.3e      %10.3f         % .4e\n', ...
        top5(r,1), top5(r,2), top5(r,3), top5(r,4));
end
fprintf('---------------------------------------------------------------\n');
