
clear; clc; close all;

% Parameters
R = 8.314;          % J/(mol·K)
T = 298.15;         % Temperature
rho = 1000;  % Density of water in kg/m^3
g = 9.81;    % Acceleration due to gravity in m/s^2
m = -0.001e6;       % soil water potential
L = 1.0;            % Phloem (sieze tube) Length
a = 1.0e-5;         % Phloem (sieze tube) raduis
u = m - rho*g*L;    % leaf water potential
k = 1.0e-12;        % permeability 



% Concentration range
c_min = 0;
c_max = 2000;
num_points = 4000;
c = linspace(c_min, c_max, num_points)';   % mol/m^3

% Dynamic viscosity 
nu = Dynamic_Viscosity_T_Sucrose_CW(T, c);
nu = nu(:);   % ensure column vector

% Compute J(c)
N = -m .* k .* L .* a.^2 .* R .* T .* (c.^2) - m .* u .* a.^2 .* k .* L .* c;
D = R .* T .* c .* (a.^3 + 4 .* k .* nu .* L.^2) + m .* (a.^3 - 4 .* k .* nu .* L.^2);
J = N ./ D;

% Find maximum J and c_opt
J_max = J(1);
c_at_Jmax = c(1);

for i = 2:length(J)
    if J(i) > J_max
        J_max = J(i);
        c_at_Jmax = c(i);
    end
end

% Plot 
figure;
plot(c, J, 'b', 'LineWidth', 2); hold on;
plot(c_at_Jmax, J_max, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % mark max point
xlabel('Concentration, c (mol/m^3)');
ylabel('Flow, J(c)');
title('J as a function of c at T = 298.15 K');
legend('J(c)', sprintf('J_{max} = %.2e at c = %.1f', J_max, c_at_Jmax), 'Location', 'best');
grid on; box on;

% Display result
fprintf('Maximum J = %.4e at c = %.2f mol/m^3\n', J_max, c_at_Jmax);
print('myFigure','-djpeg','-r300')
