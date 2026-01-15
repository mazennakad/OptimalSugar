clear; clc; close all;

% Conditions
T = 298.15;                     % Temperature [K]
c = linspace(0, 2000, 1000)';   % Concentration range [mol/m^3]

%  viscosity function
nu = Dynamic_Viscosity_T_Sucrose_CW(T, c);
nu = nu(:);  % ensure column vector

% Numerical derivative (2-point forward difference)
dnu_dc = zeros(size(nu));
for i = 1:length(c)-1
    dnu_dc(i) = (nu(i+1) - nu(i)) / (c(i+1) - c(i));
end
dnu_dc(end) = dnu_dc(end-1);  % repeat last value to keep same length

% Plot ν(c) and dν/dc(c) -----
figure('Color','w');
yyaxis left
plot(c, nu, 'b', 'LineWidth', 2);
ylabel('\nu (Pa·s)');
yyaxis right
plot(c, dnu_dc, 'r--', 'LineWidth', 2);
ylabel('d\nu/dc (Pa·s per mol·m^{-3})');

xlabel('Concentration, c (mol/m^3)');
title(sprintf('Sucrose solution: ν and dν/dc vs c at T = %.2f K', T));
legend('\nu (Pa·s)', 'd\nu/dc', 'Location', 'best');
grid on; box on;
print('myFigure','-djpeg','-r300')
