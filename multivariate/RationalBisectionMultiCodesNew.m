%% Part 2: Driver Script to Load KdV Data and Run Approximation
clear; close all; 

load KdV_sine.mat

N1 = 20; 
N2 = 20; 
eps = 1e-3; 
s = 10;

xx = t(1,1:s:end);
mm = x(1,1:s:end);
uu = usol(1:s:end,1:s:end);

xr = t(1,:);
mr = x(1,:);

T1 = xx;
T2 = mm;
d1 = length(T1);
d2 = length(T2);

T = 1:d1*d2;

B = reshape(uu.', 1, []);
f = B;

B1 = reshape(usol.', 1, []);
f1 = B1;

% Run Rational Bisection Method
tic;
[p, q, g, Err] = RationalBisectionMultiNew(f, T, T1, T2, d1, eps, N1, N2, xr, mr, f1);
toc;

fprintf('Max Error: %.8f\n', max(Err(:)));
disp('Numerator coefficients:'); disp(p);
disp('Denominator coefficients:'); disp(q);

[X, Y] = meshgrid(T1, T2);

figure;
surf(X, Y, uu);
xlabel('x'); ylabel('t');
shading interp; colorbar;
title('Original KdV Solution');
set(gca, 'FontSize', 16);
set(gcf, 'Color', 'w');
print(gcf, '-dpdf', 'KdV_solution_3D_10_points');
