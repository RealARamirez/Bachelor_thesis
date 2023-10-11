clc; clear all; close all;

%%% Breaks as collocation points
t0 = 0; tf = 1;
nbrk = 4; k = 4; iv=1;
breaks = linspace(t0,tf,nbrk);
tau = breaks;
knots = augknt(breaks, k);
colmat = spcol(knots, k, brk2knt(tau, 3));
A = [iv.*colmat(1, :); colmat(1:3:end, :) - colmat(2:3:end, :)];
B = [iv.*1; zeros(length(tau), 1)];
coef = A\B;
sp = spmak(knots, coef.');
xx = linspace(t0, tf, 1000);
figure; hold on;
plot(xx, fnval(sp, xx), 'b')
plot(xx, exp(xx), 'r--')
legend('Approximation', 'Solution')
title('Exponential approxiation')
hold off;
Results.Domain = xx;
Results.Estimation = fnval(sp, xx);
Results.Solution = exp(xx);
Results.Error = abs(exp(xx) - fnval(sp, xx));
Results.MaxError = max(Results.Error);
Results.MeanError = mean(Results.Error);
Results.StdError = std(Results.Error);
save('L1D1V_1a_Spline', '-struct', 'sp');
save('L1D1V_1a_Results', '-struct', 'Results');

%%% Breaks as collocation points
t0 = 0; tf = 1;
nbrk = 4; k = 4; iv=1;
breaks = linspace(t0,tf,nbrk);
tau = functiontau(breaks, 1);
knots = augknt(breaks, k);
colmat = spcol(knots, k, brk2knt(tau, 3));
A = [iv.*colmat(1, :); colmat(1:3:end, :) - colmat(2:3:end, :)];
B = [iv.*1; zeros(length(tau), 1)];
coef = A\B;
sp = spmak(knots, coef.');
xx = linspace(t0, tf, 1000);
figure; hold on;
plot(xx, fnval(sp, xx), 'b')
plot(xx, exp(xx), 'r--')
legend('Approximation', 'Solution')
title('Exponential approxiation')
hold off;
Results.Domain = xx;
Results.Estimation = fnval(sp, xx);
Results.Solution = exp(xx);
Results.Error = abs(exp(xx) - fnval(sp, xx));
Results.MaxError = max(Results.Error);
Results.MeanError = mean(Results.Error);
Results.StdError = std(Results.Error);
save('L1D1V_1b_Spline', '-struct', 'sp');
save('L1D1V_1b_Results', '-struct', 'Results');