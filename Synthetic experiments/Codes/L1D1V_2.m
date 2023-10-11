clc; clear all; close all;

%%% Error depending on the breaks
t0 = 0; tf = 2*pi;
k = 6; iv=1;
error = zeros(3, 3);
for i=3:5
    breaks = linspace(t0,tf,i);
    tau = functiontau(breaks, 1);;
    knots = augknt(breaks, k);
    colmat = spcol(knots, k, brk2knt(tau, 3));
    A = [iv.*colmat(1, :); iv.*colmat(2, :); colmat(1:3:end, :) + colmat(3:3:end, :)];
    B = [0; iv.*1; zeros(length(tau), 1)];
    coef = A\B;
    sp = spmak(knots, coef.');
    xx = linspace(t0, tf, 1000);
    e = abs(sin(xx) - fnval(sp, xx));
    error(i-2,1) = mean(e);
    error(i-2,2) = max(e);
    error(i-2,3) = std(e);
end
Results.Errornbrks = error;

%%% Error depending on the number of collocation points
t0 = 0; tf = 2*pi;
nbrks = 5; k = 6; iv=1;
error = zeros(4, 1);
for i=1:4
    breaks = linspace(t0,tf,nbrks);
    tau = functiontau(breaks, i);
    knots = augknt(breaks, k);
    colmat = spcol(knots, k, brk2knt(tau, 3));
    A = [iv.*colmat(1, :); iv.*colmat(2, :); colmat(1:3:end, :) + colmat(3:3:end, :)];
    B = [0; iv.*1; zeros(length(tau), 1)];
    coef = A\B;
    sp = spmak(knots, coef.');
    xx = linspace(t0, tf, 1000);
    e = abs(sin(xx) - fnval(sp, xx));
    error(i) = mean(e);
    error(i,2) = max(e);
    error(i,3) = std(e);
end
figure; hold on;
plot(xx, fnval(sp, xx), 'b')
plot(xx, sin(xx), 'r--')
legend('Approximation', 'Solution')
title('Sin approxiation')
hold off;
Results.Errortau = error;
Results.Domain = xx;
Results.Estimation = fnval(sp, xx);
Results.Solution = sin(xx);
Results.Error = abs(sin(xx) - fnval(sp, xx));
Results.MaxError = max(Results.Error);
Results.MeanError = mean(Results.Error);
Results.StdError = std(Results.Error);
save('L1D1V_2_Spline', '-struct', 'sp');
save('L1D1V_2_Results', '-struct', 'Results');