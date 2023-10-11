clc; clear all; close all;

%%% Fourth order spline
t0 = 0; tf = 2*pi;
nbrk = 5; k = 4; iv=1000;
breaks = linspace(t0,tf,nbrk);
tau = functiontau(breaks, 2);
knots = augknt(breaks, k);
colmat = spcol(knots, k, brk2knt(tau, 3));
A = [iv.*colmat(1, :); iv.*colmat(end-2, :); colmat(1:3:end, :) + colmat(3:3:end, :)];
B = [iv.*1; iv.*1; zeros(length(tau), 1)];
coef = A\B;
sp = spmak(knots, coef.');
spd1 = fnder(sp, 1);
spd2 = fnder(sp, 2);
xx = linspace(t0, tf, 1000);
figure; hold on;
plot(xx, fnval(sp, xx), 'b')
plot(xx, cos(xx), 'r--')
legend('Approximation', 'Solution')
title('Fourth order spline Cos approxiation')
hold off;
figure; hold on;
plot(xx, fnval(spd1, xx), 'b')
plot(xx, -sin(xx), 'r--')
legend('Approximation', 'Solution')
title('Cos first derivative approxiation')
hold off;
figure; hold on;
plot(xx, fnval(spd2, xx), 'b')
plot(xx, -cos(xx), 'r--')
legend('Approximation', 'Solution')
title('Cos decond derivative approxiation')
hold off;
Results.Domain = xx;
Results.Estimation = fnval(sp, xx);
Results.Estimationd1 = fnval(spd1, xx);
Results.Estimationd2 = fnval(spd2, xx);
Results.Solution = cos(xx);
Results.Solutiond1 = -sin(xx);
Results.Solutiond2 = -cos(xx);
Results.Error = abs(cos(xx) - fnval(sp, xx));
Results.Errord1 = abs(-sin(xx) - fnval(spd1, xx));
Results.Errord2 = abs(-cos(xx) - fnval(spd2, xx));
Results.MaxError = max(Results.Error);
Results.MeanError = mean(Results.Error);
Results.StdError = std(Results.Error);
Results.MaxErrord1 = max(Results.Errord1);
Results.MeanErrord1 = mean(Results.Errord1);
Results.StdErrord1 = std(Results.Errord1);
Results.MaxErrord2 = max(Results.Errord2);
Results.MeanErrord2 = mean(Results.Errord2);
Results.StdErrord2 = std(Results.Errord2);
save('L1D1V_3a_Spline', '-struct', 'sp');
save('L1D1V_3a_Results', '-struct', 'Results');


%%% Sixth order spline
t0 = 0; tf = 2*pi;
nbrk = 5; k = 6; iv=1;
breaks = linspace(t0,tf,nbrk);
tau = functiontau(breaks, 4);
knots = augknt(breaks, k);
colmat = spcol(knots, k, brk2knt(tau, 3));
A = [iv.*colmat(1, :); iv.*colmat(end-2, :); colmat(1:3:end, :) + colmat(3:3:end, :)];
B = [iv.*1; iv.*1; zeros(length(tau), 1)];
coef = A\B;
sp = spmak(knots, coef.');
spd1 = fnder(sp, 1);
spd2 = fnder(sp, 2);
xx = linspace(t0, tf, 1000);
figure; hold on;
plot(xx, fnval(sp, xx), 'b')
plot(xx, cos(xx), 'r--')
legend('Approximation', 'Solution')
title('Sixth order spline Cos approxiation')
hold off;
figure; hold on;
plot(xx, fnval(spd1, xx), 'b')
plot(xx, -sin(xx), 'r--')
legend('Approximation', 'Solution')
title('Cos first derivative approxiation')
hold off;
figure; hold on;
plot(xx, fnval(spd2, xx), 'b')
plot(xx, -cos(xx), 'r--')
legend('Approximation', 'Solution')
title('Cos decond derivative approxiation')
hold off;
Results.Domain = xx;
Results.Estimation = fnval(sp, xx);
Results.Estimationd1 = fnval(spd1, xx);
Results.Estimationd2 = fnval(spd2, xx);
Results.Solution = cos(xx);
Results.Solutiond1 = -sin(xx);
Results.Solutiond2 = -cos(xx);
Results.Error = abs(cos(xx) - fnval(sp, xx));
Results.Errord1 = abs(-sin(xx) - fnval(spd1, xx));
Results.Errord2 = abs(-cos(xx) - fnval(spd2, xx));
Results.MaxError = max(Results.Error);
Results.MeanError = mean(Results.Error);
Results.StdError = std(Results.Error);
Results.MaxErrord1 = max(Results.Errord1);
Results.MeanErrord1 = mean(Results.Errord1);
Results.StdErrord1 = std(Results.Errord1);
Results.MaxErrord2 = max(Results.Errord2);
Results.MeanErrord2 = mean(Results.Errord2);
Results.StdErrord2 = std(Results.Errord2);
save('L1D1V_3b_Spline', '-struct', 'sp');
save('L1D1V_3b_Results', '-struct', 'Results');