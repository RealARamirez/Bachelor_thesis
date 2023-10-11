clc; clear all; close all;
t0 = -1; tf = 1; nbrk = 5; k = 6;
breaks = linspace(t0,tf,nbrk);
xx = linspace(t0, tf, 10000);
tau = functiontau(breaks, 4);
knots = augknt(breaks, k);
colmat = spcol(knots, k, brk2knt(tau, (k-1)));
colmativ = spcol(knots, k, [-1 1]);
%% Eikonal: epsilon=0
A = [10000.*colmativ; 2.*colmat(2:(k-1):end, :)];
B = [[0 0].'; ones(length(tau), 1)];
coefs = A\B;
sp = spmak(knots, coefs.');
spd1 = fnder(sp, 1);
figure; hold on
plot(xx, fnval(sp, xx), 'b')
plot(xx,  1 + xx.*(xx<0) - xx.*(xx>0), 'r')
legend('Initial guess', 'Eikonal')
title('Initial guess')
hold off
for i=1:10
    A = [10000.*colmativ; fnval(spd1, tau).'.*2.*colmat(2:(k-1):end, :)];
    B = [[0 0].'; fnval(spd1, tau).'.^2+ones(length(tau), 1)];
    coefs = A\B;
    sp = spmak(knots, coefs.');
    spd1 = fnder(sp, 1);
end
Results.Domain = xx;
Results.Estimation = fnval(sp, xx);
Results.Solution = 1 + xx.*(xx<0) - xx.*(xx>0);
Results.Error = abs(1 + xx.*(xx<0) - xx.*(xx>0) - fnval(sp, xx));
Results.MaxError = max(Results.Error);
Results.MeanError = mean(Results.Error);
Results.StdError = std(Results.Error);
save('NL1D1V_2a_Spline', '-struct', 'sp');
save('NL1D1V_2a_Results', '-struct', 'Results');
figure; hold on
plot(xx, fnval(sp, xx), 'b')
plot(xx,  1 + xx.*(xx<0) - xx.*(xx>0), 'r')
legend('Approximation', 'Eikonal')
title('Epsilon 0')
hold off
%% Eikonal: epsilon=0.25
epsilon=0.25;
A = [10000.*colmativ; 2.*colmat(2:(k-1):end, :)-epsilon.*colmat(3:(k-1):end, :)];
B = [[0 0].'; ones(length(tau), 1)];
coefs = A\B;
sp = spmak(knots, coefs.');
spd1 = fnder(sp, 1);
for i=1:10
    A = [10000.*colmativ; fnval(spd1, tau).'.*2.*colmat(2:(k-1):end, :)-epsilon.*colmat(3:(k-1):end, :)];
    B = [[0 0].'; fnval(spd1, tau).'.^2+ones(length(tau), 1)];
    coefs = A\B;
    sp = spmak(knots, coefs.');
    spd1 = fnder(sp, 1);
end
Results.Domain = xx;
Results.Estimation = fnval(sp, xx);
Results.Solution = epsilon.*(log(cosh(1/epsilon))-log(cosh(xx./epsilon)));
Results.Error = abs(epsilon.*(log(cosh(1/epsilon))-log(cosh(xx./epsilon))) - fnval(sp, xx));
Results.MaxError = max(Results.Error);
Results.MeanError = mean(Results.Error);
Results.StdError = std(Results.Error);
save('NL1D1V_2b_Spline', '-struct', 'sp');
save('NL1D1V_2b_Results', '-struct', 'Results');
figure; hold on
plot(xx,  epsilon.*(log(cosh(1/epsilon))-log(cosh(xx./epsilon))), 'r')
plot(xx, fnval(sp, xx), 'b--')
legend('Approximation', 'Solution')
title('Epsilon 0.25')
hold off
%% Eikonal: epsilon=0.125
epsilon=0.125;
A = [10000.*colmativ; 2.*colmat(2:(k-1):end, :)-epsilon.*colmat(3:(k-1):end, :)];
B = [[0 0].'; ones(length(tau), 1)];
coefs = A\B;
sp = spmak(knots, coefs.');
spd1 = fnder(sp, 1);
for i=1:10
    A = [10000.*colmativ; fnval(spd1, tau).'.*2.*colmat(2:(k-1):end, :)-epsilon.*colmat(3:(k-1):end, :)];
    B = [[0 0].'; fnval(spd1, tau).'.^2+ones(length(tau), 1)];
    coefs = A\B;
    sp = spmak(knots, coefs.');
    spd1 = fnder(sp, 1);
end
Results.Domain = xx;
Results.Estimation = fnval(sp, xx);
Results.Solution = epsilon.*(log(cosh(1/epsilon))-log(cosh(xx./epsilon)));
Results.Error = abs(epsilon.*(log(cosh(1/epsilon))-log(cosh(xx./epsilon))) - fnval(sp, xx));
Results.MaxError = max(Results.Error);
Results.MeanError = mean(Results.Error);
Results.StdError = std(Results.Error);
save('NL1D1V_2c_Spline', '-struct', 'sp');
save('NL1D1V_2c_Results', '-struct', 'Results');
figure; hold on
plot(xx,  epsilon.*(log(cosh(1/epsilon))-log(cosh(xx./epsilon))), 'r')
plot(xx, fnval(sp, xx), 'b--')
legend('Approximation', 'Solution')
title('Epsilon 0.125')
hold off
%% Eikonal: epsilon=0.125 and newknt
epsilon=0.125;
knots = newknt(sp, length(breaks)+3);
breaks = knt2brk(knots);
tau = functiontau(breaks, 4);
colmat = spcol(knots, k, brk2knt(tau, (k-1)));
colmativ = spcol(knots, k, [-1 1]);
A = [10000.*colmativ; 2.*colmat(2:(k-1):end, :)-epsilon.*colmat(3:(k-1):end, :)];
B = [[0 0].'; ones(length(tau), 1)];
coefs = A\B;
sp = spmak(knots, coefs.');
spd1 = fnder(sp, 1);
for i=1:10
    A = [10000.*colmativ; fnval(spd1, tau).'.*2.*colmat(2:(k-1):end, :)-epsilon.*colmat(3:(k-1):end, :)];
    B = [[0 0].'; fnval(spd1, tau).'.^2+ones(length(tau), 1)];
    coefs = A\B;
    sp = spmak(knots, coefs.');
    spd1 = fnder(sp, 1);
end
Results.Domain = xx;
Results.Estimation = fnval(sp, xx);
Results.Solution = epsilon.*(log(cosh(1/epsilon))-log(cosh(xx./epsilon)));
Results.Error = abs(epsilon.*(log(cosh(1/epsilon))-log(cosh(xx./epsilon))) - fnval(sp, xx));
Results.MaxError = max(Results.Error);
Results.MeanError = mean(Results.Error);
Results.StdError = std(Results.Error);
save('NL1D1V_2d_Spline', '-struct', 'sp');
save('NL1D1V_2d_Results', '-struct', 'Results');
figure; hold on
plot(xx,  epsilon.*(log(cosh(1/epsilon))-log(cosh(xx./epsilon))), 'r')
plot(xx, fnval(sp, xx), 'b--')
legend('Approximation', 'Solution')
title('Epsilon 0.125 knots redistribution')
hold off
%% Eikonal: epsilon=0
for i=1:10
    A = [10000.*colmativ; fnval(spd1, tau).'.*2.*colmat(2:(k-1):end, :)];
    B = [[0 0].'; fnval(spd1, tau).'.^2+ones(length(tau), 1)];
    coefs = A\B;
    sp = spmak(knots, coefs.');
    spd1 = fnder(sp, 1);
end
Results.Domain = xx;
Results.Estimation = fnval(sp, xx);
Results.Solution = 1 + xx.*(xx<0) - xx.*(xx>0);
Results.Error = abs(1 + xx.*(xx<0) - xx.*(xx>0) - fnval(sp, xx));
Results.MaxError = max(Results.Error);
Results.MeanError = mean(Results.Error);
Results.StdError = std(Results.Error);
save('NL1D1V_2e_Spline', '-struct', 'sp');
save('NL1D1V_2e_Results', '-struct', 'Results');
figure; hold on
plot(xx,  1 + xx.*(xx<0) - xx.*(xx>0), 'r')
plot(xx, fnval(sp, xx), 'b--')
legend('Approximation', 'Eikonal')
title('Epsilon 0 with previous initial guess')
hold off