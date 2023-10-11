clc; clear all; close all;
t0 = 0; tf = pi/2; nbrk = 5; k = 6;
breaks = linspace(t0,tf,nbrk);
tau = functiontau(breaks, 4);
knots = augknt(breaks, k);
colmat = spcol(knots, k, brk2knt(tau, (k-1)));
%%%
% y'^2+ y^2 = 50/(1 + sin(2x))^2
% Ey'' + 2yk'y' + 2yky = 50/(1 + sin(2x))^2 + yk'^2 + yk^2
% w0 = 2yk; w1 = 2yk'; w3=E; b(x) = 50/(1 + sin(2x))^2 + yk'^2 + yk^2
% w0y + w1y' + w2y'' = b(x)
%%%
IG = @(x) -5.*x.^2
IGd = @(x) 5.*x
epsilon=0; w0=2.*IG(tau).'; w1=2.*IGd(tau).'; b1=50./(1 + sin(2.*tau)).^2 + IG(tau).^2 + IGd(tau).^2; b1=b1.';
A= [10000.*colmat(1, :); 10000.*colmat((length(tau)-1)*(k-1)+1, :); w0.*colmat(1:(k-1):end, :)+w1.*colmat(2:(k-1):end, :)+epsilon.*colmat(3:(k-1):end, :)];
B = [10000.*5;10000.*5; b1];
coefs = A\B;
y = spmak(knots, coefs.');
xx = linspace(t0, tf, 10000);
tolerance=1e-2;
while 1
    ytau = fnval(y, tau);
    dy = fnder(y, 1); ydtau = fnval(dy, tau);
    epsilon=0; w0=2.*ytau.'; w1=2.*ydtau.'; b1=50./(1 + sin(2.*tau)).^2 + ytau.^2 + ydtau.^2; b1=b1.';
    A= [10000.*colmat(1, :); 10000.*colmat((length(tau)-1)*(k-1)+1, :); w0.*colmat(1:(k-1):end, :)+w1.*colmat(2:(k-1):end, :)+epsilon.*colmat(3:(k-1):end, :)];
    B = [10000.*5;10000.*5; b1];
    coefs = A\B;
    z = spmak(knots, coefs.');
    maxdif = max(abs(y.coefs - z.coefs))
    y = z;
    if maxdif<tolerance, break, end
end
sp = y;
figure; hold on;
plot(xx, 5.*sqrt(1./(1+sin(2.*xx))), 'r')
plot(xx, fnval(y, xx), 'b--')
legend('Solution', 'Approximation')
title('Non-linear experiment')
hold off;
IG = @(x) 5.*sqrt(1./(1+sin(2.*x)))
IGd = @(x) -5.*cos(2.*x)./(1 + sin(2.*x)).^(3/2)
Results.Domain = xx;
Results.Estimation = fnval(sp, xx);
Results.Solution = 5.*sqrt(1./(1+sin(2.*xx)));
Results.Error = abs(5.*sqrt(1./(1+sin(2.*xx))) - fnval(sp, xx));
Results.MaxError = max(Results.Error);
Results.MeanError = mean(Results.Error);
Results.StdError = std(Results.Error);
save('NL1D1V_1_Spline', '-struct', 'sp');
save('NL1D1V_1_Results', '-struct', 'Results');