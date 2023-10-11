clc; clear all; close all;
t0 = 0; tf = 1;
nBrks = 5; k = 6; nPBB = 4;
brksu = linspace(t0,tf,nBrks);
knotsu = augknt(brksu,k);
tauu = functiontau(brksu,nPBB);
colmat = spcol(knotsu, k, brk2knt(tauu, 3));
[rows ncoefs] = size(colmat);
A1 = [colmat(1, :); colmat(1:3:end, :) - colmat(2:3:end, :)];
A2 = [zeros(1, ncoefs); colmat(1:3:end, :)];
A3 = [zeros(1, ncoefs); 4.*colmat(1:3:end, :)];
A4 = [colmat(1, :); -2.*colmat(1:3:end, :) - colmat(2:3:end, :)];
B1 = [2; zeros(length(tauu), 1)];
B2 = [-3; zeros(length(tauu), 1)];
A = [A1 A2; A3 A4];
B = [B1; B2];
coeffs = A\B;
coefs = zeros(2, ncoefs);
coefs(1, :) = coeffs(1:ncoefs);
coefs(2, :) = coeffs(1+ncoefs:end);
sp = spmak(knotsu, coefs);
N = 50; uev = linspace(t0,tf,N);
outcome = fnval(sp, uev);
Solx = @(t) exp(2.*t)+exp(-3.*t);
Soly = @(t) exp(2.*t)-4.*exp(-3.*t);
Results.Domain = uev;
Results.Estimationx = outcome(1, :);
Results.Estimationy = outcome(2, :);
Results.Solutionx = Solx(uev);
Results.Solutiony = Soly(uev);
Results.Errorx = abs(Solx(uev) - outcome(1, :));
Results.Errory = abs(Soly(uev) - outcome(2, :));
Results.MaxErrorx = max(Results.Errorx);
Results.MeanErrorx = mean(Results.Errorx);
Results.StdErrorx = std(Results.Errorx);
Results.MaxErrory = max(Results.Errory);
Results.MeanErrory = mean(Results.Errory);
Results.StdErrory = std(Results.Errory);
figure; hold on;
plot(uev, outcome(1, :), 'b')
plot(uev, Solx(uev), 'r--')
title('Solution x(t)')
xlabel('t')
ylabel('x(t)')
legend('Approximation', 'Analytical solution')
hold off;
figure; hold on;
plot(uev, outcome(2, :), 'b')
plot(uev, Soly(uev), 'r--')
title('Solution y(t)')
xlabel('t')
ylabel('y(t)')
legend('Approximation', 'Analytical solution')
hold off;
figure; hold on;
plot(outcome(1, :), outcome(2, :), 'b')
plot(Solx(uev), Soly(uev), 'r--')
xlim([min(outcome(1, :)-0.015) max(outcome(1, :))])
title('Solution y(x)')
xlabel('x(t)')
ylabel('y(t)')
legend('Approximation', 'Analytical solution')
hold off;
save('LMD1V_1_Spline', '-struct', 'sp');
save('LMD1V_1_Results', '-struct', 'Results');