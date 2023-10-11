clc; close all; clear all;
ui = 0; uf = 1200;
k = 6; nBrks = 10; nPBB = 4;
r = 0.006; C = 1e13; N0 = 1e9;
brksu = linspace(ui,uf,nBrks);
knotsu = augknt(brksu,k);
tauu = functiontau(brksu,nPBB);
colmat = spcol(knotsu, k, brk2knt(tauu, 3));
A = [colmat(1, :); colmat(2:3:end, :)];
B = [N0; r*C.*ones(length(tauu), 1)];
coefs = A\B;
sp = spmak(knotsu, coefs.');
Nk = fnval(sp, tauu).';
for i=1:5
    knotsu = newknt(sp);
    tauu = functiontau(brksu,nPBB);
    colmat = spcol(knotsu, k, brk2knt(tauu, 3));
    w = r*(log(C./Nk)-1);
    A = [colmat(1, :); colmat(2:3:end, :) - w.*colmat(1:3:end, :)];
    B = [N0; r.*Nk];
    coefs = A\B;
    sp = spmak(knotsu, coefs.');
    Nk = fnval(sp, tauu).';
end
N = @(t) C.*exp(log(N0/C).*exp(-r.*t))
ux = linspace(ui, uf, 10000);
figure; hold on;
plot(ux, N(ux), 'r')
plot(ux, fnval(sp, ux), 'b--')
legend('Ground truth', 'Approximation')
title('Gompertz')
hold off;
error = abs(fnval(sp, ux) - N(ux))./N(ux);
[max(error) mean(error) std(error)]
error1 = abs(fnval(sp, ux) - N(ux));
[max(error1) mean(error1) std(error1)]