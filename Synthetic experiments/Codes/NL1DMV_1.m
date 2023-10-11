clc; clear all; close all;
ui = 0; uf = 1; vi = 0; vf = 1; 
k = 6; nBrks = 5; nPBB = 4;
brksu = linspace(ui,uf,nBrks); brksv = linspace(vi,vf,nBrks);
knotsu = augknt(brksu,k); knotsv = augknt(brksv,k);
tauu = functiontau(brksu,nPBB); tauv = functiontau(brksv,nPBB);
colmatu = spcol(knotsu,k,brk2knt(tauu,3),'sparse');
colmatv = spcol(knotsv,k,brk2knt(tauv,3),'sparse');
colmatvi = spcol(knotsu, k, vi, 'sparse');
colmatui = spcol(knotsu, k, ui, 'sparse');
colmatu0 = colmatu(1:3:end,:); colmatv0 = colmatv(1:3:end,:);
colmatu1 = colmatu(2:3:end,:); colmatv1 = colmatv(2:3:end,:);
colmatui1 = kron(colmatu0, colmatvi);
colmatvi1 = kron(colmatui, colmatv0);
colmatu0v1 = kron(colmatu0,colmatv1);
colmatu0v0 = kron(colmatu0,colmatv0);
colmatu1v0 = kron(colmatu1,colmatv0);
IC = @(u) u;
IC1 = kron(IC(tauu).', ones(length(tauv), 1));
iv = 10000;
CI1 = tauu;
B1 = IC1;
B = [iv.*CI1.'; B1];
A = [iv.*colmatui1; colmatu0v1 + IC1.*colmatu1v0 + colmatu0v0];
coeffs = A \ B;
coeffs = reshape(coeffs, [length(knotsu)-k length(knotsv)-k]);
sp = spmak({knotsu, knotsv}, coeffs);
spd1 = fnder(sp, [0 1]);
for i=1:4
    IC1 = fnval(sp, {tauu, tauv}); IC1=reshape(IC1, [], 1);
    IC1d = fnval(spd1, {tauu, tauv}); IC1d=reshape(IC1d, [], 1);
    B1 = IC1.*IC1d;
    B = [iv.*CI1.'; B1];
    A = [iv.*colmatui1;  colmatu0v1 + IC1.*colmatu1v0 + IC1d.*colmatu0v0];
    coeffs = A \ B;
    max(abs(reshape(coeffs, [], 1) - reshape(sp.coefs, [], 1)))
    coeffs = reshape(coeffs, [length(knotsu)-k length(knotsv)-k]);
    sp = spmak({knotsu, knotsv}, coeffs);
    spd1 = fnder(sp, [0 1]);
end
N = 50; uev = linspace(ui,uf,N); vev = linspace(vi,vf,N);
[Uevg,Vevg] = ndgrid(uev,uev);Uev = Uevg(:); Vev=Vevg(:);
Results.Domainu = Uevg;
Results.Domainv = Vevg;
Results.Estimation = fnval(sp, {uev, vev});
Results.Estimationdx = fnval(spd1, {uev, vev});
Results.Solution = Vevg./(Uevg+1);
Results.Solutiondx = 1./(Uevg+1);
Results.Error = abs(Vevg./(Uevg+1) - fnval(sp, {uev, vev}));
Results.Errordx = abs(1./(Uevg+1) - fnval(spd1, {uev, vev}));
Results.MaxError = max(Results.Error, [], 'all');
Results.MeanErrordx = mean(Results.Error, 'all');
Results.StdError = std(Results.Error, [], 'all');
Results.MaxErrordx = max(Results.Errordx, [], 'all');
Results.MeanError = mean(Results.Errordx, 'all');
Results.StdErrordx = std(Results.Errordx, [], 'all');
save('NL1DMV_1_Spline', '-struct', 'sp');
save('NL1DMV_1_Results', '-struct', 'Results');
figure;
surf(Uevg,Vevg,Vevg./(Uevg+1));
title('Real')
xlabel('t-axis')
ylabel('x-axis')
figure;
surf(Uevg,Vevg,fnval(sp, {uev, vev}));
title('Approximation');
xlabel('t-axis')
ylabel('x-axis')
figure;
surf(Uevg,Vevg,1./(Uevg+1));
title('Real')
xlabel('t-axis')
ylabel('x-axis')
figure;
surf(Uevg,Vevg,fnval(spd1, {uev, vev}));
title('Approximation');
xlabel('t-axis')
ylabel('x-axis')