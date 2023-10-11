clc; clear all; close all;
t0 = 0.01; tf = 1;
nBrks = 10; k = 6; nPBB = 4;
brksu = linspace(t0,tf,nBrks); brksv = linspace(t0,tf,nBrks);
knotsu = augknt(brksu,k); knotsv = augknt(brksv,k);
tauu = functiontau(brksu,nPBB); tauv = functiontau(brksv,nPBB);
[Ugr,Vgr] = meshgrid(tauu,tauv); Ugr=Ugr(:); Vgr=Vgr(:);
colmatu = spcol(knotsu, k, brk2knt(tauu, 3), 'sparse');
colmatv = spcol(knotsv, k, brk2knt(tauv, 3), 'sparse');
colmat = kron(colmatv(1:3:end, :), colmatu(1:3:end, :));
B1 = sqrt(Ugr.^2 + Vgr.^2);
B2 = atan(Vgr./Ugr);
B = [B1 B2];
coef = colmat\B;
coef = reshape(coef.', [2 length(knotsu)-k length(knotsv)-k]);
sp = spmak({knotsu, knotsv}, coef);
N = 50; uev = linspace(t0,tf,N); vev = linspace(t0,tf,N); [Uevg,Vevg] = meshgrid(uev,uev);
outcome = fnval(sp, {uev, vev});
Results.Domainu = Uevg;
Results.Domainv = Vevg;
Results.Estimationx = outcome(1, :);
Results.Estimationy = outcome(2, :);
Results.Solutionx = sqrt(Uevg.^2 + Vevg.^2);
Results.Solutiony = atan(Vevg./Uevg);
Results.Errorx = abs(sqrt(Uevg.^2 + Vevg.^2) - squeeze(outcome(1, :, :)));
Results.Errory = abs(atan(Vevg./Uevg)- squeeze(outcome(2, :, :)));
Results.MaxErrorx = max(Results.Errorx, [], 'all');
Results.MeanErrorx = mean(Results.Errorx, 'all');
Results.StdErrorx = std(Results.Errorx, [], "all");
Results.MaxErrory = max(Results.Errory, [], 'all');
Results.MeanErrory = mean(Results.Errory, 'all');
Results.StdErrory = std(Results.Errory, [], 'all');
figure; grid on;
surf(Uevg,Vevg, squeeze(outcome(1, :, :)));
title('Approximation')
xlabel('x-axis')
ylabel('y-axis')
figure; grid on;
surf(Uevg,Vevg, squeeze(outcome(2, :, :)));
title('Approximation')
xlabel('x-axis')
ylabel('y-axis')
figure; grid on;
surf(Uevg,Vevg, sqrt(Uevg.^2 + Vevg.^2));
title('Solution')
xlabel('x-axis')
ylabel('y-axis')
figure; grid on;
surf(Uevg,Vevg, atan(Vevg./Uevg));
title('Solution')
xlabel('x-axis')
ylabel('y-axis')
save('LMDMV_1_Spline', '-struct', 'sp');
save('LMDMV_1_Results', '-struct', 'Results');