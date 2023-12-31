clc; clear all; close all;
t0 = 0; tf = 1;
nBrks = 10; k = 6; nPBB = 4;
brksu = linspace(t0,tf,nBrks); brksv = linspace(t0,tf,nBrks);
knotsu = augknt(brksu,k); knotsv = augknt(brksv,k);
tauu = functiontau(brksu,nPBB); tauv = functiontau(brksv,nPBB);
[Ugr,Vgr] = meshgrid(tauu,tauv); Ugr=Ugr(:); Vgr=Vgr(:);
colmatu = spcol(knotsu, k, brk2knt(tauu, 3), 'sparse');
colmatv = spcol(knotsv, k, brk2knt(tauv, 3), 'sparse');
colmat = kron(colmatv(1:3:end, :), colmatu(1:3:end, :));
B1 = Ugr; B2 = Vgr;
B3 = ones(length(Ugr), 1) - Ugr;
B = [B1 B2 B3];
B = [B1 B2 B3];
coef = colmat\B;
coef = reshape(coef.', [3 length(knotsu)-k length(knotsv)-k]);
sp = spmak({knotsu, knotsv}, coef);
uev = linspace(t0, tf, 50); vev = linspace(t0, tf, 50); [Uevg,Vevg] = meshgrid(uev,uev);
outcome = fnval(sp, {uev, vev});
figure;
surf(squeeze(outcome(1, :, :)), squeeze(outcome(2, :, :)), squeeze(outcome(3, :, :)))
title('Approximation');
xlabel('x-axis')
ylabel('y-axis')
Results.Domainu = uev;
Results.Domainv = vev;
Results.Estimationx = squeeze(outcome(1, :, :));
Results.Estimationy = squeeze(outcome(2, :, :));
Results.Estimationz = squeeze(outcome(3, :, :));
Results.Solutionx = Uevg;
Results.Solutiony = Vevg;
Results.Errorx = abs(Uevg - squeeze(outcome(1, :, :)));
Results.Errory = abs(Vevg - squeeze(outcome(2, :, :)));
Results.Errorz = abs(ones(size(Uevg)) - Uevg - squeeze(outcome(3, :, :)));
Results.MaxErrorx = max(Results.Errorx, [], 'all');
Results.MeanErrorx = mean(Results.Errorx, 'all');
Results.StdErrorx = std(Results.Errorx, [], 'all');
Results.MaxErrory = max(Results.Errory, [], 'all');
Results.MeanErrory = mean(Results.Errory, 'all');
Results.StdErrory = std(Results.Errory, [], 'all');
Results.MaxErrorz = max(Results.Errorz, [], 'all');
Results.MeanErrorz = mean(Results.Errorz, 'all');
Results.StdErrorz = std(Results.Errorz, [], 'all');
save('LMDMV_2_Spline', '-struct', 'sp');
save('LMDMV_2_Results', '-struct', 'Results');