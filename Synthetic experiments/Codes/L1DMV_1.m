clc; clear all; close all;
ui = 0; uf = 1; vi = 0; vf = 1; 
k = 6; nBrks = 10; nPBB = 4;
brksu = linspace(ui,uf,nBrks); brksv = linspace(vi,vf,nBrks);
knotsu = augknt(brksu,k); knotsv = augknt(brksv,k);
tauu = functiontau(brksu,nPBB); tauv = functiontau(brksv,nPBB);
[Ugr,Vgr] = meshgrid(tauu,tauv); Ugr=Ugr(:); Vgr=Vgr(:);
colmatu = spcol(knotsu,k,brk2knt(tauu,3),'sparse');
colmatv = spcol(knotsv,k,brk2knt(tauv,3),'sparse');
colmatu0 = colmatu(1:3:end,:); colmatv0 = colmatv(1:3:end,:);
colmatu1 = colmatu(2:3:end,:); colmatv1 = colmatv(2:3:end,:);
colmatu2 = colmatu(3:3:end,:); colmatv2 = colmatv(3:3:end,:);
colmatu0v0 = kron(colmatu0,colmatv0); colmatu0v1 = kron(colmatu0,colmatv1); colmatu0v2 = kron(colmatu0,colmatv2);
colmatu1v0 = kron(colmatu1,colmatv0); colmatu1v1 = kron(colmatu1,colmatv1); colmatu1v2 = kron(colmatu1,colmatv2);
colmatu2v0 = kron(colmatu2,colmatv0); colmatu2v1 = kron(colmatu2,colmatv1); colmatu2v2 = kron(colmatu2,colmatv2);
% Initial Condition
colmatu_ic1 = spcol(knotsu,k,tauu,'sparse');
colmatv_ic1i = spcol(knotsu,k,ui,'sparse');
colmatv_ic1f = spcol(knotsu,k,uf,'sparse');
colmatuv_ic1i = kron(colmatu_ic1,colmatv_ic1i);
colmatuv_ic1f = kron(colmatu_ic1,colmatv_ic1f);
colmatu_ic2i = spcol(knotsu,k,vi,'sparse');
colmatu_ic2f = spcol(knotsu,k,vf,'sparse');
colmatv_ic2 = spcol(knotsu,k,tauv,'sparse');
colmatuv_ic2i = kron(colmatu_ic2i,colmatv_ic2);
colmatuv_ic2f = kron(colmatu_ic2f,colmatv_ic2);
CI = [zeros(size(colmatuv_ic1i,1),1);zeros(size(colmatuv_ic1f,1),1);zeros(size(colmatuv_ic2i,1),1);zeros(size(colmatuv_ic2f,1),1)];
CI((2*length(tauv)+1):(3*length(tauv)), :) = (1 - tauv).*sin(5.*tauv);
CI = [zeros(size(colmatuv_ic1i,1),1);zeros(size(colmatuv_ic1f,1),1);zeros(size(colmatuv_ic2i,1),1);zeros(size(colmatuv_ic2f,1),1)];
CI((2*length(tauv)+1):(3*length(tauv)), :) = (1 - tauv).*sin(5.*tauv);
B1 = 10.*(Ugr-1).*cos(5.*Vgr) - 25.*(Ugr-1).*(Vgr-1).*sin(5.*Vgr);
B = [CI;B1];
A = [colmatuv_ic1i;colmatuv_ic1f;colmatuv_ic2i;colmatuv_ic2f;colmatu2v0+colmatu0v2];
coeffs = A \ B;
coef =  reshape(coeffs, [length(knotsu)-k length(knotsv)-k]);
sp = spmak({knotsu, knotsv}, coef);
N = 50; uev = linspace(ui,uf,N); vev = linspace(vi,vf,N); [Uevg,Vevg] = meshgrid(uev,uev);
Results.Domainu = uev;
Results.Domainv = vev;
Results.Estimation = fnval(sp, {uev, vev});
figure; grid on;
surf(Uevg,Vevg,fnval(sp, {uev, vev}));
title('Approximation')
xlabel('x-axis')
ylabel('y-axis')
%% LAPLACE: Teórica
syms u v z(u,v)
ui = 0; uf = 1; vi = 0; vf = 1; 
CI = 1; CId=0;
solp = (1-u)*(1-v)*sin(5*v);
N = 1e3; uev = linspace(ui,uf,N); vev = linspace(vi,vf,N);
[Uev,Vev] = meshgrid(uev,uev);
h1 = figure; hold on; grid on;
fsurf(solp,[ui,uf,vi,vf]);
title('Solution')
%% Results
N = 50; uev = linspace(ui,uf,N); vev = linspace(vi,vf,N);
[Uev,Vev] = meshgrid(uev,uev);
Results.Solution = (1-Uev).*(1-Vev).*sin(5.*Vev);
Results.Error = abs(Results.Solution - Results.Estimation);
Results.MaxError = max(Results.Error, [], 'all');
Results.MeanError = mean(Results.Error, 'all');
Results.StdError = std(Results.Error, [], 'all');
save('L1DMV_1_Spline', '-struct', 'sp');
save('L1DMV_1_Results', '-struct', 'Results');