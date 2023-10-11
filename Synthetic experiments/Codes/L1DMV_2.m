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
colmatv_ic1 = spcol(knotsu,k,[ui,uf],'sparse');
colmatuv_ic1 = kron(colmatu_ic1,colmatv_ic1);
colmatu_ic2 = spcol(knotsu,k,[ui,uf],'sparse');
colmatv_ic2 = spcol(knotsu,k,tauv,'sparse');
colmatuv_ic2 = kron(colmatu_ic2,colmatv_ic2);
CI = [zeros(size(colmatuv_ic1,1),1);zeros(size(colmatuv_ic2,1),1)];
B1 = sin(pi*Ugr).*sin(pi*Vgr);
B = [CI;B1];
A = [colmatuv_ic1;colmatuv_ic2;colmatu2v0+colmatu0v2];
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
z_u2 = diff(z,u,2); z_v2 = diff(z,v,2);
EQ = syms + z_v2 == sin(pi^u)*sin(pi^u);
ui = 0; uf = 1; vi = 0; vf = 1; 
CI = 1; CId=0;
solp = -sin(pi*u)*sin(pi*v)/2/pi^2;
N = 1e3; uev = linspace(ui,uf,N); vev = linspace(vi,vf,N);
[Uev,Vev] = meshgrid(uev,uev);
h1 = figure; hold on; grid on;
fsurf(solp,[ui,uf,vi,vf]);
title('Solution')
%% Results
N = 50; uev = linspace(ui,uf,N); vev = linspace(vi,vf,N);
[Uev,Vev] = meshgrid(uev,uev);
Results.Solution = -(sin(pi.*Uev).*sin(pi.*Vev))./2/pi^2;
Results.Error = abs(Results.Solution - Results.Estimation);
Results.MaxError = max(Results.Error, [], 'all');
Results.MeanError = mean(Results.Error, 'all');
Results.StdError = std(Results.Error, [], 'all');
save('L1DMV_2_Spline', '-struct', 'sp');
save('L1DMV_2_Results', '-struct', 'Results');