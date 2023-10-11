clc; close all; clear all;
Tumor = pcread('Tumour_1.ply');
nbrks = 9; k = 6;
x = Tumor.XLimits; y = Tumor.YLimits; z = Tumor.ZLimits;
xx = Tumor.Location(:, 1) - mean(x); yy = Tumor.Location(:, 2) - mean(y); zz = Tumor.Location(:, 3) - mean(z);
[azimuth,elevation,r] = cart2sph(xx,yy,zz); [atau Ia] = sort(azimuth); [etau Ie] = sort(elevation);
abreaks = [atau(1); atau(1:end-1)+diff(atau); atau(end)]; abreaks = [abreaks(1:(ceil(length(abreaks)/(nbrks-1))):end-1); abreaks(end)];
ebreaks = [etau(1); etau(1:end-1)+diff(etau); etau(end)]; ebreaks = [ebreaks(1:(ceil(length(ebreaks)/(nbrks-1))):end-1); ebreaks(end)];
aknots = augknt(abreaks, k); eknots = augknt(ebreaks, k);
colmatas = spcol(aknots, k, atau); colmata(Ia, :) = colmatas; colmata = repmat(colmata, [1 nbrks+4]);
colmates = spcol(eknots, k, etau); colmate(Ie, :) = colmates; colmate = kron(colmate, ones(1, nbrks+4));
A = [colmata.*colmate];
coef = A\r; coef = reshape(coef, [nbrks+4 nbrks+4]);
sp = spmak({aknots, eknots}, coef);
outcome = fnval(sp, [azimuth elevation].');
error1 = abs((r-outcome.')./r);
[mean(error1) std(error1) max(error1)]
Results.DomainAzimuth = azimuth;
Results.DomainElevation = elevation;
Results.Estimation = fnval(sp, [azimuth elevation].');
Results.Solution = r;
Results.Error = abs((r-outcome.')./r);
Results.MaxError = max(error1);
Results.StdError = std(error1);
Results.MeanError =mean(error1);
save('Tumour_1_Spline', '-struct', 'sp');
save('Tumour_1_Results', '-struct', 'Results');
[outx outy outz] = sph2cart(azimuth, elevation, outcome.'); out = [outx+mean(x) outy+mean(y) outz+mean(z)];
ptCloud = pointCloud(out);
figure;
pcshow(Tumor, colormap='w')
title('Tumour 1 pre-operative model')
set(gcf,'color','w');
figure;
pcshow(ptCloud)
title('Tumour 1 Spline outcome')
Representation = stlread('Tumour_1.stl');
T = Representation.ConnectivityList;
Rx = Representation.Points(:,1);
Ry = Representation.Points(:,2);
Rz = Representation.Points(:,3);
figure;
trimesh(T, Rx, Ry, Rz)
hold on;
scatter3(out(:,1), out(:,2), out(:,3), 'r')
legend('Ground truth mesh', 'Obtained points')
title('Tumour 1')
hold off;