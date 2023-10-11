clc; close all; clear all;
Data = load('workspace_IMUs_sqt1.mat');
t = Data.Lshin_array_c(:,1)';
gyr = Data.Lshin_array_c(:,2:4)'/180*pi;
gt = Data.Lshin_rot_eul(:,[3,2,1]);
k = 6; nPBB = 7;
brks = diff(t)./2;
brks = [t(1)-brks(1) t(1:end-1)+brks t(end)+brks(end)];
brks = [brks(1) brks(nPBB:nPBB:end-1) brks(end)];
knots = augknt(brks, k);
colmat = spcol(knots, k, brk2knt(t, 3));
[rows ncoefs] = size(colmat); iv=10000;
omega1 = gyr(1, :).'; omega2 = gyr(2, :).'; omega3 = gyr(3, :).';
A = [colmat(1, :); colmat(2:3:end, :)];
coef = A\[0 0 0; gyr.'];
sp = spmak(knots, coef.');
spd1 = fnder(sp, 1);
outcome = fnval(sp, t);
outcomed1 = fnval(spd1, t);
for i=1:5
    yaw = outcome(1, :).'; yawd1 = outcomed1(1, :).';
    pitch = outcome(2, :).'; pitchd1 = outcomed1(2, :).';
    roll = outcome(3, :).'; rolld1 = outcomed1(3, :).';
    A1 = [iv.*colmat(1, :); colmat(2:3:end, :)];
    A2 = [zeros(1, ncoefs); -cos(pitch).*rolld1.*colmat(1:3:end, :)];
    A3 = [zeros(1, ncoefs); -sin(pitch).*colmat(2:3:end, :)];
    A4 = [zeros(1, ncoefs); (cos(pitch).*rolld1.*cos(yaw) - sin(yaw).*pitchd1).*colmat(1:3:end, :)];
    A5 = [iv.*colmat(1, :); cos(yaw).*colmat(2:3:end, :) - sin(yaw).*rolld1.*sin(pitch).*colmat(1:3:end, :)];
    A6 = [zeros(1, ncoefs); sin(yaw).*cos(pitch).*colmat(2:3:end, :)];
    A7 = [zeros(1, ncoefs); -(cos(yaw).*pitchd1 + sin(yaw).*cos(pitch).*rolld1).*colmat(1:3:end, :)];
    A8 = [zeros(1, ncoefs); -cos(yaw).*sin(pitch).*rolld1.*colmat(1:3:end, :)-sin(yaw).*colmat(2:3:end, :)];
    A9 = [iv.*colmat(1, :); cos(yaw).*cos(pitch).*colmat(2:3:end, :)];
    B1 = [0; omega1 - cos(pitch).*rolld1.*pitch];
    B2 = [0; omega2 + (cos(pitch).*rolld1.*cos(yaw) - sin(yaw).*pitchd1).*yaw - sin(yaw).*rolld1.*sin(pitch).*pitch];
    B3 = [0; omega3 - (cos(yaw).*pitchd1 + sin(yaw).*cos(pitch).*rolld1).*yaw - cos(yaw).*sin(pitch).*rolld1.*pitch];
    A = [A1 A2 A3; A4 A5 A6; A7 A8 A9];
    B = [B1; B2; B3];
    coeffs = A\B;
    coefs = zeros(3, ncoefs);
    coefs(1, :) = coeffs(1:ncoefs);
    coefs(2, :) = coeffs(1+ncoefs:2*ncoefs);
    coefs(3, :) = coeffs(1+2*ncoefs:end);
    convergence = [abs(coefs(1, :) - sp.coefs(1, :)); abs(coefs(2, :) - sp.coefs(2, :)); abs(coefs(3, :) - sp.coefs(3, :))];
    [mean(convergence(1, :)) mean(convergence(2, :)) mean(convergence(3, :))]
    [max(convergence(1, :)) max(convergence(2, :)) max(convergence(3, :))]
    sp = spmak(knots, coefs);
    spd1 = fnder(sp, 1);
    outcome = fnval(sp, t);
    outcomed1 = fnval(spd1, t);
end
save('Sqt_orientation_Spline', '-struct', 'sp');
Results.t = t;
Results.yawgt = gt(:, 1).';
Results.pitchgt = gt(:, 2).';
Results.rollgt = gt(:, 3).';
Results.yaw = outcome(1, :);
Results.pitch = outcome(2, :);
Results.roll = outcome(3, :);
Results.erroryaw = abs(gt(:, 1).' - outcome(1, :));
Results.errorpitch = abs(gt(:, 2).' - outcome(2, :));
Results.errorroll = abs(gt(:, 3).' - outcome(3, :));
Results.erroryawmean = mean(Results.erroryaw);
Results.errorpitchmean = mean(Results.errorpitch);
Results.errorrollmean = mean(Results.errorroll);
Results.erroryawmax = max(Results.erroryaw);
Results.errorpitchmax = max(Results.errorpitch);
Results.errorrollmax = max(Results.errorroll);
Results.erroryawstd = std(Results.erroryaw);
Results.errorpitchstd = std(Results.errorpitch);
Results.errorrollstd = std(Results.errorroll);
figure; hold on;
plot(t, gt(:, 1), 'r');
plot(t, outcome(1, :), 'b--');
legend('Ground truth', 'Approximation')
xlabel('time (s)')
ylabel('rotation (rad)')
title('Yaw')
hold off;
figure; hold on;
plot(t, gt(:, 2), 'r');
plot(t, outcome(2, :), 'b--');
legend('Ground truth', 'Approximation')
xlabel('time (s)')
ylabel('rotation (rad)')
title('Pitch')
hold off;
figure; hold on;
plot(t, gt(:, 3), 'r');
plot(t, outcome(3, :), 'b--');
legend('Ground truth', 'Approximation')
xlabel('time (s)')
ylabel('rotation (rad)')
title('Roll')
hold off;
vector = [outcome(1, :); outcome(2, :); outcome(3, :)];
rotZYX = eul2rotm(vector.');
acc = Data.Lshin_array_c(:, 5:7).';
acc1 = zeros(size(acc));
for i=1:length(vector)
    acc1(:, i) = acc(:, i)'*rotZYX(:, :, i).';
end
gravity = 9.8;
acc1 = gravity.*acc1;
A = [colmat(1, :); colmat(2, :); colmat(3:3:end, :)];
C = [[A; zeros(size(A)); zeros(size(A))], [zeros(size(A)); A; zeros(size(A))], [zeros(size(A)); zeros(size(A)); A]];
coef1 = ridge([0; 0; acc1(1,:).'], A, 0.5, 0);
coef2 = ridge([0; 0; acc1(2,:).'], A, 5, 0);
coef3 = ridge([0; 0; acc1(3,:).'], A, 50, 0);
coef = [coef1(2:end, :)-std(coef1(2:end, :)) coef2(2:end, :)-std(coef2(2:end, :)) coef3(2:end, :)-std(coef3(2:end, :))];
sp = spmak(knots, coef.');
spd1 = fnder(sp, 1);
outcome = fnval(sp, t);
outcomed1 = fnval(spd1, t);
gt = Data.Lshin_pos;
Results.zgt = gt(:, 1).';
Results.ygt = gt(:, 2).';
Results.xgt = gt(:, 3).';
Results.z = outcome(1, :);
Results.y = outcome(2, :);
Results.x = outcome(3, :);
Results.errorz = abs(gt(:, 1).' - outcome(1, :));
Results.errory = abs(gt(:, 2).' - outcome(2, :));
Results.errorx = abs(gt(:, 3).' - outcome(3, :));
Results.errorzmean = mean(Results.errorz);
Results.errorymean = mean(Results.errory);
Results.errorxmean = mean(Results.errorx);
Results.errorzmax = max(Results.errorz);
Results.errorymax = max(Results.errory);
Results.errorxmax = max(Results.errorx);
Results.errorzstd = std(Results.errorz);
Results.errorystd = std(Results.errory);
Results.errorxstd = std(Results.errorx);
save('Sqt_pos_Spline', '-struct', 'sp');
save('Sqt_Results', '-struct', 'Results');
figure; hold on;
plot(t, gt(:, 1), 'r');
plot(t, outcome(1, :), 'b--');
legend('Ground truth', 'Approximation')
xlabel('time (s)')
ylabel('distance (m)')
title('Z axis')
hold off;
figure; hold on;
plot(t, gt(:, 2), 'r');
plot(t, outcome(2, :), 'b--');
legend('Ground truth', 'Approximation')
xlabel('time (s)')
ylabel('distance (m)')
title('Y axis')
hold off;
figure; hold on;
plot(t, gt(:, 3), 'r');
plot(t, outcome(3, :), 'b--');
legend('Ground truth', 'Approximation')
xlabel('time (s)')
ylabel('distance (m)')
title('X axis')
hold off;