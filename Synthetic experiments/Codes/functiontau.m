function y = functiontau(breaks, points)
    syms x
    subinterval = diff(breaks);
    roots = vpasolve(legendreP(points,x) == 0);
    tau = (breaks(2:end).*ones(points, 1) + breaks(1:end-1).*ones(points, 1) + subinterval.*roots)/2;
    tau = double(tau);
    tau = reshape(tau, 1, []);
    y = sort([breaks(1) tau breaks(end)]);
end