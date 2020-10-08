% hw1-2d poly regression
n = 100;
xs = arrayfun(@(j) cos((2*j-1)*pi/(2*n)), 1:n);
ys = arrayfun(@(x) x^2*exp(x), xs);

% points are in xs, ys - fit to third deg poly
poly = polyfit(xs, ys, 3);

% function
syms x;
polyf = str2sym("1.2195*x^3+1.5403*x^2-0.0558*x-0.0692");

% compute infinity norm
pxs = -1:0.002:1;
pys = arrayfun(@(x) x^2*exp(x), pxs);
pps = arrayfun(@(xv) subs(polyf,x,xv), pxs);
plot(pxs, pys);
hold on;
plot(pxs, pps);
legend('f(x)', 'p_3^+(x)');

% do element wise subtraction then a pass of abs
ns = arrayfun(@(i) abs(pys(i)-pps(i)), 1:length(pxs));
max(ns)