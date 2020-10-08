n = 4;
xs = arrayfun(@(j) cos((2*j-1)*pi/(2*n)), 1:n);
ys = arrayfun(@(x) x^2*exp(x), xs);

% lagrange interp
% where g(x) = g0 + g1 * x + g2 * x^2 ... or
% linear combination of Lagrange Basis Polys

% g(x) = Sum gmi * li(n)
% li(n) := Product(m /= i) (x - x_m)/(x_i - x_m)

% use symbolic math
syms x;

g_tempsum = str2sym('0');
for i = 1:length(xs)
    li_tempprod = str2sym('1');
    for j = 1:length(xs)
        if i == j
            continue
        end
        li_tempprod = li_tempprod * ...
            str2sym(strcat('(x-', num2str(xs(j)), ')/(', num2str(xs(i)-xs(j)), ')'));
    end
    g_tempsum = g_tempsum + ys(i) * li_tempprod;
    li_tempprod;
end

poly = expand(g_tempsum)

% plotting code is a little more intricate
pxs = -1:0.002:1;
pys = arrayfun(@(x) x^2*exp(x), pxs);
pps = arrayfun(@(xv) subs(poly,x,xv), pxs);
plot(pxs, pys);
hold on;
plot(pxs, pps);
legend('f(x)', 'p_3(x)');


%(b) Calculate infinity norm
% do element wise subtraction then a pass of abs
ns = arrayfun(@(i) abs(pys(i)-pps(i)), 1:length(pxs));
max(ns)



% try spline interpolation
sys = spline(xs,ys,pxs);
hold off;
plot(pxs, pys);
hold on;
plot(pxs, sys);
legend('f(x)', 'p_3^+(x)');
nns = arrayfun(@(i) abs(sys(i)-pys(i)), 1:length(pxs));
max(nns)