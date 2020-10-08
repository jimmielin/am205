n = [1 2 3 4 5];
gm = [1 1 2 6 24];
hm = arrayfun(@(g) log(g), gm);

% lagrange interp
% where g(x) = g0 + g1 * x + g2 * x^2 ... or
% linear combination of Lagrange Basis Polys

% g(x) = Sum gmi * li(n)
% li(n) := Product(m /= i) (x - x_m)/(x_i - x_m)

% use symbolic math
syms x;

g_tempsum = str2sym('0');
for i = 1:length(n)
    li_tempprod = str2sym('1');
    for j = 1:length(n)
        if i == j
            continue
        end
        li_tempprod = li_tempprod * ...
            str2sym(strcat('(x-', num2str(n(j)), ')/(', num2str(n(i)-n(j)), ')'));
    end
    g_tempsum = g_tempsum + gm(i) * li_tempprod;
    li_tempprod;
end

g = expand(g_tempsum)

% now do the same but for hm
tempsum = str2sym('0');
for i = 1:length(n)
    li_tempprod = str2sym('1');
    for j = 1:length(n)
        if i == j
            continue
        end
        li_tempprod = li_tempprod * ...
            str2sym(strcat('(x-', num2str(n(j)), ')/(', num2str(n(i)-n(j)), ')'));
    end
    tempsum = tempsum + hm(i) * li_tempprod;
    li_tempprod;
end

p = expand(tempsum)
h = exp(p);

xs = 1:0.01:5;
gs = arrayfun(@(xp) subs(g,x,xp), xs);
hs = arrayfun(@(xp) subs(h,x,xp), xs);
gms = arrayfun(@(xp) gamma(xp), xs);
plot(xs, gms, 'LineWidth', 1);
hold on;
plot(xs, gs, 'LineWidth', 1);
plot(xs, hs, 'LineWidth', 1);
xlabel('n', 'Fontsize', 15);
ylabel('y', 'Fontsize', 15);
set(gca,'Fontsize', 15);

legend('gamma(x)', 'g(x)', 'h(x)')

% compute relative errors and their max
gs_error_max = 0;
hs_error_max = 0;
for i = 1:length(xs)
    gs_error = abs(gs(i) - gms(i))/gms(i);
    hs_error = abs(hs(i) - gms(i))/gms(i);
    if gs_error > gs_error_max
        gs_error_max = gs_error;
    end
    if hs_error > hs_error_max
        hs_error_max = hs_error;
    end
end

gs_error_max
hs_error_max