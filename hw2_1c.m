% solve generalized
% (4x-y)^4+x^4-1=0
% x^4+y^4-1=0

% newton method
% number of iterations n:
n = 16;
f = @(d)[d(1)^4+d(2)^4-1; (4*d(1)-d(2))^4+d(1)^4-1];
df = @(d)[4*d(1)^3, 4*d(2)^3; 16*(4*d(1)-d(2))^3+4*d(1)^3, 4*(d(2)-4*d(1))^3];

% make an initial guess
d = [-1; 1];
for i = 1:n
    dd = -df(d)\f(d); % solve for the delta...
    d = d + dd;
end
disp("End result of f(d) (should be close to zero):")
f(d)
disp("End d:")
d

% make implicit plot (=0) for 4-norm
fimplicit(@(x,y) (x.^4 + y.^4).^(1/4) - 1);
hold on
fimplicit(@(x,y) ((4*x-y).^4 + x.^4).^(1/4) - 1);
legend(["||b||_4=1" "||Ab||_4=1"], 'Location', 'northwest')

% plot points
hold off;
xs = [0, 0, -1/sqrt(5), 1/sqrt(5), -1/2, 1/2, 1/17^(1/4), -1/17^(1/4)];
ys = [-1,1, -2/sqrt(5), 2/sqrt(5), -1, 1, 2/17^(1/4), -2/17^(1/4)];
scatter(xs, ys)
hold on;