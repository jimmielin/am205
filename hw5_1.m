%% gradient descent ...
% steepest descent line search

x0 = [2, 1];
lxs = [x0(1)];
lys = [x0(2)];
lzs = [f(x0(1), x0(2))];
for i = 1:2000%2000
    % maximum of 2000 iterations ...
    % generate grad
    ngrad = (-1) * gradF(x0(1), x0(2));
    fn = f(x0(1), x0(2));
    
    % optimize for stepsize ... line search
    fn_mini = @(a) f(x0(1)+a*ngrad(1), x0(2)+a*ngrad(2));
    a = fminsearch(fn_mini, 0);
    
    % update x0
    x0 = [x0(1)+a*ngrad(1), x0(2)+a*ngrad(2)];
    fn = f(x0(1), x0(2));
    
    lxs = [lxs x0(1)];
    lys = [lys x0(2)]; 
    lzs = [lzs fn];
    
    fprintf("n = %d, x1, x2 = %f, %f, fnw = %f\n", i, x0(1), x0(2), fn);

    % check abs step size ...
    if i > 1 && (abs(x0(1) - lxs(i-1)) < 1e-8 || abs(x0(2) - lys(i-1)) < 1e-8)
        break
    end
end

%% Newton's method
x0 = [2, 1];
lxs = [x0(1)];
lys = [x0(2)];
lzs = [f(x0(1), x0(2))];
for i = 1:15%2000
    % maximum of 2000 iterations ...
    % generate grad
    ngrad = (-1) * gradF(x0(1), x0(2));
    fn = f(x0(1), x0(2));
    
    % hessian
    x = x0(1);
    y = x0(2);
    hx = [-400*(y-x^2)+800*x^2+2 -400*x; -400*x 200];
    
    % now linear solve for vector s...
    s = hx \ transpose(ngrad);
    
    % update x0
    x0 = [x0(1)+s(1), x0(2)+s(2)];
    fn = f(x0(1), x0(2));
    
    lxs = [lxs x0(1)];
    lys = [lys x0(2)]; 
    lzs = [lzs fn];
    
    fprintf("n = %d, x1, x2 = %f, %f, fnw = %f\n", i, x0(1), x0(2), fn);

    % check abs step size ...
    if i > 1 && (abs(x0(1) - lxs(i-1)) < 1e-8 || abs(x0(2) - lys(i-1)) < 1e-8)
        break
    end
end

%% now plot ...
opengl hardware;
[xs, ys] = peaks(30);
fs = 100*(ys-xs.^2)^2 + (1-xs).^2;
s = surf(xs, ys, fs, 'FaceAlpha',0.5);
s.EdgeColor = 'none';
title("Surface plot of Rosenbrock Function, starting [2, 1]");
colorbar

hold on;
plot3(lxs, lys, lzs, 'Color', '#FF0000', 'LineWidth', 2);

function f = f(x, y)
    f = 100*(y-x^2)^2 + (1-x)^2;
end

function gradF = gradF(x, y)
    gradF = [2*(200*x^3 - 200*x*y + x - 1), 200*(y-x^2)];
end