% params:

xs = 0:0.001:3;
ys = arrayfun(@(t) Sx(t), xs);
zs = arrayfun(@(t) 2/sqrt(3)*sin(2*pi/3*t), xs);

plot(xs, ys, 'LineWidth', 1)
hold on;
plot(xs, zs, 'LineWidth', 1)
legend('s_x(t)', 'sine truth')

% (c)
hold off;
ms = arrayfun(@(t) Sy(t), xs);
ns = arrayfun(@(t) 2*cos(2*pi*t/3), xs);
plot(xs, ms, 'LineWidth', 1);
hold on;
plot(xs, ns, 'LineWidth', 1);
legend('s_y(t)', 'cosine truth');

% (d) parametric curve
hold off;
Xt = arrayfun(@(v) sqrt(3)/2*v, zs);
Yt = arrayfun(@(v) v/2, ns);
X = arrayfun(@(v) sqrt(3)/2*v, ys);
Y = arrayfun(@(v) v/2, ms);
plot(X, Y, 'LineWidth', 1);
legend('sqrt(3)/2*s_x(t) and 1/2*s_y(t) parametric');
A = polyarea(X, Y); % 
fprintf(1, 'A = %.10f\n', A) 

% (e) parameteric curve
ts = 0:0.001:3;
g = 1/(Sy(0) - Sy(3/2));
wtx = arrayfun(@(t) g*sqrt(3)*(2/sqrt(3)*sin(2*pi/3*t)-2/sqrt(3)*sin(2*pi/3*(t-3/2))), ts);
wty = arrayfun(@(t) g*(2*cos(2*pi*t/3)-2*cos(2*pi*(t-3/2)/3)), ts);
wx = arrayfun(@(t) g*sqrt(3)*(Sx(t)-Sx(t-3/2)), ts);
wy = arrayfun(@(t) g*(Sy(t)-Sy(t-3/2)), ts);
%plot(wx, wy, 'LineWidth', 1);
%legend('w_x(t), w_y(t) parameteric');
%A2 = polyarea(wx, wy);
%fprintf(1, 'A2 = %.10f\n', A2)


function S = Sx(t)

    % patch periodic interval with a kludge
    while t < 0
        t = t + 3;
    end
    
    while t > 3
        t = t - 3;
    end

    A = [0 1 -1];
    B = [2 -1 -1];
    C = [0 -3 3];
    D = [-1 2 -1];
    
    if t >= 0 && t < 1
        idx = 1;
    elseif t < 2
        idx = 2;
    elseif t <= 3
        idx = 3;
    end
    
    S = A(idx) + B(idx)*(t-idx+1) + ...
        C(idx)*(t-idx+1)^2 + ...
        D(idx)*(t-idx+1)^3;
end


function S = Sy(t)
    % patch periodic interval with a kludge
    while t < 0
        t = t + 3;
    end
    
    while t > 3
        t = t - 3;
    end
    A = [2 -1 -1];
    B = [0 -3 3];
    C = [-6 3 3];
    D = [3 0 -3];
    
    if t >= 0 && t < 1
        idx = 1;
    elseif t < 2
        idx = 2;
    elseif t <= 3
        idx = 3;
    end
    
    S = A(idx) + B(idx)*(t-idx+1) + ...
        C(idx)*(t-idx+1)^2 + ...
        D(idx)*(t-idx+1)^3;
end