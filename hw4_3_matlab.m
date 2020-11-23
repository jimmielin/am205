xmesh = linspace(-1, 1, 101);
solinit = bvpinit(xmesh, @guess);

sol = bvp4c(@bvpfcn, @bcfcn, solinit);
plot(sol.x, sol.y);

function dydx = bvpfcn(x,y)
    dydx = zeros(2,1);
    dydx = [y(2)
            exp(y(1))];
end

function res = bcfcn(ya,yb)
    res = [ya(1)
           yb(1)];
end

function g = guess(x)
    g = [0
         0];
end