% integrate using composite trapezoid rule
% not generalizing

%% (a)
% Run the integration for multiple n...
nmax = 50;
Intt = zeros(1, nmax);
abser = zeros(1, nmax);
for n = 1:nmax
    Intt(n) = In(n);
    abser(n) = abs(Intt(n) - 8*pi/9); % take abs. val
end

hs = arrayfun(@(n) pi/3/n, 1:nmax);

loglog(hs, abser+0.0000000001); % make a slight transformation for zero-data
xlabel("h");
ylabel("absolute error");
xlim([min(hs) 1.05]);
ylim([10^(-6) 2]);
grid on
title("Absolute Error of I_n f(x) over [0, \pi/3]")

hold on
% find the maximum infinity norm of f''
% which is ~ 0.69325 pos and ~ -16 at 0, so its 16
bnd = arrayfun(@(n) (pi/3/n)^2*pi/36*16, 1:nmax);
loglog(hs, bnd)
legend(["Abs. Err", "E(h)"], "location", "northwest")


% function Intg = In(n)
%     Intg = 0.0;
%     % n = 50;
%     ll = 0;
%     hh = pi/3;
% 
%     h = (hh-ll)/n;
%     for i = 1:(n-1)
%         Intg = Intg + 2 * f(ll + i*h);
%     end
%     Intg = Intg + f(ll) + f(hh);
%     Intg = Intg * h / 2;
% end
% 
% function f = f(x)
%    f = 1/(5/4 - cos(x)); 
% end

%% (b) 
nmax = 50;
Intt = zeros(1, nmax);
abser = zeros(1, nmax);
for n = 1:nmax
    Intt(n) = In(n);
    abser(n) = abs(Intt(n) - 8*pi/3); % take abs. val
end

hs = arrayfun(@(n) 2*pi/n, 1:nmax);

loglog(hs, abser+0.0000000001); % make a slight transformation for zero-data
xlabel("h");
ylabel("absolute error");
xlim([min(hs) 6.5]);
grid on
title("Absolute Error of I_n f(x) over [0, 2\pi]")

function Intg = In(n)
    Intg = 0.0;
    % n = 50;
    ll = 0;
    hh = 2*pi;

    h = (hh-ll)/n;
    for i = 1:(n-1)
        Intg = Intg + 2 * f(ll + i*h);
    end
    Intg = Intg + f(ll) + f(hh);
    Intg = Intg * h / 2;
end

function f = f(x)
   f = 1/(5/4 - cos(x)); 
end
