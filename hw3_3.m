% adaptive integration using 3-pt gaussian primitive

% note: need to output # of intervals and total accumulated error estimate
% acc. error estimate is the sum of each interval

%% (a)
%fs(1, 2, 1/3)
%g(2, 0)
phi = 0.95;

xs = -0.5:0.001:0.5;
gs = arrayfun(@(x) g(x, phi), xs);

plot(xs, gs);
legend([strcat("g(x; ", num2str(phi), ")")]);

%% (b)
%result = int_adaptive(@f3, 0, 1, 0.000001);
%fprintf("Result = %f, intervals = %d, sumerror = %f\n", ...
%        result.value, result.intervals, result.error)

phis = 0.01:0.01:0.99;
is = zeros(1, length(phis));
itvs = zeros(1, length(phis));

for p = 1:length(phis)
    phi = phis(p);
    % g_lambda = @(x) g(x, phi);
    i_param = int_adaptive(@(x) g(x, phi), -0.5, 0.5, 0.000001);
    
    is(p) = i_param.value;
    itvs(p) = i_param.intervals;
end

plot(phis, is);
legend(["I(\phi)"]);
xlabel("\phi");

pause
plot(phis, itvs);
legend(["Intervals(\phi)"]);
xlabel("\phi");
ylabel("Integration intervals required");

% subroutine signature:
%  int_adaptive(@f, a, b, tol) 
%  => (for first call, set prev_error to arbitrarily large)
%  terminates when current_error < tol * (b-a)
%  current_error is a trial of 
%  returns: struct { value, intervals, error }

function iaout = int_adaptive(f, a, b, T)
    % try the primitive integration
    p1_val  = gauss3(f, a, b);
    p1_test = gauss3(f, a, (a+b)/2) + gauss3(f, (a+b)/2, b);
    
    % fprintf("p1_val = %f, p1_test = %f\n", p1_val, p1_test)
    err = abs(p1_test - p1_val);
    if err > T*(b-a)
        % fprintf("=> larger than threshold, branching out\n\n")
        
        % this would be so much easier with haskell foldl
        resLeft = int_adaptive(f, a, (a+b)/2, T);
        resRight = int_adaptive(f, (a+b)/2, b, T);
        out.value = resLeft.value + resRight.value;
        out.intervals = resLeft.intervals + resRight.intervals;
        out.error = resLeft.error + resRight.error;
    else
        out.value = p1_val;
        out.intervals = 1;
        out.error = err;
    end
    iaout = out;
end

% integrator primitive
function v = gauss3(f, a, b)
    k = (b-a)/2;
    m = (b+a)/2;
    v = 5/9*f(-k*sqrt(3/5)+m) + 8/9*f(m) + 5/9*f(k*sqrt(3/5)+m);
    v = k*v;
end

% can pass this function as function handle @f1
function y = g(x, phi)
    k = 0;
    y = 0.0;
    % determine phi^k target
    n_target = 64;
    for n = 1:64
        if phi^n < 10^-16
            n_target = n;
            break;
        end
    end
    
    y = fs(x, n_target, phi);
end

function y = fs(x, k, phi)
    if k == 0
        y = abs(x);
    else
        y = abs(fs(x, k-1, phi) - phi^k);
    end
end