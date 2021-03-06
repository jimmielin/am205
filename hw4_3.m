n = 101; % data: U0, U1, U2, ..., U100
h = 2/(n-1);

% the below arrays exclude padding
U_prev = zeros(n-2, 1);
U_next = zeros(n-2, 1);
J = zeros(n-2, n-2);
F = zeros(n-2, 1);

% construct the Jacobian... and solve with Newton
% fill the trivial terms first
for i = 1:(n-2)
    if i > 1
        J(i,i-1) = 1/h^2;
    end
    
    if i < (n-2)
        J(i,i+1) = 1/h^2;
    end
end

for iter = 1:500 % max iterations...
    % construct jac and fn
    for i = 1:(n-2)
        J(i,i) = -2/h^2 - exp(U_next(i));
        
        if i == 1
            F(i) = (U_next(2) - 2*U_next(1))/h^2 - exp(U_next(1));
        elseif i == (n-2)
            F(i) = (-2*U_next(n-2) + U_next(n-3))/h^2 - exp(U_next(n-2));
        else
            F(i) = (U_next(i+ 1)-2*U_next(i)+U_next(i-1))/h^2 - exp(U_next(i));
        end
    end
    
    % copy data
    U_prev(:) = U_next(:);
    
    % do not invert the matrix... it is really easy to get into numerical
    % instability mess, just solve the unknown
    % J(U_next - U_prev) = -F(U_prev)
    F = -F;
    Usolve = J \ F;
    U_next = Usolve + U_prev;
    
    % are we there yet?
    DeltaU = U_next - U_prev;
    relstep = norm(DeltaU)/norm(U_next);
    
    fprintf("iter = %d, relstep = %f\n", iter, relstep);
    
    if relstep <= 1e-10
        break
    end
end

%% plot with approximation
hs = arrayfun(@(i) -1+i*h, 0:(n-1));
vs = [0; U_next; 0];

plot(hs, vs);