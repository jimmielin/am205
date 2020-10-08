% Lx = b lower triangular solve in B space

%L = [1 0 0 0; 0 1 0 0;  1 1 1 0; 1 0 1 1];
%U = [1 0 1 0; 0 1 1 1;  0 0 1 0; 0 0 0 1];
%b = [1 1 1 0];

%A = [1 1 0 1; 1 1 0 1; 0 1 0 0; 0 1 1 1];


%a
[P, L, U] = LUpivotB(a);

% perform a test to make sure PA = LU
isequal(bin_mul(P,a), bin_mul(L,U))

% simply solve lower-triangular Ly=Pb first
pb = bin_mul(P, b);
y = fsolve(L, transpose(pb));

% now solve Ux=y
x = transpose(rsolve(U, y));

% now test
isequal(bin_mul(a, x), b)

% and output
x


%fsolve(L, b)
%rsolve(U, b)

% auxiliary functions for B-set
function s = Bsum(a, b)
    if (a ~= b); s = 1; else; s = 0; end
end

% LU decomposition with pivoting for B-set
function [P, L, U] = LUpivotB(A)
    [m, n] = size(A);
    L = eye(n);
    P = eye(n);
    U = A; % start with L, P being identity and U = A
    for k = 1:(m-1) % for each row
        pivot = max(abs(U(k:m, k))); % can just find 1 if needed
        for j = k:m
            if(abs(U(j,k)) == pivot)
                idx = j;
                break;
            end
        end
        U([k,idx],k:m) = U([idx,k],k:m);
        L([k,idx],1:(k-1)) = L([idx,k],1:(k-1));
        P([k,idx],:) = P([idx,k],:);
        % k, U, L, P
        for j = (k+1):m
            % note binary division
            if U(k,k) == 0
                j, k, idx
                U
                error("division by zero")
            end
            L(j,k) = U(j,k); % division is always equal /U(k,k)
            for l = k:m
                U(j,l) = Bsum(U(j,l), L(j,k)*U(k,l));
            end
        end
    end
end

% solve Lx=b
function x = fsolve(L, b)
    if size(L, 1) ~= size(b, 2) || size(b, 1) ~= 1
        error('matrix sizes incorrect')
    end
    
    x = zeros(1, size(L, 1));
    for i = 1:size(L, 1)
        % for each row do forward substitution
        num = b(i);
        for j = 1:(i-1)
            % num = num - L(i,j)*x(j);
            num = Bsum(num, L(i,j)*x(j));
        end
        den = L(i,i);
        % now perform division
        if den == 0
            error('division by zero')
        end
        % x(i) = num/den;
        x(i) = num; % divide by 1 is always 1
    end
end

% solve Ux=b
function x = rsolve(U, b)
    % upper triangular reverse substitution
    if size(U, 1) ~= size(b, 2) || size(b, 1) ~= 1
        error('matrix sizes incorrect')
    end
    
    x = zeros(1, size(U, 1));
    for i = size(U, 1):-1:1
        num = b(i);
        den = U(i,i);
        for j = i+1:size(U,1)
            num = Bsum(num, U(i,j)*x(j));
        end
        if den == 0
            error('division by zero')
        end
        x(i) = num; % divide by 1 is always 1
    end
end

% this is from the hw description
function A=bin_mul(L,U)
    [lnr,lnz]=size(L);
    [unr,unz]=size(U);
    if lnz ~= unr
        error('matrix dimension dismatch')
    end

    A=zeros(lnr,unz);

    for i=1:lnr
        for j=1:unz
            for k=1:lnz
                A(i,j)=bitxor(A(i,j),bitand(L(i,k),U(k,j)));
            end
        end
    end
end