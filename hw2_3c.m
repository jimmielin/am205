% generate matrix A first
% handy conversion from i, j to numeric: num = 7*i+j+1
% from num to i,j, i = idivide(num-1, 7); j = mod(num-1, 7);

kerSize = zeros(9,9);

for m=1:9
    for n=1:9
        %m = 8; % vertical (rows)
        %n = 4; % horizontal (cols)
        s = m*n;
        A = zeros(s,s);
        for num = 1:s
            i = fix((num-1)/n);
            j = mod(num-1, n);
            % i, j in [0, 6]
            % affecting indices are i+/-1, j+/- total 4 and self
            % update matrix vectors as
            % b(num) = \sum^m [A(num, m) x(m)], m is converted from
            % affecting indices

            % just some copy paste...
            r = n*i+j+1; A(num,r) = 1;
            % fprintf("i, j, num, r = %d %d %d %d\n", i, j, num, r)
            if (i ~= 0); r = n*(i-1)+j+1; A(num,r) = 1; end
            if (j ~= 0); r = n*i+j; A(num,r) = 1; end
            if (i ~= (m-1)); r = n*(i+1)+j+1; A(num,r) = 1; end
            if (j ~= (n-1)); r = n*i+j+2; A(num,r) = 1; end
        end

        ks = size(null2(A));
        k1 = ks(1);
        k2 = ks(2);
        fprintf("n = %d, m = %d, kersizes = %d, %d\n", n, m, k1, k2)
        kerSize(n,m) = size(null2(A), 2);
    end
end

kerSize


%Computes the null space of a binary matrix A over GF(2)
function [NullSpace]=null2(A)
    A=mod(A,2);
    %number of constraints:
    m=size(A,1);
    %number of variables:
    n=size(A,2);
    %number of independent constraints:
    r=gfrank(A,2);
    %Take care of the trivial cases:
    if (r==n)
        NullSpace=[];
    elseif (r==0)
        NullSpace=eye(n,n);
    end
    %Add one constraint at a time.
    %Maintain a matrix X whose columns obey all constraints examined so far.
    %Initially there are no constraints:
    X=eye(n,n);
    for i=1:m
        y=mod(A(i,:)*X,2);
        % identify 'bad' columns of X which are not orthogonal to y
        % and 'good' columns of X which are orthogonal to y
        GOOD=[X(:,find(not(y)))];
        %convert bad columns to good columns by taking pairwise sums
        if (nnz(y)>1)
          BAD=[X(:,find(y))];
          BAD=mod(BAD+circshift(BAD,1,2),2);
          BAD(:,1)=[];
        else
            BAD=[];
        end
        X=[GOOD,BAD];
    end%for i
    NullSpace=X;
end