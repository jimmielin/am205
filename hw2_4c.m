if false
    ns = 10:10:200;
    reles = [];

    for i = 1:length(ns)
        % this will be in a loop later
        n = ns(i);

        x = ones(n, 1);
        Gn = generate_g(n);
        b = Gn * x; % right-hand side vector (n x n * n x 1 = n x 1)

        % solve Gn xhat = b using LU factorization
        % (obviously the answer is Gn xhat = Gn x => xhat = x)
        % no need to implement, just use matlab \
        xhat = Gn \ b;

        % calculate 2-norm relative error
        rele = norm(x - xhat)/norm(x); % always /1

        fprintf("n = %d, rele = %f\n", n, rele);

        reles = [reles, rele];
    end

    plot(ns, reles);
    xlabel("n");
    ylabel("2-norm reelative error");
end

Gn = generate_g(8);


%% generate matrix
function Gn = generate_g(n)
    Gn = eye(n); % first generate eye
    % rightmost column is 1
    Gn(:,n) = 1;
    % fill lower triangular with -1
    for i = 2:n
        Gn(i,1:(i-1)) = -1;
    end
end