% lu factorization

%% time measurement
if false
    ns = 10:10:1200;
    times = [];

    for i = 1:length(ns)
        tic
        generate_g(ns(i));
        ts = toc;
        fprintf("n = %d, ts = %f\n", ns(i), ts);
        times = [times, ts];
    end
end

%% plotting results
plot(ns, times)
xlabel("n");
ylabel("time");


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