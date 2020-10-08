%% read in every image as S_j
files = dir('am205_leaves_534x420/main/*.png');
n = 420;
m = 534;
S = [];
L = size(files, 1);
A = zeros(n*m*3, L);
Sjv_sum = zeros(n*m*3, 1);
ptr = 1;
for file = files'
    Sj = imread(strcat('am205_leaves_534x420/main/', file.name));
    % Sj is 280x356x3 unit8 (0~255) for H x W x RGB
    % compress into vector 3mn
    Sjv = zeros(n*m*3, 1);
    Sjv(1:n*m) = Sj(:,:,1);         % r
    Sjv(n*m+1:2*n*m) = Sj(:,:,2);   % g
    Sjv(n*m*2+1:3*n*m) = Sj(:,:,3); % b
    
    % convert to cmyk (0 - 255 range orig. new 0 - 1 range)
    Sjv(:) = (255.0 - Sjv(:)) ./ 255.0;
    
    % for later averaging
    Sjv_sum = Sjv_sum + Sjv;
    
    % write into A
    A(:,ptr) = Sjv;
    ptr = ptr + 1;
end

%% compute average leaf vector
Sjv_sum = Sjv_sum / L;

% rearrange into image format (temporary in Sj)
Sj = zeros(n, m, 3);  % handle as double so (1, 1, 1) makes sense in matlab
Sj(:,:,1) = reshape(Sjv_sum(1:n*m), [n, m]);
Sj(:,:,2) = reshape(Sjv_sum(n*m+1:2*n*m), [n, m]);
Sj(:,:,3) = reshape(Sjv_sum(n*m*2+1:3*n*m), [n, m]);

%% draw
Sj_rgb = 1 - Sj;
image(Sj_rgb)

%% (b) assemble matrix A (3mn x L) where each column A(:,j) = Sj - Savg
% Sjv_sum is now "helpfully" wrong named, it is actually the average...
for i = 1:L
    % previous step loop already wrote into A, now just need to sub the avg
    A(:,i) = A(:,i) - Sjv_sum;
end

% perform reduced singular value decomposition s.t.
[U,S,V] = svd(A, 'econ');

%% (b) locate and plot
% for j = 1,2,3
jmax = 3;
cj = zeros(jmax);
dj = zeros(jmax);
upj = zeros(jmax, size(A, 1));
unj = zeros(jmax, size(A, 1));
for j = 1:jmax
    cj = min(U(:,j)); dj = max(U(:,j));
    % scaled pos, neg components
    for i = 1:size(A, 1)
        upj(j,i) = max(0, U(i,j)/dj);
        unj(j,i) = max(0, U(i,j)/cj); % neg/neg =pos
    end
    
    % optionally plot
    if j <= 3
        % rearrange into image and plot
        Sj = zeros(n, m, 3); % handle as double
        Sj(:,:,1) = reshape(upj(j,1:n*m), [n, m]);
        Sj(:,:,2) = reshape(upj(j,n*m+1:2*n*m), [n, m]);
        Sj(:,:,3) = reshape(upj(j,n*m*2+1:3*n*m), [n, m]);
        Sj_rgb = 1 - Sj;
        hold off
        image(Sj_rgb)
        hold on
        title(strcat("u^P_{:,", num2str(j), "}"));
        pause
        
        Sj(:,:,1) = reshape(unj(j,1:n*m), [n, m]);
        Sj(:,:,2) = reshape(unj(j,n*m+1:2*n*m), [n, m]);
        Sj(:,:,3) = reshape(unj(j,n*m*2+1:3*n*m), [n, m]);
        Sj_rgb = 1 - Sj;
        hold off
        image(Sj_rgb)
        hold on
        title(strcat("u^N_{:,", num2str(j), "}"));
        pause
    end
end

%% (c) projection operator
% define image vector T
T = A(:, 33) + Sjv_sum; % random image, :,number
k = 16;
P = Sjv_sum; % this is the average vec starting point...
for j = 1:k
    P = P + transpose(U(:,j)) * (T - Sjv_sum) * U(:,j);
    j
end
% format and plot
Sj(:,:,1) = reshape(P(1:n*m), [n, m]);
Sj(:,:,2) = reshape(P(n*m+1:2*n*m), [n, m]);
Sj(:,:,3) = reshape(P(n*m*2+1:3*n*m), [n, m]);
Sj_rgb = 1 - Sj;
hold off
image(Sj_rgb)
hold on
title(strcat("P(T,k), Tidx = 32, k = ", num2str(k)));


%% (d) compute projections
% store projection for each image in new matrix...
Pjs = zeros(3*m*n, L); % each projection stored in a column vector in Pjs mtrx
for i = 1:L % note offset by 1 since we are 1-idxd
    % define image vector T
    T = A(:, i) + Sjv_sum; % random image, :,number
    k = 8;
    P = Sjv_sum; % this is the average vec starting point...
    for j = 1:k
        P = P + transpose(U(:,j)) * (T - Sjv_sum) * U(:,j);
    end
    % copy in
    Pjs(:,i) = P;
end

%% (d2) compute distances
dis = zeros(L,1);
for i = 1:L % for each image
    dis(i) = 1/(m*n) * (norm(Sjv_sum - Pjs(:,i)))^2;
    fprintf("i = %d, dis = %f\n", i, dis(i))
end
[miv, mii] = min(dis)
[mxv, mxi] = max(dis)

% plot images...
T = A(:, mii) + Sjv_sum;
Sj(:,:,1) = reshape(T(1:n*m), [n, m]);
Sj(:,:,2) = reshape(T(n*m+1:2*n*m), [n, m]);
Sj(:,:,3) = reshape(T(n*m*2+1:3*n*m), [n, m]);
Sj_rgb = 1 - Sj;
% truncation processing
Sj_rgb(Sj_rgb < 0) = 0;
Sj_rgb(Sj_rgb > 1) = 1;
hold off
image(Sj_rgb)
hold on
title(strcat("S_j, k = 8, min j idx = ", num2str(mii-1)));
pause

T = A(:, mxi) + Sjv_sum;
Sj(:,:,1) = reshape(T(1:n*m), [n, m]);
Sj(:,:,2) = reshape(T(n*m+1:2*n*m), [n, m]);
Sj(:,:,3) = reshape(T(n*m*2+1:3*n*m), [n, m]);
Sj_rgb = 1 - Sj;
% truncation processing
Sj_rgb(Sj_rgb < 0) = 0;
Sj_rgb(Sj_rgb > 1) = 1;
hold off
image(Sj_rgb)
hold on
title(strcat("S_j, k = 8, max j idx = ", num2str(mxi-1)));

%% (e) read in extra leaves
files = dir('am205_leaves_534x420/extra/*.png');
n = 280;
m = 356;
LE = size(files, 1);
AE = zeros(n*m*3, LE);
ptr = 1;
for file = files'
    Sj = imread(strcat('am205_leaves_534x420/extra/', file.name));
    % Sj is 280x356x3 unit8 (0~255) for H x W x RGB
    % compress into vector 3mn
    Sjv = zeros(n*m*3, 1);
    Sjv(1:n*m) = Sj(:,:,1);         % r
    Sjv(n*m+1:2*n*m) = Sj(:,:,2);   % g
    Sjv(n*m*2+1:3*n*m) = Sj(:,:,3); % b
    
    % convert to cmyk (0 - 255 range orig. new 0 - 1 range)
    Sjv(:) = (255.0 - Sjv(:)) ./ 255.0;
    
    % write into AE
    AE(:,ptr) = Sjv;
    ptr = ptr + 1;
end

% now calculate distances again
PjsE = zeros(3*m*n, LE); % each projection stored in a column vector in Pjs mtrx
for i = 1:LE % note offset by 1 since we are 1-idxd
    % define image vector T
    T = AE(:, i); % random image, :,number
    k = 8;
    P = Sjv_sum; % this is the average vec starting point...
    for j = 1:k
        P = P + transpose(U(:,j)) * (T - Sjv_sum) * U(:,j);
    end
    % copy in
    PjsE(:,i) = P;
end

dis = zeros(LE,1);
for i = 1:LE % for each image
    dis(i) = 1/(m*n) * (norm(Sjv_sum - PjsE(:,i)))^2;
    fprintf("iextra = %d, dis = %f\n", i, dis(i))
end