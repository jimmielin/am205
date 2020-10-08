% read concs
fileID = fopen('measured-concs.txt','r');
formatSpec = '%f %f %f';
sizeArray = [3 Inf];
array = fscanf(fileID, formatSpec, sizeArray);
fclose(fileID);

xvars = array(1,:);
yvars = array(2,:);
cvars = array(3,:); % concs

% 7x7 regular lattice covering -3,...3. set up coords
% we have 200 data points

% construct the matrix to solve in least squares sense for
% Ax = b
% first, concs can be multd by 16*pi to reduce error
cvarssc = cvars .* (16*pi);
% now we just need to fit where
% A_{ij} = exp(-d_{ij}^2/16) where i=1,2,...,200, j=1,2...49
% generate the matrix...

A = zeros(200, 49);
for i = 1:length(cvars)
    for j = 1:49
        % number the lattice like so
        % 1 2 3  4  5  6  7
        % 8 9 10 11 12 13 14
        % ...
        % ...
        %                 49
        % s.t. x coord = j%7 = 1,2,..,0,1,2,...
        xcoord = mod(j, 7);
        if xcoord == 0
            xcoord = 7;
        end
        % y coord = floor(j/7)+1
        ycoord = floor(j/7)+1;
        if xcoord == 7
            ycoord = ycoord - 1;
        end
        % now offset coordinates by 4 since 1 is actually -3
        xcoord = xcoord - 4;
        ycoord = ycoord - 4;
        
%         if i == 1
%             fprintf("i=%d, j=%d\n", xcoord, ycoord)
%         end
        
        A(i,j) = exp(-((xcoord-xvars(i))^2 + (ycoord-yvars(i))^2)/16);
    end
end

% target vector transpose
cvarst = transpose(cvarssc);

% least square sense solution
x = A\cvarst;
x