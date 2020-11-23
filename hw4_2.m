% HW4-2
% fill arrays
map = load("pierce.txt");

js = 100;
ks = 200;
h = 36.6;
c = 3.43*1e4;
p0 = 10;
w = 100*pi;
dt = h/(2*c);

% time
t0 = 0.000;
te = 2.500; % 1.005
ns = ceil((te-t0)/dt); % dt ~ 0.0005

% P stored in cells for horiz, vert, time
% indexes, not distances (k,j,t)
% the indices here are NOT the Pn,j,k = (kh, (99-j)h, ndt) coords!
P = zeros(ks, js, ns);

% fill initial values, note that ns are 1-indexed
% and everything else is too, so shift one
P(:,:,1) = 0.0;
P((15+1):(18+1), (57+1):(60+1), 2) = p0 * sin(w * dt); % memorder KJN

% timestep loop
for np1 = 3:ns
   n = np1 - 1; % n = 0,1,2... only run from 2,...
   % but np1 is used for idx, as idx runs from 1:... (loop from 3,...)
   % so n=0,1,2,...,ns-1 => np1=1,2,...ns
   
   % advance a time step for every point in the grid
   for jp1 = 1:js % true index plus one one
        j = jp1-1;
        for kp1 = 1:ks
            k = kp1-1;
            
            val = 0.0; % initial val
            % enforce throughout within S ...
            if k >= 15 && k <= 18 && j >= 57 && j <= 60
                val = p0 * sin(w * n * dt);
            else
                %fprintf("\n\n np1,j,k=%d,%d,%d\n", np1, j, k);
                
                val = c^2*dt^2/(h^2) * ...
                        (P_acc(map, P, k, j+1, np1-1, k, j) + ...
                         P_acc(map, P, k, j-1, np1-1, k, j) + ...
                         P_acc(map, P, k-1, j, np1-1, k, j) + ...
                         P_acc(map, P, k+1, j, np1-1, k, j) - ...
                         4*P_acc(map, P, k, j, np1-1, k, j)) + ...
                      2*P_acc(map, P, k, j, np1-1, k, j) - ...
                        P_acc(map, P, k, j, np1-2, k, j);
            end
            
            %P(kp1,99-j+1,np1) = val; % memorder KJN
            P(kp1,jp1,np1) = val; % memorder KJN
        end
   end
   
   %tmp_plot = transpose(P(:,:,np1));
   %custom_plot1('test.png',tmp_plot,map,-10.0,10.0,3,dt*n);
end


%% analyse the data...
% (c)
cData = reshape(P(73+1,35+1,:),1,[]);
gData = reshape(P(109+1,61+1,:),1,[]);
mData = reshape(P(188+1,91+1,:),1,[]);
ts = t0:dt:te;
plot(ts, cData);
hold on;
plot(ts, gData);
plot(ts, mData);
xlim([0 1]);
title("p(t) for each person");
legend(["C" "G" "M"]);
xlabel("t [s]");
ylabel("p(t)");

for i = 1:size(ts,2)
   if abs(cData(i)) > 0.001
       fprintf("for C, first |p(t)|>0.001 is t = %f\n", ts(i));
       break
   end
end

for i = 1:size(ts,2)
   if abs(gData(i)) > 0.001
       fprintf("for G, first |p(t)|>0.001 is t = %f\n", ts(i));
       break
   end
end

for i = 1:size(ts,2)
   if abs(mData(i)) > 0.001
       fprintf("for M, first |p(t)|>0.001 is t = %f\n", ts(i));
       break
   end
end

% a handy access function to handle bdy conditions,
% returning the correct value
function res = P_acc(map, P, k, j, np1, curr_k, curr_j)
    ka = k;
    %ja = (99-j);
    ja = j;
    
    % handle oob - simply return zero
    if j <= 0 || k <= 0 || j >= 100 || k >= 200
        res = 0.0;
        return
    end
    
    %fprintf("np1,j,k,cj,ck=%d,%d,%d,%d,%d\n", np1, j, k, curr_j, curr_k);
    
    if abs(curr_k - k) == 1 && curr_j == j && map(j+1, k+1) == 1
        ka = curr_k;
        %ja = (99-curr_j);
        ja = curr_j;
    end
    
    if abs(curr_j - j) == 1 && curr_k == k && map(j+1, k+1) == 1
        ka = curr_k;
        ja = curr_j;
        %ja = (99-curr_j);
    end % treat as equivalent
    
    res = P(ka+1, ja+1, np1);
end