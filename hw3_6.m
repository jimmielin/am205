% is_collide(-1, -1, -2, 1, 0.5) % x0, y0, x1, y1, R

%% (c)

nums = 2500;

num_moon = 0;
num_earth = 0;
num_nothing = 0;

% plot ...
cla reset;
xlim([-1.4 2]);
ylim([-1.2 1.2]);
viscircles([0 0], 0.02, 'Color', [0.0745 0.6235 1.0000], 'LineWidth', 0.5);
hold on
viscircles([1 0], 0.005, 'Color', [0.2 0.2 0.2], 'LineWidth', 0.5);
title("Trajectories (colliding only)")

for i = 1:nums
    % below logic goes in loop
    x0 = [1.0798, 0];
    x1 = [1.0802, -0.0189];
    dt0 = 0.02;

    % randomize
    x0 = x0 + [normrnd(0, 0.002), normrnd(0, 0.002)];
    x1 = x1 + [normrnd(0, 0.002), normrnd(0, 0.002)];

    v0 = (x1-x0)/dt0;

    % now solve ode...
    odeset("RelTol", 1e-10, "AbsTol", 1e-10);

    % init_vec
    % init_dvec
    init_vec = [x0(1) x0(2) v0(1) v0(2)];

    [t, vec] = ode45(@rhs, 0:0.001:10, init_vec);
    
    % check if collision, by iterating through length of vecs
    is_col = 0;
    for j = 1:(size(vec, 1)-1)
       if is_collide(vec(j,1), vec(j,2), vec(j+1,1), vec(j+1,2), 0.02)
           num_earth = num_earth + 1;
           fprintf("=> it will hit the Earth!!! %d\n", i)
           if num_earth == 1
               scatter(vec(:,1), vec(:,2), 0.1)
           end
           
           is_col = 1;
           
           break
       end
       
       if is_collide_moon(vec(j,1), vec(j,2), vec(j+1,1), vec(j+1,2), 0.005)
           num_moon = num_moon + 1;
           fprintf("=> it will hit the Moon! %d\n", i)
           vec(j:j+1,:)
           if num_moon == 1
               scatter(vec(:,1), vec(:,2), 0.1)
           end
           
           is_col = 1;
           
           break
       end
           
       % note, here is a caveat:
       % if it collides with BOTH the earth and the moon, we only compute
       % the first trajectory, of course, because then the asteroid is
       % gone.
    end
    
    if ~is_col
        num_nothing = num_nothing + 1;
    end
    
    % scatter(vec(:,1), vec(:,2), 0.1)
end

%viscircles([0 0], 0.02, 'Color', [0.0745 0.6235 1.0000], 'LineWidth', 1);
%viscircles([1 0], 0.005, 'Color', [0.2 0.2 0.2], 'LineWidth', 1);
xlabel("x");
ylabel("y");

% compute probability...
fprintf("Num Total = %d, Earth = %d, Moon = %d, Nothing = %d\n", ...
        nums, num_earth, num_moon, num_nothing);
fprintf("Probabilities:\n => Earth: %f,\n => Moon: %f\n", ...
        num_earth/nums, num_moon/nums);

% create the target ode...
function dvec_dt = rhs(t, xv)
    % mu param
    mu = 0.01;
    
    % convert
    x = xv(1); y = xv(2); u = xv(3); v = xv(4);
    
    % xv = x, y, u, v
    dx = u;
    dy = v;
    du = v + (x - mu) - (1 - mu)*x/(x^2+y^2)^(3/2) - ...
         mu*(x-1)/((x-1)^2+y^2)^(3/2);
    dv = -u + y - (1 - mu)*y/(x^2+y^2)^(3/2) - mu*y/((x-1)^2+y^2)^(3/2);
    
    dvec_dt = [dx; dy; du; dv];
end

% calculate slope-intercept
function res = is_collide(x0, y0, x1, y1, R)
    % calculate line slope (no check of div/zero)
    if x1-x0 == 0
        k = inf;
        b = x0; % x-intercept
    else
        k = (y1-y0)/(x1-x0);
        b = y0 - k*x0;
    end
    [xout, yout] = linecirc(k, b, 0, 0, R);
    if isnan(xout(1,1))
        res = false;
    else
        if (xout(1) >= min(x0,x1) && xout(1) <= max(x0,x1)) || ...
           (xout(2) >= min(x0,x1) && xout(2) <= max(x0,x1))
            res = true;
        else
            res = false;
        end
    end
end

function res = is_collide_moon(x0, y0, x1, y1, R)
    % calculate line slope (no check of div/zero)
    if x1-x0 == 0
        k = inf;
        b = x0; % x-intercept
    else
        k = (y1-y0)/(x1-x0);
        b = y0 - k*x0;
    end
    [xout, yout] = linecirc(k, b, 1, 0, R);
    if isnan(xout(1,1))
        res = false;
    else
        if (xout(1) >= min(x0,x1) && xout(1) <= max(x0,x1)) || ...
           (xout(2) >= min(x0,x1) && xout(2) <= max(x0,x1))
            res = true;
        else
            res = false;
        end
    end
end