function [KP_distance, max_vel, min_vel, mean_vel] = KPPropagation_061621(L, A, C, p, age)
    %%% Inputs
    %
    %   L : Vector of upstream distance along stream channel
    %   A : Vector of drainage areas along stream channel
    %
    %   C : Parameter 1 in model (needs to be tuned)
    %   p : Parameter 2 in model (needs to be tuned)
    %   age : Age of capture, essentially the time period the model runs

    % Initialize variables
    cum_time = 0;
    n = 2; %Initialize n at 2 to account for using channel outlet in first calculation.
    vel = zeros(1,length(L));


    % For each step between pixels, the duration required to propagate that
    % distance must be calculated. The following needs to be looped the number
    % of pixels in the stream segment, or until the age is reached.
    %while (
    for n = 2:length(L) % Initialize n at 2 to account for using channel outlet in first calculation.
        if cum_time < age
            dt = (L(n)-L(n-1))*((C*A(n)^p)^(-1));
            cum_time = cum_time + dt;
            vel(n) = (L(n)-L(n-1))/dt;
        else
            n = n - 1;
            break
        end
    end

    KP_distance = L(n);
    max_vel = max(vel(2:n)); %Ignore first element since this is zero.
    min_vel = min(vel(2:n));
    mean_vel = mean(vel(2:n));

end
