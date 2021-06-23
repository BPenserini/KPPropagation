function KP_distance = KPPropagation_051721(L, A, C, p, age)
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


    % For each step between pixels, the duration required to propagate that
    % distance must be calculated. The following needs to be looped the number
    % of pixels in the stream segment, or until the age is reached.
    %while (
    for n = 2:length(L) % Initialize n at 2 to account for using channel outlet in first calculation.
        if cum_time < age
            dt = (L(n)-L(n-1))*((C*A(n)^p)^(-1));
            cum_time = cum_time + dt;
        else
            n = n - 1;
            break
        end
    end

    KP_distance = L(n);

end
