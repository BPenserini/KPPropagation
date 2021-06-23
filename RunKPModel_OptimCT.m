function [KP_distance, Modeled_KP_Locations, sum_of_misfit, OptimParamCombos] = RunKPModel_OptimCT(Obs_Distances, StreamSgmnts, A, Cmin, Cmax, Tmin, Tmax, p) %, output_filename)
%function [KP_distance, sum_of_misfit, C_optimal, p_optimal] = RunKPModel_051721(Obs_Distances, StreamSgmnts, A, Cmin, Cmax, pmin, pmax, age)

%   Inputs:
%   Obs_Distances - N x 3 matrix with KP elevation, StreamNumber, and
%   Upstream Length from outlet
%   StreamSgmnts - Cell array of STREAMobjs of streams with knickpoints. Output
%   from SegmentPicker in TAK
%   A - GRIDobj of drainage area for region of interest
%   Cmin - Minimum value of Parameter 1 in model
%   Cmax - Maximum value of Parameter 1 in model
%   Tmin - Minimum value of the age of capture in model
%   Tmax - Maximum value of the age of capture in model
%   p - Parameter 2 in the model. For n~1, this value is m and is set as
%   a known value which is one of the inputs
%   output_filename - string that is the name of the output file with x, y,
%   and stream number of modeled knickpoints. Need file extension (i.e.
%   '.csv')

%%% Desired Inputs/Outputs
%
%   L : Vector of upstream distance along stream channel
%   A : Vector of drainage areas along stream channel
%   X : Vector of X locations along stream channel
%   Y : Vector of Y locations along stream channel
%   Z : Vector of channel elevations

C = logspace(log10(Cmin), log10(Cmax), 50); %Generate a linearly spaced array of C values to test, based on input
T = linspace(Tmin, Tmax, 10); %Generate a linearly spaced array of T values to test, based on input

% Pre-allocate size of output matrix KP_distance
KP_distance = zeros(length(StreamSgmnts), length(C) * length(T));

for i = 1:length(C)
    for j = 1:length(T)
        
        % Loop over all elements within the StreamSgments cell array, extract
        % vectors of upstream length and area, then input that into the KP
        % propagation model to determine a modeled propagation distance based on
        % the input model parameters.

        for n = 1:length(StreamSgmnts)

            % Call VectorExtractor to output vectors of area and upstream length
            [L_vec, A_vec] = VectorExtractor(StreamSgmnts{n}, A);

            % Call KPPropagation_051721 to model knickpoint propagation and output
            % knickpoint propagation distance over the designated time
            KP_distance(n,((i-1)*10 + j)) = KPPropagation_051721(L_vec, A_vec, C(i), p, T(j));

        end
    end
end


%%% Calculate misfit

% Only retain rows in KP_distance that correspond to a Stream Number
% in Obs_Distances

% Initialize matrix of squared differences for all Streams in Obs_Distances

Misfit = zeros(size(Obs_Distances,1),length(C) * length(T));
OptimParamCombos = zeros(length(C),2);

% Calculate the sum of squared differences in upstream length for each
% parameter combination. Should be one value per column (i.e. 1 x C*p
% vector)

for i = 1:size(Obs_Distances, 1)
    for j = 1:size(KP_distance, 2)
        
        Misfit(i, j) = (Obs_Distances(i,3) - KP_distance(Obs_Distances(i,2),j))^2; % ensures only the StreamNumbers for streams in the Observed_KPs file are compared
        
    end
end
    
sum_of_misfit = sum(Misfit,1); %Produces row vector that is sum of columns in Misfit (i.e. the sum of squared differences for each parameter combo)

for i = 1:length(C)
    optimal = find(sum_of_misfit(((i-1)*10 + 1):(i*10)) == min(sum_of_misfit(((i-1)*10 + 1):(i*10)))); %Finds the T with the lowest misfit for a given C value. Subindexed array represents misfits for all Ts run for a given C.
    
    OptimParamCombos(i,1) = C(i);
    OptimParamCombos(i,2) = T(optimal);
end


% optimal = find(sum_of_misfit == min(sum_of_misfit));
% 
% C_optimal = C(floor(optimal/10))
% 
% if mod(optimal,10) == 0
%     T_optimal = T(10)
% else
%     T_optimal = T(mod(optimal,10))
% end
    
%%%Create matrix for modeled knickpoint locations (X, Y, and StreamNumber)
Modeled_KP_Locations = zeros(size(Obs_Distances,1),3);

% Loop over all STREAMobjs in StreamSgmnts that have an observed knickpoint
% for comparison
for i = 1:size(Obs_Distances,1)
    
    % Identify the modeled upstream propagation distance for the particular
    % stream in this loop iteration
    modeled_distance = KP_distance(Obs_Distances(i,2), optimal);
    
    % Use find() to identify the index of the knickpoint location in the
    % distance NAL within the STREAMobj of this loop iteration
    KP_index = find(modeled_distance == StreamSgmnts{Obs_Distances(i,2)}.distance);
    
    % Use that index to extract X and Y location of the modeled knickpoint
    X_KP = StreamSgmnts{Obs_Distances(i,2)}.x(KP_index);
    Y_KP = StreamSgmnts{Obs_Distances(i,2)}.y(KP_index);
    
    % Save X, Y, and StreamNumber as a row in a N x 3 matrix, which can
    % then be used to bring data into GIS
    Modeled_KP_Locations(i,:) = [X_KP Y_KP Obs_Distances(i,2)]; 
    
end

%writematrix(Modeled_KP_Locations,output_filename);

end