function [ModelKPDistances, OptimParamCombos, Misfit, sum_of_misfit] = RunKPModel_OptimKT_061621(ObsKPDistances, StreamSgmnts, A, m, Kmin, Kmax, Tcap_min, Tcap_max, OutputFilenamePrefix)

%   Inputs:
%   Obs_Distances - N x 3 matrix with KP elevation, StreamNumber (corresponds to column number in StreamSgmnts), and
%   Upstream Length from outlet
%   StreamSgmnts - Cell array of STREAMobjs of streams with knickpoints. Output
%   from SegmentPicker in TAK
%   A - GRIDobj of drainage area for region of interest
%   m - Area exponent from stream-power law, assuming n = 1. Best fitting
%   value shall be used.
%   Kmin - Minimum value of Parameter 1 in model, which is K assuming n=1
%   Kmax - Maximum value of Parameter 1 in model, which is K assuming n=1
%   Tcap_min - Minimum value of the age of capture in model (in years)
%   Tcap_max - Maximum value of the age of capture in model (in years)

%   OutputFilenamePrefix - Text string that will be used to create output
%   '.csv' files. The Tcap value used will be appended to
%   OutputFilenamePrefix in the file name. 
%   Ex of final file name: 'OutputFilenamePrefix_Tcap.csv'


%%% Desired Inputs/Outputs
%
%   KP_Distance - Matrix with the modeled propagation distance for each
%   stream for each parameter combination examined. Number of rows should
%   equal number of individual streams modeled and number of columns should
%   be number of parameter combos (i.e. length(K)*length(Tcap))
%
%   OptimParamCombos - Cell array of best fitting K (1st column), T_cap (2nd column) pairs for the ranges
%   supplied as input. Should also include the sum of misfit for each
%   best fitting combination (column 3) and the index of that misfit value
%   in the sum_of_misfit array (column 4).

%   Will want to have a function that can then generate an output file with
%   x, y, and stream number of modeled knickpoints and also the propagated
%   distance for each knickpoint (although this can be extracted from x, y,
%   data.

K = logspace(log10(Kmin), log10(Kmax), 150); %Generate a log spaced array of K values to test, based on input
Tcap = linspace(Tcap_min, Tcap_max, 19); %Generate a linearly spaced array of capture age (Tcap) values to test, based on input

% Pre-allocate size of output matrix KP_distance
ModelKPDistances = zeros(length(StreamSgmnts), length(K) * length(Tcap));
ModelMaxVels = zeros(length(StreamSgmnts), length(K) * length(Tcap));
ModelMinVels = zeros(length(StreamSgmnts), length(K) * length(Tcap));
ModelMeanVels = zeros(length(StreamSgmnts), length(K) * length(Tcap));

% Loop produces a length(StreamSgmnts) x length(K)*length(Tcap) matrix of
% modeled knickpoint propagation distances upstream. This is the core of
% the KP propagation model.
% for i = 1:length(K)
for i = 1:length(Tcap)
%     for j = 1:length(Tcap)
    for j = 1:length(K)
        
        % Loop over all elements within the StreamSgments cell array, extract
        % vectors of upstream length and area, then input that into the KP
        % propagation model to determine a modeled propagation distance based on
        % the input model parameters.

        for n = 1:length(StreamSgmnts)

            % Call VectorExtractor to output vectors of area and upstream
            % length for the stream being modeled
            [L_vec, A_vec] = VectorExtractor(StreamSgmnts{n}, A);

            % Call KPPropagation_061621 to model knickpoint propagation and output
            % knickpoint propagation distance over the designated time.
            % Value should be a float.
%             ModelOutput = zeros(1,4); %Initialize vector to store output from KP_Propagation_061621
            [ModelKPDistances(n,((i-1)*length(K) + j)), ...
                ModelMaxVels(n,((i-1)*length(K) + j)), ...
                ModelMinVels(n,((i-1)*length(K) + j)), ...
                ModelMeanVels(n,((i-1)*length(K) + j))] = ...
                KPPropagation_061621(L_vec, A_vec, K(j), m, Tcap(i));
            
%             ModelKPDistances(n,((i-1)*length(K) + j)) = ModelOutput(1);
%             ModelMaxVels(n,((i-1)*length(K) + j)) = ModelOutput(2);
%             ModelMinVels(n,((i-1)*length(K) + j)) = ModelOutput(3);
%             ModelMeanVels(n,((i-1)*length(K) + j)) = ModelOutput(4);
%             ModelKPDistances(n,((i-1)*length(K) + j)) = KPPropagation_051721(L_vec, A_vec, K(j), m, Tcap(i));

        end
    end
end

% Now that all of the KP distances have been modeled for all parameter
% combinations (i.e. all combos of K and Tcap), need to calculate misfit
% for each modeled case.

%%% Calculate misfit

% Only retain rows in KP_distance that correspond to a Stream Number
% in Obs_Distances

% Initialize matrix of squared differences for all Streams in Obs_Distances

Misfit = zeros(size(ObsKPDistances,1),length(K) * length(Tcap)); %Initiates a matrix with height equal to the number of observed knickpoints and width equal to the number of modeled K,Tcap pairs to record the differences between modeled and observed KP upstream distance 
OptimParamCombos = cell(length(Tcap),4); %Intiates a cell array to record the optimal combinations of K and Tcap. This enables there to be multiple best fit K values, in the case that optimal is an array with multiple equally good fits (Note: this may likely occur when there's a very poor fit).

% Calculate the sum of squared differences in upstream length for each
% parameter combination. Should be one value per column (i.e. 1 x K*T
% vector)

% Loop to first calculate the squared difference between observed and
% estimated upstream propagation distance for each KP and K, Tcap
% combination
for i = 1:size(ObsKPDistances, 1)
    for j = 1:size(ModelKPDistances, 2)
        
        Misfit(i, j) = (ObsKPDistances(i,3) - ModelKPDistances(ObsKPDistances(i,2),j))^2; % ensures only the StreamNumbers for streams in the Observed_KPs file are compared
        
    end
end

% Sum misfit for all streams with KPs
sum_of_misfit = sum(Misfit,1); %Produces row vector that is sum of columns in Misfit (i.e. the sum of squared differences for each parameter combo)

% Now want to figure out which K value results in the minimum misfit for
% each modeled Tcap value
for i = 1:length(Tcap)
    optimal = find(sum_of_misfit(((i-1)*length(K) + 1):(i*length(K))) == min(sum_of_misfit(((i-1)*length(K) + 1):(i*length(K))))); %Finds the index of the K with the lowest misfit for a given Tcap value in the K array. Subindexed array represents misfits for all Tcaps run for a given K.
    
    OptimParamCombos{i,1} = K(optimal);
    OptimParamCombos{i,2} = Tcap(i);
    OptimParamCombos{i,3} = sum_of_misfit((i-1)*length(K) + optimal); %Finds the sum_of_misfit corresponding to the best fitting K for a given Tcap
    OptimParamCombos{i,4} = (i-1)*length(K) + optimal; % Index of best fitting K for a given Tcap in the sum_of_misfit array. Corresponds to the column with the optimal parameter combo in the KP_distance matrix, as well.
    
end

% Code to output .csv files that can be imported into GIS.

GenerateKPModelOutputFiles_OptimKT_061621(OptimParamCombos, StreamSgmnts, ObsKPDistances, ModelKPDistances, ModelMaxVels, ModelMinVels, ModelMeanVels, OutputFilenamePrefix);

end