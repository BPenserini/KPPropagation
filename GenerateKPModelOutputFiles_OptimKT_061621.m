function GenerateKPModelOutputFiles_OptimKT_061621(OptimParamCombos, StreamSgmnts, ObsKPDistances, ModelKPDistances, ModelMaxVels, ModelMinVels, ModelMeanVels, OutputFilenamePrefix)

% Inputs
%   OptimParams - N x 2 cell array consisting of optimal K, Tcap
%   combinations, as prescribed by RunKPModel_OptimKT_060921. First column
%   is best fitting K value for Tcap in second column
%
%   StreamSgmnts - Cell array of STREAMobjs of streams with knickpoints. 
%   Output from SegmentPicker in TAK
%
%   ObsKPDistances - N x 3 matrix with observed KP elevation, corresponding StreamNumber (column in StreamSgmnts), and
%   corresponding upstream length from outlet.
%
%   ModelDistances - M x N matrices with modeled knickpoint propagation
%   distances and max, min, and mean horizontal knickpoint velocities. Each of the M rows represents an instance of StreamSgmnts
%   run through the KP propagation model, while each of the N rows
%   represents a different parameter combination used to model KP
%   propagation. See Code in RunKPModel_OptimKT_061621 for details as to
%   which columns refer to which parameter combos.
%
%   OutputFilenamePrefix - Text string that will be used to create output
%   '.csv' files. The Tcap value used will be appended to
%   OutputFilenamePrefix in the file name. 
%   Ex of final file name: 'OutputFilenamePrefix_Tcap.csv'

% Outputs
%
%
%
%
%
%

for i = 1:size(OptimParamCombos,1)
    
    %%%Create matrix for modeled knickpoint locations (X, Y, upstream propagation distance, max velocity, min velocity, mean velocity, and StreamNumber)
    ModeledKPLocations = zeros(size(ObsKPDistances,1),7);
    
    % Loop over all STREAMobjs in StreamSgmnts that have an observed knickpoint
    % for comparison
    for j = 1:size(ObsKPDistances,1)
        
        % Identify the modeled upstream propagation distance for the particular
        % stream in this loop iteration
        L_KP = ModelKPDistances(ObsKPDistances(j,2),OptimParamCombos{i,4});
        max_vel = ModelMaxVels(ObsKPDistances(j,2),OptimParamCombos{i,4});
        min_vel = ModelMinVels(ObsKPDistances(j,2),OptimParamCombos{i,4});
        mean_vel = ModelMeanVels(ObsKPDistances(j,2),OptimParamCombos{i,4});
        
        % Use find() to identify the index of the knickpoint location in the
        % distance NAL within the STREAMobj of this loop iteration
        KP_index = find(L_KP == StreamSgmnts{ObsKPDistances(j,2)}.distance);
        
        % Use that index to extract X and Y location of the modeled knickpoint
        X_KP = StreamSgmnts{ObsKPDistances(j,2)}.x(KP_index);
        Y_KP = StreamSgmnts{ObsKPDistances(j,2)}.y(KP_index);
        
        % Save X, Y, upstream propagation distance and StreamNumber as a row in a N x 3 matrix, which can
        % then be used to bring data into GIS
        ModeledKPLocations(j,:) = [X_KP Y_KP L_KP max_vel min_vel mean_vel ObsKPDistances(j,2)];     
    
    end  
    
    T = array2table(ModeledKPLocations);
    T.Properties.VariableNames(1:7) = {'X', 'Y', 'L_m_Upstream','Max_Vel','Min_Vel', 'Mean_Vel', 'StreamNumber'};
    writetable(T, [OutputFilenamePrefix '_' num2str(OptimParamCombos{i,2}) '.csv']);
        
end

end