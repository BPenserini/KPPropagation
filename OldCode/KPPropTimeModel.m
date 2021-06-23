function [dt] = KPPropTimeModel(A,Sc,C,p)
%For Debugging
%function [dt] = KPPropTimeModel(DEM,FD,A,S,DEMc,KP,C,p)

% Turn this into a function where the inputs are KP location, initiation
% location, C value/vector, and p value/vector

% Inputs:
%
% DEM - GRIDobj of digital elevation data. Needs to be resampled to whole
%   number using 'resample_grid' tag in MakeStreams (TAK) prior to
%   supplying as input
%
% FD - FLOWobj of flow directions, output by MakeStreams (TAK).
%
% S - STREAMobj of stream network, output by MakeStreams (TAK). This
%   function will treat the outlet of S as the desired downstream extent of
%   the channel segment that is extracted.
%
% A - GRIDobj of flow accumulation (i.e. drainage area) for the basin of
%   interest
%
% Sc - STREAMobj that is clipped to be the channel segment between
%   initiation point and knickpoint.
%
% DEMc - Conditioned DEM, output of ConditionDEM (TAK). Do not use as DEM
%   input.
%
% KP - 1x3 (for now) matrix of channel heads with x coord, y coord, and pick number
%   as the columns.
%
% C - vector containing different values of the C parameter for
%   area-dependent knickpoint migration.
%
% p -  vector containing different values of the p parameter for
%   area-dependent knickpoint migration.


% Output:
%
% dt - Matrix of simulated knickpoint propagation times. Units are
%   dependent on calibration and selection of C and p. Row index
%   corresponds to C index and column index corresponds to p index.



% [Sc] = SegmentPicker(DEM,FD,A,S,1,...
%     'conditioned_DEM',DEMc,...
%     'direction','down',...
%     'method','prev_picks',...
%     'ref_concavity',0.45,...
%     'picks',KP);

A_nal = getnal(Sc, A); % Extracts a node-attribute list for drainage area along stream segment.

A_inorder = zeros(numel(A_nal),1); % Initializes vector for drainage area along stream segment.
dist_inorder = zeros(numel(Sc.distance),1); % Initializes vector for interpixel distances.
dt = zeros(numel(C),numel(p)); % Initializes matrix to store simulated propogation times.


for r = numel(Sc.ix):-1:1 % Reorders the NALs so that they correctly represent pixels flowing into each other
    
    A_inorder(abs(1+numel(Sc.ix)-r)) = A_nal(Sc.ix(r));
    dist_inorder(abs(2+numel(Sc.ix)-r)) = Sc.distance(Sc.ix(r));
end

A_inorder = A_inorder(1:numel(A_inorder)-1); 
dist_inorder = dist_inorder(2:numel(dist_inorder))- ...
    dist_inorder(1:numel(dist_inorder)-1);




for i = 1:size(dt,1)    % Loops through rows of the dt matrix, each row represents a different C value
    
    for j = 1:size(dt,2) % Loops through columns of the dt matrix, each column represents a different p value
       
        for r = 1:numel(dist_inorder) % Loops through model, calculating travel time between each pixel.
            
            dt(i,j) = dt(i,j) +  ...
                dist_inorder(r)*(C(i)*(A_inorder(r)^p(j)))^(-1); % Area dependent propogation equation
       
        end
    end
end
end















% function [dt] = KPPropTimeModel(C, p, distance, area)
% % Turn this into a function where the inputs are KP location, initiation
% % location, C value/vector, and p value/vector
% 
% dt = zeros(numel(C),numel(p)); % initializes matrix to store simulated propogation times.
% 
% for i = 1:size(dt,1)    % Loops through rows of the dt matrix, each row represents a different C value
%     
%     for j = 1:size(dt,2) % Loops through columns of the dt matrix, each column represents a different p value
%        
%         for r = 1:numel(distance) % Loops through model, calculating travel time between each pixel.
%             
%             dt(i,j) = dt(i,j) +  ...
%                 distance(r)*(C(i)*(area(r)^p(j)))^(-1); % Area dependent propogation equation
%        
%         end
%     end
% end
% end