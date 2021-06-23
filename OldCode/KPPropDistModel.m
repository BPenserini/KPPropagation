%function [coords, prop_length, A_nal, A_inorder, dist_inorder, Sc, coord_index, dt] = KPPropDistModel(DEM,FD,A,S,DEMc,CH,C,p,age)
%For debugging

function [coords] = KPPropDistModel(DEM,FD,A,S,DEMc,CH,C,p,age)
% Want to output predicted coordinates of knickpoint. 

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
% S - STREAMobj that is clipped to be the channel segment between
%   initiation point and knickpoint.
%
% DEMc - Conditioned DEM, output of ConditionDEM (TAK). Do not use as DEM
%   input.
%
% CH - 1x3 (for now) matrix of channel heads with x coord, y coord, and pick number
%   as the columns. Each row is output of 'streampoi' with 'channel_heads'
%   selected, plus an integer value for the third column. (See
%   SegmentPicker help)
%
% C - vector containing different values of the C parameter for
%   area-dependent knickpoint migration.
%
% p -  vector containing different values of the p parameter for
%   area-dependent knickpoint migration.
%
% age - Matrix of number in years representing the presumed age of the knickpoint(s).
%   Serves as a flag to exit out of the loop that propagates the
%   knickpoint. Dimensions are numel(C) x numel(p)


% Output:
%
% coords - Simulated knickpoint propagation location(s), based on input  


% TRY TO UPDATE SegmentPicker SO IT INCLUDES CHANNEL SEGMENT THAT GOES
% UPSTREAM OF KNICKPOINT. CURRENTLY THE INPUT COORDINATES NEED TO BE OF THE
% CHANNEL HEAD WITH THE KNICKPOINT.

coords = zeros(size(CH,1),3);

for k = 1:size(CH,1)

    [Sc] = SegmentPicker(DEM,FD,A,S,1,... % 1 is a placeholder for the basin number, will need to change if want to take multiple inputs
        'conditioned_DEM',DEMc,...
        'direction','down',...
        'method','prev_picks',...
        'ref_concavity',0.45,...
        'picks',CH(k,:));

    A_nal = getnal(Sc, A); % Extracts a node-attribute list for drainage area along stream segment.

    A_inorder = zeros(numel(A_nal),1); % Initializes vector for drainage area along stream segment.
    dist_inorder = zeros(numel(Sc.distance),1); % Initializes vector for interpixel distances.
    dt = zeros(numel(C),numel(p)); % Initializes matrix to store simulated propagation times.
    prop_length = 0; % Initializes value to store total length KP propagated.

    for r = numel(Sc.ix):-1:1 % Reorders the NALs so that they correctly represent pixels flowing into each other

        A_inorder(abs(1+numel(Sc.ix)-r)) = A_nal(Sc.ix(r));
        dist_inorder(abs(2+numel(Sc.ix)-r)) = Sc.distance(Sc.ix(r));
    end

    A_inorder = A_inorder(1:numel(A_inorder)-1); 
    dist_inorder = dist_inorder(2:numel(dist_inorder))- ...
        dist_inorder(1:numel(dist_inorder)-1);

     % Loops through columns of the dt matrix, each column represents a different p value

    for i = 1:size(dt,1)    % Loops through rows of the dt matrix, each row represents a different C value

        for j = 1:size(dt,2) % Loops through columns of the dt matrix, each column represents a different p value

            for r = 1:numel(dist_inorder) % Loops through model, calculating travel time between each pixel.

                time_inc = dist_inorder(r)*(C(i)*(A_inorder(r)^p(j)))^(-1);

                if (dt(i,j) + time_inc) >= age % Break out of loop if dt has exceeded the specified age.
                    break
                else
                    dt(i,j) = dt(i,j) + time_inc;
                    prop_length = prop_length + dist_inorder(r);
                    coord_index = find(A_nal == A_inorder(r));
                end

                %dt(i,j) = dt(i,j) +  ...
                %    dist_inorder(r)*(C(i)*(A_inorder(r)^p(j)))^(-1); % Area dependent propogation equation
            end
            coords(k,:) = [Sc.x(coord_index) Sc.y(coord_index) coord_index];
        end
    end
end
end



% Similar to above, but assumes singular C,p values
%
% for r = 1:numel(dist_inorder) % Loops through model, calculating travel time between each pixel.
% 
%     time_inc = dist_inorder(r)*(C*(A_inorder(r)^p))^(-1);
% 
%     if (dt + time_inc) >= age % Break out of loop if dt has exceeded the specified age.
%         break
%     else
%         dt = dt + time_inc;
%         prop_length = prop_length + dist_inorder(r);
%         coord_index = find(A_nal == A_inorder(r));
%     end
% 
%     %dt(i,j) = dt(i,j) +  ...
%     %    dist_inorder(r)*(C(i)*(A_inorder(r)^p(j)))^(-1); % Area dependent propogation equation
% end
% 
% coords = [Sc.x(coord_index) Sc.y(coord_index) coord_index];
% end









% function [dt] = RunKPModel(C, p, distance, area)
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