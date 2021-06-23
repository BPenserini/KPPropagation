function [dt] = RunKPModel(S, A, C, p)
% Turn this into a function where the inputs are KP location, initiation
% location, C value/vector, and p value/vector

% Inputs:
%
% S - STREAMobj that is clipped to be the channel segment between
%   initiation point and knickpoint.
%
% A - GRIDobj of flow accumulation (i.e. drainage area) for the basin of
%   interest
%
% C - vector containing different values of the C parameter for
%   area-dependent knickpoint migration.
%
% p -  vector containing different values of the p parameter for
%   area-dependent knickpoint migration.


A_nal = getnal(S, A); % Extracts a node-attribute list for drainage area along stream segment.

A_inorder = zeros(numel(A_nal),1); % Initializes vector for drainage area along stream segment.
dist_inorder = zeros(numel(S.distance),1); % Initializes vector for interpixel distances.
dt = zeros(numel(C),numel(p)); % Initializes matrix to store simulated propogation times.


for r = numel(S.ix):-1:1 % Reorders the NALs so that they correctly represent pixels flowing into each other
    
    A_inorder(abs(1+numel(S.ix)-r)) = A_nal(S.ix(r));
    dist_inorder(abs(2+numel(S.ix)-r)) = S.distance(S.ix(r));
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