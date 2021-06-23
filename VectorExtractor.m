function [L, A] = VectorExtractor(S, A_grid)%, DEM)

% Inputs:
%
% S - STREAMobj that is clipped to be the channel along which there is a knickpoint.
%
% A_grid - GRIDobj of flow accumulation (i.e. drainage area) for the basin
% containing S

A_nal = getnal(S, A_grid); % Extracts a node-attribute list for drainage area along stream segment.
L_nal = S.distance;
%Z_nal = getnal(S,DEM);
A = zeros(numel(A_nal),1); % Initializes vector for drainage area along stream segment.
L = zeros(numel(S.distance),1); % Initializes vector for interpixel distances.
%Z = zeros(numel(Z_nal),1);

for r = numel(S.ix):-1:1 % Reorders the NALs so that they correctly represent pixels flowing into each other
    
    A(abs(2+numel(S.ix)-r)) = A_nal(S.ix(r));
    L(abs(2+numel(S.ix)-r)) = L_nal(S.ix(r));
    %Z(abs(2+numel(S.ix)-r)) = Z_nal(S.ix(r));
end

%Need to account for ix only being of giver nodes and add the values for
%the outlet
A(1) = max(A_nal);
L(1) = min(L_nal);
%Z(1) = min(Z_nal);

%Convert A to m^2
A = A*30*30;

end