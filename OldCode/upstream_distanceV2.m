
%distance = zeros(numel(St.x),1); 

% for r = numel(St.ix):-1:1
%     distance(St.ix(r)) = distance(St.ixc(r)) + ...
%         sqrt((St.x(St.ixc(r))-St.x(St.ix(r)))^2 + (St.y(St.ixc(r))-St.y(St.ix(r)))^2);
% end

%area_nal = getnal(S_seg,A);

A_inorder = zeros(numel(area_nal),1);
distance_inorder = zeros(numel(S_seg.distance),1);

for r = numel(S_seg.ix):-1:1
    A_inorder(abs(1+numel(S_seg.ix)-r)) = area_nal(S_seg.ix(r));
    distance_inorder(abs(2+numel(S_seg.ix)-r)) = S_seg.distance(S_seg.ix(r));
end

A_inorder = A_inorder(1:numel(A_inorder)-1);
interpixel_dist = distance_inorder(2:numel(distance_inorder))- ...
    distance_inorder(1:numel(distance_inorder)-1);

%d = sqrt((St.x(St.ix(numel(St.ix)))-St.x(St.ixc(numel(St.ixc))))^2 +...
%    (St.y(St.ix(numel(St.ix)))-St.y(St.ixc(numel(St.ixc))))^2)