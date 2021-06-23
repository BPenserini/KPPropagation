
%distance = zeros(numel(St.x),1); 

% for r = numel(St.ix):-1:1
%     distance(St.ix(r)) = distance(St.ixc(r)) + ...
%         sqrt((St.x(St.ixc(r))-St.x(St.ix(r)))^2 + (St.y(St.ixc(r))-St.y(St.ix(r)))^2);
% end

A_inorder = zeros(numel(area_nal)-1,1);
distance_inorder = zeros(numel(distance_nal)-1,1);

for r = numel(S_seg.ixc):-1:1
    A_inorder(abs(1+numel(S_seg.ixc)-r)) = area_nal(S_seg.ixc(r));
    distance_inorder(abs(1+numel(S_seg.ixc)-r)) = distance_nal(S_seg.ixc(r));
end

%d = sqrt((St.x(St.ix(numel(St.ix)))-St.x(St.ixc(numel(St.ixc))))^2 +...
%    (St.y(St.ix(numel(St.ix)))-St.y(St.ixc(numel(St.ixc))))^2)