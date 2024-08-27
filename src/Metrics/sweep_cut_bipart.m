function [best_cluster,best_conductance] = sweep_cut_bipart(x,W,Lap,normalized)
% [best_cluster,best_conductance] = sweep_cut(x,Lap);
% Retrieves a partition with a sweep_cut on x
%
% INPUT:
% x: vector where the cut will be performed
% Lap: Graph Laplacian
% normalized: 1/0 scale x by the Degree
%
% OUTPUT:
% best_cluster: best partition achieved
% best_conductance: best Cheeger ratio achieved


if normalized
    d = (diag(Lap));
    [xs, idx_s] = sort(x.*d,'descend');
else
    [xs, idx_s] = sort(x,'descend');
end

for i = 1:length(xs) - 1
    indic    = zeros(length(xs),1);
    indic(idx_s(1:i)) = 1;
    sweep(i) = Compute_Cheeger_Cut(Lap,indic);
    % sweep(i) = computeCCutValue(indic,W,1);
end
[best_conductance, min_idx] = min(sweep);

best_cluster = zeros(length(xs),1);
best_cluster(idx_s(1:min_idx)) = 1;

end % function