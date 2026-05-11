function [f1] = compute_f1measure(cluster_comp, cluster_gt)
% Compute the F1-score between two clusters
% Input:
%  cluster_comp: computed cluster 
%  cluster_gt: ground truth cluster
% Output:
% f1: F1-score


f1 = 0;

for i = 1:size(cluster_gt, 1)

    precison = size(intersect(cluster_gt(i,:), cluster_comp),2) /size(cluster_comp,2);

    recall = size(intersect(cluster_gt(i,:), cluster_comp),2) /size(cluster_gt(i,:),2);

    if precison == 0 && recall == 0

        f1 = max(f1, 0);

    else

        f1 = max(f1, 2*precison*recall / (precison + recall));

    end

end

