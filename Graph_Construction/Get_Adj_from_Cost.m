function [W] = Get_Adj_from_Cost(C)
% Create the adjacency matrix W from a distance (cost) matrix C
% sigma is taken as the mean distance of the points
% in the respective row

n     = length(C(:,1));
sigma = zeros(n,1);

fprintf('---------------------------------------------\n');
fprintf('Adjacency matrix based on Gaussian similarity\n');
fprintf('---------------------------------------------\n');


for i = 1:n
     [distance,index] = sort(C(:,i));
     [row,col,values] = find(distance);
     sigma(i) = mean(values);
end
% % single sigma for all
% sigma_mean = mean(sigma);
% sigma_mean
% W = exp(-2*C.^2 ./ (sigma_mean.^2));

% % \sigma_i is the max between \sigma_i and \sigma_j
Sigma_mat = repmat(sigma,1,n);
Sigma_mat = max(Sigma_mat, Sigma_mat');
W = exp(-1.6*C.^2 ./ (Sigma_mat.^2));

% % \Sigma is \Sigma_i \odot \Sigma_j
% Sigma_mat = repmat(sigma,1,n);
% W = exp(-4*C.^2 ./ (Sigma_mat.*Sigma_mat'));



% Replace the 1s with 0s
W(W==1) = 0;
W = sparse(W);

end