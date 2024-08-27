function [A] = adjacency_to_incidence(T)
% Convert the adjacency matrix of a graph to its incidence
% Input
% T             : the adjacency matrix (n_nodes x n_nodes)
% Output
% A             : the incidence matrix (n_edges x n_edges)


[adj_rows,adj_cols,~] = find(T);

n = length(adj_rows);

rows = zeros(n,1);
cols = zeros(n,1);
vals = zeros(n,1);

j = 1;
k = 1;

for i = 1:n
    
    if adj_cols(i) < adj_rows(i)
       
       rows(j) =  k;
       rows(j+1) = k;
       
       cols(j) = adj_cols(i);
       cols(j+1) = adj_rows(i);
       
       vals(j) = -1;
       vals(j+1) = 1;
       
       j = j + 2;
       
       k = k + 1;
       
    end
    
end

j = j - 1;

k = k - 1;

A = sparse(rows(1:j),cols(1:j),vals(1:j),k,adj_cols(n));

end