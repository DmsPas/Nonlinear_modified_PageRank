% Main execution script for the real-world experiment in the paper
%
% "Nonlinear Modified PageRank Problem for Local Graph Partitioning",
% by Costy Kodsi and Dimosthenis Pasadakis.
%
% available online at XXX.
% 
% This program is a free software: you can redistribute it and/or modify
% it under the terms of the Apache-2.0 license.

% Add paths and rng
clear all;close all;
rng(1991);
addpaths_NPR;
%%
cases = {
    % --------------------- %
    %    Orbis Graphs       %
    % --------------------- %
    'Orbis_Km';
    % 'Orbis_Days';
    % 'Orbis_Cost';
    };
fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
fprintf('|                                                    |\n');
fprintf('|               Nonlinear modified Pagerank          |\n');
fprintf('|                                                    |\n');
fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n\n');
%% Parameters
print_level          = 0;    % 0/1/2 for no print/low/high verbosity
p_levels             = [1.95; 1.9; 1.8; 1.7; 1.6; 1.55; 1.5; 1.45]; % levels of p
beta                 = 0.01; % beta = (1 - teleportation constant)/teleportation constant
num_trials           = 3;    % number of different starting nodes
norm_Lap             = 2;    % 0/1/2 for combinatorial/norm. symmetric/norm. random walk Laplacian
write_output_to_file = false;   % write output to file

if write_output_to_file == true
    diary NPR_Orbis_diary.txt;
end

%% Print cases and statistics
nc = length(cases);
maxlen = 0;
for c = 1:nc
    if length(cases{c}) > maxlen
        maxlen = length(cases{c});
    end
end

num_cases = 0;
for c = 1:nc
    fprintf('.');
    temp = load(cases{c});
    %     if ~isfield(temp,'coords')
    %         temp.coords = [];
    %     end
    Graph_n_labels(c) = temp;
    num_cases = num_cases + 1;
end

node_cluster_ratio = zeros(nc,1);

fprintf('\n\n Report Cases - %2d %12s %9s %10s %20s\n',num_cases,'Nodes','Edges','Clusters','Clusters/Nodes');
fprintf('------------------------------------------------------\n');
for c = 1:nc
    % get number of vertices
    params.numberOfVertices = size(Graph_n_labels(c).W,1);
    % get number of edges
    params.numberOfEdges    = nnz(Graph_n_labels(c).W)/2;
    % get number of clusters
    params.numberOfClusters = size(unique(Graph_n_labels(c).label),1);
    % calculate the cl/node ratio
    node_cluster_ratio(c)   = params.numberOfClusters/params.numberOfVertices;
    % spacing
    spacers = repmat('.', 1, maxlen+3-length(cases{c}));
    % print headers
    fprintf('%s %s %18d %10d %10d %20.3f\n',cases{c},spacers,params.numberOfVertices, ...
        params.numberOfEdges,params.numberOfClusters,node_cluster_ratio(c));
end
fprintf('------------------------------------------------------\n\n');

% Initialize variables for results
RCCut_all           = zeros(num_trials,nc);
FSCORE_all          = zeros(num_trials,nc);
best_p_all          = zeros(num_trials,nc);
time_all            = zeros(num_trials,nc);
start_nodes_all     = zeros(num_trials,nc);
best_cluster_all    = zeros(params.numberOfVertices,num_trials);
best_card_all       = zeros(params.numberOfVertices,num_trials);
dist_from_seed_all  = zeros(params.numberOfVertices,num_trials);
labels_dist_all     = zeros(params.numberOfVertices,num_trials);
dist_from_seed_all1 = zeros(params.numberOfVertices,num_trials);
%% Loop over cases
for orig_ID = 1:nc
    c = orig_ID;
    spacers          = repmat('.', 1, maxlen+3-length(cases{c}));
    W                = Graph_n_labels(c).W;
    C                = Graph_n_labels(c).C;
    labels_full      = Graph_n_labels(c).label;
    coords           = Graph_n_labels(c).coords;
    % ensure label numbering starts from 1
    if min(labels_full) == 0
        labels_full = labels_full + 1;
        fprintf('Labels start from 0. Adding 1.\n');
    elseif min(labels_full) == 1
        fprintf('Labels start from 1.\n');
    else
        fprintf('Labelling error, min(labels) = %f\n',min(labels_full));
    end
    % Read graph information
    n_clusters      = size(unique(Graph_n_labels(c).label),1);
    n_nodes         = size(W,1);
    n_edges         = nnz(W)/2;
    % Starting nodes are fixed in this experiment
    % They correspond to the cities of [Constantinople, Londinium, Rome]
    start_nodes    = [126,232,324]'; 

    start_nodes_all(:,orig_ID) = start_nodes;

    % Create graph Laplacian and Incidence matrices
    [L,Deg,vw]  = CreateLapl(W,norm_Lap);
    % Build incidence matrix
    B           = adjacency_to_incidence(W);
    % Moore-penrose pseudoinverse, B^+
    B_pinv      = pinv(full(B));

    % M = \beta D + L, used in the objective
    M           = beta.*Deg + L;
    % N = (\beta D + L) * B^+, used in the objective and Jacobian
    N           = M * B_pinv;
    % Used to estimate an initial solution
    Pinv_N_B    = pinv(N * B);

    for tr = 1:length(start_nodes)
        % Initialize RCCut and plevels
        RCCut_p_levs  = zeros(length(p_levels));
        Fscore_p_levs = zeros(length(p_levels));
        % Current starting node
        start_node   = start_nodes(tr);

        % Initialize starting vector
        s             = zeros(n_nodes,1);
        s(start_node) = 1.0;
        s             = beta .* s;

        % Initialize solution vector
        x = Pinv_N_B * s;

        % Shortest paths on graph
        Gr              = graph(W);
        dist_from_seed  = distances(Gr,start_node);
        dist_from_seed_all(:,tr) = dist_from_seed;
        % Furthest node from seed 
        [~,vrt_idx]     = max(dist_from_seed_all(:,tr));
        % in case the node numbering of the graph is not sequential        
        [vrt_max,~]     = ind2sub([n_nodes n_nodes],vrt_idx);
        
        % Shortest paths on graph
        Gr1              = graph(C);
        dist_from_seed1  = distances(Gr1,start_node);
        dist_from_seed_all1(:,tr) = dist_from_seed1;
        % Sort the nodes based on distance on graph
        [sorted_dists,sorted_nodes] = sort(dist_from_seed_all1(:,tr));

        % Set options for the Lev-Mar optimization
        opt.tol_1       = 1.0e-5;
        opt.tol_2       = 1.0e-5;
        opt.max_iter    = 500;
        opt.vrt_max     = vrt_max;
        opt.vrt_max_val = 1e-12;
        opt.lambda_c    = 1.0e-3;

        % FUNCTION CALL FOR P-LEVELS
        tStart_NonlinPR = tic; % Initialize timer
        [best_p, best_cond, best_Fscore,best_cluster,best_card,labels_dist_curr] = ...
            NPR_levels_Orbis(start_node,x,W,L,N,B,s,sorted_nodes,p_levels,opt, print_level);
        best_p_all(tr,orig_ID)    = best_p;
        RCCut_all(tr,orig_ID)     = best_cond;
        FSCORE_all(tr,orig_ID)    = best_Fscore;
        time_all(tr,orig_ID)      = toc(tStart_NonlinPR);
        best_cluster_all(:,tr)    = best_cluster;
        best_card_all(tr,orig_ID) = best_card;
        labels_dist_all(:,tr)     = labels_dist_curr;
    end % node trials

end % num cases

for orig_ID = 1:nc
    W = Graph_n_labels(orig_ID).W;
    n_nodes     = size(W,1);
    n_edges      = nnz(W)/2;
    fprintf('=============\n');
    fprintf('%8s, nodes:%2d, edges:%2d, time(sec):%f\n',cases{orig_ID},...
        n_nodes, n_edges, sum(time_all(:,orig_ID)));
    fprintf('=============\n');
    fprintf('%s %10s %7s %10s %12s %12s %12s \n', 'trial', 's_node', 'best_p',...
        'RCCut', 'Cardinality','F-score','time(sec)')
    for tr = 1:num_trials
        fprintf('%2d %10d %5f %15f %10f %10f %10f  \n', tr, start_nodes_all(tr,orig_ID),  ...
            best_p_all(tr,orig_ID), RCCut_all(tr,orig_ID), best_card_all(tr,orig_ID), ...
            FSCORE_all(tr,orig_ID),time_all(tr,orig_ID));
    end
end

% Print mean results
fprintf('=============\n');
fprintf('%10s\n','Mean results');
fprintf('%s %10s %10s %12s \n', 'Case', 'RCCut', 'F-score', 'time(sec)')
fprintf('=============\n');
for orig_ID = 1:nc
    fprintf('%2s %10f %10f %10f \n', cases{orig_ID}, mean(RCCut_all(:,orig_ID)), ...
        mean(FSCORE_all(:,orig_ID)), mean(time_all(:,orig_ID)));
end
fprintf('=============\n');

% Visualize the results on map
Visualize_Orbis_Results(W,labels_full,labels_dist_all,coords,best_cluster_all);