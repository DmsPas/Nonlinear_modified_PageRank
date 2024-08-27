% Main execution script for the paper
%
% "Nonlinear Modified PageRank Problem for Local Graph Partitioning",
% by Costy Kodsi and Dimosthenis Pasadakis.
%
% available online at XXX.
% 
% This program is a free software: you can redistribute it and/or modify
% it under the terms of the Apache-2.0 license.

% Add paths and rng
clear all; close all; rng(1991);
addpaths_NPR;
%%
cases = {
% % Synthetic LFR graphs
'LFR_10.mat';
% 'LFR_12.mat';
% 'LFR_14.mat';
% % Synthetic Gauss graphs
% 'Gauss2_55_10NN';
% 'Gauss5_55_10NN';
% 'Gauss8_55_10NN';
};
fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
fprintf('|                                                    |\n');
fprintf('|               Nonlinear modified Pagerank          |\n');
fprintf('|                                                    |\n');
fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n\n');
%% Parameters
print_level          = 2;     % 0/1/2 for no print/low/high verbosity
p_levels             = [1.95; 1.9; 1.8; 1.7; 1.6; 1.5; 1.45]; % levels of p
% beta = (1 - teleportation constant)/teleportation constant
beta               = 0.01;      % for the LFR cases
% beta                 = 0.0001;% for the Gauss cases
num_trials           = 10;      % number of different starting nodes
norm_Lap             = 2;       % 0/1/2 for combinatorial/norm. symmetric/norm. random walk Laplacian
write_output_to_file = false;   % write output to file

if write_output_to_file == true
    diary NPR_diary.txt
end

fprintf('Normalized Laplacian is: %d\n',norm_Lap);
fprintf('beta parameter       is: %f\n',beta);

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
RCCut_all        = zeros(num_trials,nc);
FSCORE_all       = zeros(num_trials,nc);
best_p_all       = zeros(num_trials,nc);
time_all         = zeros(num_trials,nc);
start_nodes_all  = zeros(num_trials,nc);

%% Loop over cases
for orig_ID = 1:nc
    c = orig_ID;
    spacers          = repmat('.', 1, maxlen+3-length(cases{c}));
    W                = Graph_n_labels(c).W;
    labels           = Graph_n_labels(c).label;
    % ensure label numbering starts from 1
    if min(labels) == 0
        labels = labels + 1;
        fprintf('Labels start from 0. Adding 1.\n');
    elseif min(labels) == 1
        fprintf('Labels start from 1.\n');
    else
        fprintf('Labelling error, min(labels) = %f\n',min(labels));
    end

    % Read graph information
    n_clusters      = size(unique(Graph_n_labels(c).label),1);
    n_nodes         = size(W,1);
    n_edges         = nnz(W)/2;
    % Get random starting nodes
    start_nodes    = randi(n_nodes,num_trials,1);
    start_nodes_all(:,orig_ID) = start_nodes;

    % Create graph Laplacian and Incidence matrices
    [L,Deg,vw]  = CreateLapl(W,norm_Lap);
    % Build incidence matrix
    B           = adjacency_to_incidence(W);
    % Moore-penrose inverse, B^+
    B_pinv      = pinv(full(B));

    % % M = \beta D + L, used in the objective
    M           = beta.*Deg + L;

    % N = (\beta D + L) * B^+, used in the objective and Jacobian
    N           = M * B_pinv;

    for tr = 1:length(start_nodes)
        % Initialize RCCut and plevels
        RCCut_p_levs  = zeros(length(p_levels));
        Fscore_p_levs = zeros(length(p_levels));
        % Current starting node
        start_node    = start_nodes(tr);
        % Initialize starting vector
        s             = zeros(n_nodes,1);
        s(start_node) = 1.0;
        s             = beta .* s;

        fprintf('------------------------------------------\n');
        fprintf('Seed node: %d\n',start_node);
        fprintf('------------------------------------------\n');
        % Select furthest node to remove col from Jacobian of problem
        if isfield(Graph_n_labels(c),'coords')

            % Graph_n_labels(c).coords = Graph_n_labels(c).data;
            coords   = Graph_n_labels(c).coords;

            % Calculare pairwise distances from starting node
            vrt_dist = pdist2(coords(start_node,:),coords,'euclidean')';

            % Select furthest vertex
            [~,vrt_max]  = max(vrt_dist);
            fprintf('Furthest coords node id:%d\n',vrt_max);
        else
            % find the furthest node from start_node
            Gr    = graph(W);
            Dist  = distances(Gr,start_node);
            [~,vrt_idx] = max(Dist(:));
            [vrt_max,~] = ind2sub([n_nodes n_nodes],vrt_idx);
            fprintf('Furthest geodesic node id:%d\n',vrt_max);
        end

        % Initialize solution vector by solving a linear system
        x = lsqminnorm(N * B, s);

        % Assign a small starting value assigned to furthest geodesic
        x(vrt_max) = 1e-12;

        % Set options for the Lev-Mar optimization
        opt.tol_1       = 1.0e-7;
        opt.tol_2       = 1.0e-7;
        opt.max_iter    = 140;
        opt.vrt_max     = vrt_max;
        opt.vrt_max_val = 1e-12;
        opt.lambda_c    = 1.0e-3;

        % Function call for the NPR p-levels
        tStart_NonlinPR = tic; 
        [best_p, best_cond, best_Fscore,best_cluster,x] = ...
            NPR_levels(start_node,x,W,L,N,B,s,labels,p_levels,opt,print_level);
        best_p_all(tr,orig_ID) = best_p;
        RCCut_all(tr,orig_ID)  = best_cond;
        FSCORE_all(tr,orig_ID) = best_Fscore;
        time_all(tr,orig_ID)   = toc(tStart_NonlinPR);

    end % node trials

end % num cases

% Collect all results
for orig_ID = 1:nc
    W = Graph_n_labels(orig_ID).W;
    n_nodes     = size(W,1);
    n_edges      = nnz(W)/2;
    fprintf('=============\n');
    fprintf('%8s, nodes:%2d, edges:%2d, time(sec):%f\n',cases{orig_ID},...
        n_nodes, n_edges, sum(time_all(:,orig_ID)));
    fprintf('=============\n');
    fprintf('%s %10s %7s %10s %12s %12s \n', 'trial', 's_node', 'best_p',...
        'RCCut', 'F-score','time(sec)')
    for tr = 1:num_trials
        fprintf('%2d %10d %5f %15f %10f %10f  \n', tr, start_nodes_all(tr,orig_ID),  ...
            best_p_all(tr,orig_ID), RCCut_all(tr,orig_ID), ...
            FSCORE_all(tr,orig_ID),time_all(tr,orig_ID));
    end
end

% Print mean results
fprintf('=============\n');
fprintf('%10s\n','Mean results');
fprintf('%s %20s %10s %10s %10s %12s \n', 'Case', 'RCCut', 'F-score','std-RCCut','std-Fscore' ,'time(sec)')
fprintf('=============\n');
for orig_ID = 1:nc
    fprintf('%2s %10f %10f %10f %10f %10f \n', cases{orig_ID}, mean(RCCut_all(:,orig_ID)), ...
        mean(FSCORE_all(:,orig_ID)), std(RCCut_all(:,orig_ID)), std(FSCORE_all(:,orig_ID)),mean(time_all(:,orig_ID)));
end
fprintf('=============\n');