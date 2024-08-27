function [best_p, best_cond, best_Fscore, best_cluster, x] = NPR_levels(start_node, x, W, L, N, B, s, labels, p_levels, opt,print_level)
    % Run Nonlinear Pagerank for different p-levels, and select the best solution best on the minimum conductance.
    % 
    % Input
    % start_node : the seed vertec,
    % x          : starting vector x. Has a small starting value at the furthest geodesic 
    % W          : adjacency matrix of the graph
    % L          : graph Laplacian
    % N          : N = (\beta D + L) * B^+, used in the objective and Jacobian computation 
    % B          : incidence matrix
    % s          : vector for the reformulation of the PR problem. See eq. 3.2 in paper.
    %             Has value \beta*1.0 at the seed node
    % labels     : labels of the nodes
    % p_levels   : the levels of p at which the optimization problem is solved
    % opt        : struct with options for the Levenberg-Marquardt optimization
    % print_level: 0/1/2 for no print/low/high verbosity
    % Output
    % best_p      : best p-level
    % best_cond   : lowest conductance
    % best_Fscore : associated F-score
    % best_cluster: cluster indicator vector
    % x           : solution vector

    % Initialize vectors to store results
    RCCut_plevs = zeros(length(p_levels), 1);
    Fscore_plevs = zeros(length(p_levels), 1);
    best_cluster_plevs = zeros(size(L, 1), length(p_levels));

    % Loop over p-levels
    for p_iter = 1:length(p_levels)

        p_curr = p_levels(p_iter);
        if print_level >= 1
            fprintf('------------------------------------------\n');
            fprintf('\nNonlinear PageRank solver at p = %.2f\n\n', p_curr);
        end
        % Get the objective function handle
        fun = @(x)nonlinPR_diff(x, N, B, s, p_curr);

        % Adapted DTU toolbox, using the modified Jacobian matrix
        opts.tau = 1e-3;
        opts.tolg = 1e-7;
        opts.tolx = 1e-7;
        opts.maxeval = 140;
        opts.vrt_max = opt.vrt_max;
        [x, info, perf] = marquardt_modJac(fun, x, opts);
        % Store the last level
        x = x(:, end);

        % Print the Levenberg-Marquardt results
        STATS = [perf.f; perf.ng; perf.mu];
        if print_level >= 2
            fprintf('%5s %15s %15s\n', 'Obj', '|grad|', 'mu');
            fprintf('%5e %5e  %5e\n', STATS);
        end
        % Sweep-cut and conductance
        norm_by_D = 0; % 0/1 for normalized/unormalized sweep cut.
        [best_cluster, best_conductance] = sweep_cut_bipart(x, W, L, norm_by_D);
        best_cluster_plevs(:, p_iter) = best_cluster;
        RCCut_plevs(p_iter) = best_conductance;

        % F-score based on best cluster
        gt_cluster = labels(start_node);
        [gt_nodes, ~] = find(labels == gt_cluster);
        sol_nodes = find(best_cluster ~= 0);
        f1_curr   = compute_f1measure(sol_nodes', gt_nodes');
        Fscore_plevs(p_iter) = f1_curr;

        if print_level >= 1
            fprintf('\nNonlinear PageRank solver at p = %.2f, RCCut = %f, Fscore = %f\n\n', ...
                p_curr, best_conductance, f1_curr);
        end
    end % p levels

    % Collect the result with the lowest conductance from all p-levels
    [best_cond, best_cond_idx] = min(RCCut_plevs);
    best_Fscore = Fscore_plevs(best_cond_idx);
    best_p = p_levels(best_cond_idx);
    best_cluster = best_cluster_plevs(:, best_cond_idx);

end %function
