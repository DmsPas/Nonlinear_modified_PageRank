function [h, mu] = geth_mod(A, J, r, g, mu)
    % Solve  (A + mu*I)h = -g  with possible adjustment of  mu

    % Version 10.02.17.

    % This file is part of the "immoptibox" package and is distributed under
    % the 3-Clause BSD License. A separate license file should be provided as
    % part of the package.
    %
    % Adapted and modified for the Nonlinear modified PageRank method.

    mA = max(abs(A(:)));

    if mA == 0 % zero matrix
        h = g / (-mu);
        return
    end

    % fprintf('In vanilla LM from immoptibox\n');
    % Factorize with check of pos. def.
    n = size(A, 1);
    chp = 1;

    while chp
        [R, chp] = chol(A + mu * speye(n));

        if chp == 0 % check for near singularity
            chp = rcond(R) < 1e-15;
        end

        if chp
            mu = max(10 * mu, eps * mA);
        end

    end

    % Solve  (R'*R)h = -g
    h = R \ (R' \ (-g));

end %function
