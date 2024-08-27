function [func, jac] = nonlinPR_diff(x,N,B,s,p_curr)
% Objective and jacobian evaluation for the Nonlinear modified Pagerank
% Input
% x: solution vector
% N: (\beta D + L) * B^+, used in the objective and Jacobian  
% B: Incidence matrix
% s: \beta*s, a vector heavily skewed on the starting vertex
% p_curr: the current level of p
% Output
% func: objective evaluation
% jac: jacobian evaluation

% % The full objective reads
% % f(x) = B^+ (|Bx|^{p-2} \otimes Bx)

% Power |Bx|^q, and q/2
q      = p_curr - 2;
q_half = q/2;

% Compute vector ε + B * x, and select ε
Bx      = B * x;
Bx_abs  = abs(Bx);
% Fixed value of \epsilon throughout experiments
epsilon = 1e-11;  % all small- and mid-scale cases
% epsilon = 1e-6; % only large Gauss and USPS
print_statistics = 0;
if print_statistics
    fprintf('Median of vector |Bx| is %e, mean is %e and epsilon=%e \n',median(Bx_abs), mean(Bx_abs),epsilon);
end

% Compute (ε+|Bx|^2)^{(p-2)/2}
% Bx_q = abs(Bx).^q;
Bx_q = (epsilon+abs(Bx).^2).^q_half;


% The full objective
% func = N*(Bx_q.*(ε +Bx)) - s;
func = N*(Bx_q.*(epsilon + Bx)) - s;

% Power |Bx|^r
r = p_curr - 4;

% Compute vector w = (p-2) * (|Bx|^r \odot Bx) 
t_vec = (epsilon+Bx.^2).^(0.5);

% Sparse matrix multiplication (tmp = Bx .* w + Bx_q;)
v_diag = t_vec.^q + q * (Bx.^2) .* (t_vec.^r); 

% Build a sparse Jacobian
num_elem = length(v_diag);
tmp_mat  = sparse(1:num_elem, 1:num_elem,v_diag, num_elem, num_elem);
jac      = N * (tmp_mat * B);

end % FUNCTION