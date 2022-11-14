function [A_red, b_red] = nontrivial_constraints(A, b, tol)
%NONTRIVIAL_CONSTRAINTS Identify nontrivial linear constraints of the
%polytope A x <= b
%Perform this by iterative linear programming in yalmip
%may require multiple iterations (numeric difficulties)


if nargin < 3
    tol = 1e-8;
end

[m, d] = size(A);

w = sdpvar(d, 1);
cost = sdpvar(d, 1);

objective = cost'*w;
cons = A*w <= b;


linopts = sdpsettings('solver', 'mosek',...
                'mosek.MSK_DPAR_INTPNT_TOL_DFEAS', 1e-10,...
                'mosek.MSK_DPAR_INTPNT_TOL_PFEAS', 1e-10);
Lopt = optimizer(cons, objective, linopts, cost, w);


nontrivial = logical(zeros(m, 1));

% tic
% tol = 1e-8;
for i = 1:m
    %solve linprog in current direction
    if ~nontrivial(i)
        [w_curr, error_optA, statusA, dual_curr] = Lopt(full(-A(i, :)'));
        
        active_cons = dual_curr{1} > tol;
%         n_active = nnz(active_cons);
        nontrivial(active_cons) = true;
    end
end

% toc

% nnz(nontrivial)

A_red = A(nontrivial, :);
b_red = b(nontrivial, :);


end
