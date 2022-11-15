
%estimate the peak-to-peak gain of this stable system
n = 3;

A = [-0.519976096997155,0.278447513676141,1.140045359795860;0.190031751239087,-1.373949725354972,1.126891575522600;0.452036570909624,0.851505132275759,-7.137876097704595];

%system definition
%xdot = Ax + Ew
%z    = Cx + Fw
p = 1; %number of external inputs
q = n;
C = eye(n);
E = eye(n, p);
F = zeros(q, p);

xi = sdpvar(n, 1);
gamma = sdpvar(1, 1);

poscon = [zeros(n, 1); gamma*ones(n, 1)] - [A, E; C, F]*[xi; ones(p, 1)];
delta = 1e-3;
% cons  = ([poscon; xi] >= delta);
cons  = ([poscon] >= delta);
cons = [cons; xi >= delta];
sol = optimize(cons, gamma);

if ~sol.problem
    xi_rec = value(xi)
    gamma_rec = value(gamma)
else
    disp('infeasible');
end
