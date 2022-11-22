%Example 1 of https://ieeexplore.ieee.org/document/8618689
% A = diag(-[1,2,3,4]);
A = diag(-[1,2,3,4]);
B = [-1 1 0 0 0 0; 0 -1 -1 1 0 0; 1 0 1 -1 -1 1; 0 0 0 0 1 -1];

n = 4; 
m = 6;


%% Declare Variables

delta = 1e-3;

%design
y = sdpvar(n, 1);
% S = sdpvar(m, n);
s = sdpvar(m, 1);
S = blkdiag(s(1), s(2:3), s(4:5), s(6));
gamma = sdpvar(1, 1);
%parameters of p2p problem

p = 4; %number of external inputs
% q=n+m; %number of controlled outputs
q = n+m;
C = [eye(n);zeros(m, n)];
D = [zeros(n, m); eye(m)];
% D = 0.1*ones(n, m);
% C = [eye(n); zeros(m, n)]; %like an LQR setup for (C, D)
% C = rand(q, n);
% D = zeros(q, m);
% D = [zeros(n, m); eye(m)*0.5];
E = eye(n, p);
F = zeros(q, p);

%% form constraints

%basic constraints
% cons = [y >= delta; sum(y)==1];
cons = [y>=delta];
X = diag(y);

%condition (6.2) of https://ieeexplore.ieee.org/document/8618689
%equation (12) of https://arxiv.org/pdf/1204.3554.pdf
Abar = (A*X + B*S);
Cbar = (C*X + D*S);
term_x = -Abar *ones(n, 1) - E*ones(p, 1);
term_z = gamma*ones(q, 1) - Cbar*ones(n, 1) - F*ones(p, 1);

poscon_x =  reshape((1-eye(n)).*Abar, [], 1);
poscon_z = reshape(Cbar, [], 1);
% poscon_z = 0;
cons = [cons; [term_x; term_z] >= delta; [poscon_x; poscon_z] >= 0];

%% solve program
opts = sdpsettings;

sol = optimize(cons,gamma, opts);


%% recover
if sol.problem
    disp('Infeasible')
else

    y_rec = value(y);
    v_rec = 1./y_rec;

    S_rec = value(S);
    K_rec = S_rec*diag(v_rec);


    sys_clp = A + B*K_rec
    e_clp = eig(sys_clp)'
    gamma_rec = value(gamma)
    
    out_clean = struct('y', y_rec, 'v', v_rec, 'S', S_rec, 'K', K_rec,...
        'A', A, 'B', B, 'gamma', gamma);
end