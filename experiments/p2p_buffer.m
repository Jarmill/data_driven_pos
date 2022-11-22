%Example 1 of https://ieeexplore.ieee.org/document/8618689

SOLVE = 1;
SAMPLE = 0;
PLOT = 0;

A = diag(-[1,2,3,4]);
B = [-1 1 0 0 0 0; 0 -1 -1 1 0 0; 1 0 1 -1 -1 1; 0 0 0 0 1 -1];

n = 4; 
m = 6;

%generate data

PS = possim_cont(n, m);
sys = struct('A', A, 'B', B);


T = 50;
traj = PS.sample_slope(T, sys);

%matrices
e = 1;
param = struct;
param.C = [eye(n); zeros(m, n)];
param.D = [zeros(n, m); eye(m)*0.5];
param.E = eye(n, e);
param.F = zeros(n+m, e);

ST = pos_p2p_cont(traj, param);

pall_true = reshape([sys.A, sys.B], [], 1);
check_poly = ST.poly.d - ST.poly.C*pall_true;

