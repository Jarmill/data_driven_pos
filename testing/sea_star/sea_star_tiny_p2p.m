SOLVE = 1;
SAMPLE = 1;
PLOT = 1;

rng(30, 'twister')
load('sea_star_tiny.mat', 'Sys', 'n', 'm');
% n = 3; m = 2;  T = 5; %works
% %other tests don't work. why?
% rng(60, 'twister')
% n = 5; m = 4; T = 90;

nG = n;
mG = m;

SPARSE = 0;

S_spy = true(sum(mG),sum(nG));
ncurr = 0;
mcurr = 0;
for i = 1:length(nG)

    nind = ncurr + (1:nG(i));
    mind = mcurr + (1:mG(i));
    S_spy(mind, nind) = false;

    ncurr = ncurr + nG(i);
    mcurr = mcurr + mG(i);
end



sys = struct('A', Sys.globalA, 'B', Sys.globalB);
n = size(Sys.globalA, 2);
m = size(Sys.globalB, 2);
p = size(Sys.globalC, 2);
PS = possim_cont(n, m);

%% define variables
v = sdpvar(n, 1, 'full');
delta = 1e-4;
S = sdpvar(m, n);
X = diag(v);
gam = sdpvar(1,1);

A_true = Sys.globalA;
B_true = Sys.globalB;
C = Sys.globalC;
D = Sys.globalD;
% E = ones(n);
% E = ones(n, 1);
contam = 3;
E = eye(n, sum(nG(1:contam)));

q = size(E, 2);
p2p_term = gam - (C*X + D*S)*ones(p, 1);

M = metzler_indexer(n);
Acl_true = A_true*X + B_true*S;
M_term = reshape(Acl_true(M), [], 1);
out_term = reshape(C*X+D*S, [], 1);
rej_term = E*ones(q, 1) + (A_true*X + B_true*S)*ones(n, 1);

if SPARSE
    sparse_term = S(reshape(S_spy, [], 1));
else
    sparse_term = [];
end
% cons = [v>=delta; rej_term<= -delta; out_term>=0; p2p_term>=delta; M_term >=0; S==0];

% ; reshape(A_true*X, [], 1)>=0

cons = [v>=delta; rej_term <= -delta; out_term >=0; p2p_term >= delta; ...
     M_term>=0; sparse_term==0];

opts = sdpsettings('solver', 'mosek');
sol = optimize(cons, gam, opts);

if sol.problem==0
    v_rec = value(v);
    gam_rec = value(gam);
    S_rec = value(S);
    K_rec = S_rec*diag(1./v_rec);


    Acl_rec = A_true + B_true*K_rec;

    fprintf('gam=%0.4d \n', gam_rec)
else
    disp('infeasible!')
end

% Aspy = (Sys.globalA~=0);
% Bspy = (Sys.globalB~=0);

% T = 90;
% 
% % T = 50;
% traj = PS.sample_slope(T, sys);
% 
% dopts = data_opts;
% dopts.nontrivial = 0;
