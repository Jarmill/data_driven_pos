SOLVE = 1;
SAMPLE = 0;
PLOT = 0;

%system
rng(30, 'twister')
n = 3; m = 2;  T = 12; %works


PS = possim_cont(n, m);


sys = PS.rand_sys(1.4);


%matrices
e = 1;
param = struct;
param.C = [eye(n); zeros(m, n)];
param.D = [zeros(n, m); eye(m)*0.5];
param.E = eye(n, e);
param.F = zeros(n+m, e);

traj = PS.sample_slope(T, sys);

ST = pos_p2p_cont(traj, param);

pall_true = reshape([sys.A, sys.B], [], 1);
check_poly = ST.poly.d - ST.poly.C*pall_true;

%% solve
if SOLVE
out = ST.run();

%% recover and evaluate
if ~out.sol.problem
sys_clp_true = sys.A + sys.B*out.K;
sblock = sys.A*diag(out.y) + sys.B*out.S
-ones(1, n)*sblock
%why does this have an off-diagonal element
eig_clp = eig(sys_clp_true)'
else
    disp('infeasible')

end
end