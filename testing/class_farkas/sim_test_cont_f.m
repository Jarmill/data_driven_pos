rng(30, 'twister')
n = 3;
m = 2;  
T = 20;

% n = 6;
% m = 8;
% T = 40;
PS = possim_cont(n, m);


sys = PS.rand_sys(1.4);


% T = 50;
traj = PS.sim(T, sys);

ST = posstab_cont_f(traj);

pall_true = reshape([sys.A, sys.B], [], 1);
check_poly = ST.poly.d - ST.poly.C*pall_true;

%% solve
out = ST.run();

%% recover and evaluate
if ~out.sol.problem
sys_clp_true = sys.A + sys.B*out.K;
sblock = sys.A*diag(out.y) + sys.B*out.S
-ones(1, n)*sblock
%why does this have an off-diagonal element
eig_clp = eig(sys_clp_true)'

end