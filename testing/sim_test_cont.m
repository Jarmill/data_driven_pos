rng(30, 'twister')

PS = possim_cont(3, 2);


sys = PS.rand_sys(1.4);

T = 20;
traj = PS.sim(T, sys);

ST = posstab_cont(traj);

out = ST.stab();

%% recover and evaluate
sys_clp_true = sys.A + sys.B*out.K
skew_clp_true = sys.A*diag(out.y) + sys.B*out.S;
%why does this have an off-diagonal element
eig_clp = abs(eig(sys_clp_true))'

pall_true = reshape([sys.A, sys.B], [], 1);
check_poly = ST.poly.d - ST.poly.C*pall_true;
