rng(30, 'twister')

PS = possim(3, 2);


sys = PS.rand_sys(1.4);

T = 20;
traj = PS.sim(T, sys);

ST = posstab(traj);

out = ST.stab();

%% recover and evaluate
sys_clp_true = sys.A + sys.B*out.K
eig_clp = abs(eig(sys_clp_true))'



