rng(40, 'twister')

n= 3;
m =2;

PS = possim(n, m);
PS.sampler.w = @() normalize(randn(PS.n,1), 1, 'norm', 2);

sys = PS.rand_sys(1.2);

T = 30;
traj = PS.sim(T, sys);

ST = posstab(traj);

out = ST.run();

%% recover and evaluate
if out.sol.problem
    disp('infeasible')
else
    sys_clp_true = sys.A + sys.B*out.K;
    eig_clp = abs(eig(sys_clp_true))'
end




