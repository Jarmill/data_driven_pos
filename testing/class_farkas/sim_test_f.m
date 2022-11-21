rng(40, 'twister')

% n= 3;
% m =3;
% T = 30;
n = 7;
m = 7;
T = 100;

PS = possim(n, m);


sys = PS.rand_sys(1.2);


traj = PS.sim(T, sys);

ST = posstab(traj);
% [vars] = ST.make_vars()
% po = ST.poly_stab(vars)
% [cons, vars] = ST.make_program()
out = ST.run();

ST_f = posstab_f(traj);
% [vars] = ST.make_vars()
% po = ST.poly_stab(vars)
% [cons, vars] = ST.make_program()
out_f = ST_f.run();



% 
%% recover and evaluate
if out_f.sol.problem
    disp('infeasible')
else
    sys_clp_true = sys.A + sys.B*out.K;
    eig_clp = abs(eig(sys_clp_true))';
    sys_clp_true_f = sys.A + sys.B*out_f.K;
    eig_clp_f = abs(eig(sys_clp_true_f))';
    
    [eig_clp; eig_clp_f]
end




