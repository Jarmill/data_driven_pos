rng(40, 'twister')

% n= 3;
% m =3;
% T = 30;

n= 2;
m =2;
T = 12;

% n = 7;
% m = 7;
% T = 100;

PS = possim(n, m);


sys = PS.rand_sys(1.2, 1);
dopts = data_opts;
dopts.nontrivial = 1;
dopts.pos_A = 0;
dopts.pos_B = 0;

pall = reshape([sys.A, sys.B],[],1);

traj = PS.sim(T, sys);

%choose the parameters
p = 2; %number of external inputs
q = n+m;
param = struct;
param.C = [eye(n);zeros(m, n)];
param.D = [zeros(n, m); eye(m)];
param.E = eye(n, p);
param.F = zeros(q, p);


ST_f = posstab_f(traj, dopts);
% 
PT_f = pos_p2p_f(traj, param, dopts);


poly_check = PT_f.poly.d - PT_f.poly.C*pall;
%% run
out_f_stab = ST_f.run();
out_f = PT_f.run();


% 
%% recover and evaluate
if out_f.sol.problem
    disp('infeasible')
else
    sys_clp_true_f = sys.A + sys.B*out_f.K
    eig_clp_f = abs(eig(sys_clp_true_f))'
    out_f.gamma
end

