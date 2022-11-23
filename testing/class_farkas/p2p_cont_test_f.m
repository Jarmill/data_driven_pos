rng(40, 'twister')

n= 3;
m =2;
T = 20;  % p2p of 6.4539
%T = 30; % p2p of 5.0182
% T = 50;  % p2p of 4.4967
% T = 80;  % p2p of  4.0619
% T = 110; % p2p of  4.0252
% T = 140;   % p2p of  3.9825

% n = 7;
% m = 7;
% T = 100;

PS = possim_cont(n, m);


% sys = PS.rand_sys(1.2, 1);
%this system remains infeasible for p2p, but can be stabilized.
%I will need to find a system to use for continuous-time p2p testing.
sys = struct('A', [-0.2, 0.2, 0.2; 0.4, -0.7, 0.2; 0, 0.8, -3],...
    'B', [-0.4, 0.5; 0.2, -0.8; -1, 2]);

dopts = data_opts;
dopts.nontrivial = 1;
dopts.pos_A = 1;
dopts.pos_B = 0;

pall = reshape([sys.A, sys.B],[],1);

traj = PS.sample_slope(T, sys);

%choose the parameters
p = 2; %number of external inputs
q = n+m;
param = struct;
param.C = [eye(n);zeros(m, n)];
param.D = [zeros(n, m); eye(m)];
param.E = eye(n, p);
param.F = zeros(q, p);


ST_f = posstab_cont_f(traj, dopts);
% 
PT_f = pos_p2p_cont_f(traj, param, dopts);


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
    eig_clp_f = eig(sys_clp_true_f)'
    out_f.gamma
    
    out_f.K
end

