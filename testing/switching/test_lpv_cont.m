rng(42, 'twister')

SOLVE = 1;
%
n= 2;
m =2;
L = 3;

T = 10; %works
% T = 20; %works

Th_vert = [ 1  1 1 1;
           -1 -1 1 1;
           1  -1 1 -1];

th_scale = [1; 1; 0.7];
th_trans = [0; 0 ;0.2];
Th_vert = diag(th_scale)*Th_vert + th_trans;
% %does not work
% rng(50, 'twister')
% n= 3;
% m =2;
% L = 3;
% 
% Th_vert = [ 1  1 1 1;
%            -1 -1 1 1;
%            1  -1 1 -1];
% Th_vert = diag([1; 0; 0.1])*Th_vert + [0; 0; 0.5];

% Th_vert = diag([0.5; 0.5])*Th_vert + [1; 1];




Nv = size(Th_vert, 2);

PS = possim_lpv_cont(n, m, 0.1, L);

PS.sampler.th = @() [1;  (th_scale(2:3)).*(2*rand(PS.L-1,1)-1) + th_trans(2:3)];

sys = PS.rand_sys(1.3);


traj = PS.sim(T, sys);

dopts = data_opts;
ST = posstab_lpv_cont_f(traj, Th_vert, dopts);

sys_combined = [];
for i = 1:L
    sys_combined = [sys_combined, sys.A{i}];
end
sys_combined = [sys_combined, sys.B];
pall_true = reshape(sys_combined, [], 1);
% pall_true = reshape([sys.A, sys.B], [], 1);
check_poly = ST.poly.d - ST.poly.C*pall_true;

Avert = cell(Nv, 1);
for i = 1:Nv
    Avert{i} = 0;
    for k = 1:L
        Avert{i} = Avert{i} + sys.A{k}*Th_vert(k, i);
    end
end


%% perform a stabilizing control using gains u = K x

%% solve
if SOLVE
out = ST.run();

%% recover and evaluate
    if ~out.sol.problem

        sys_clp_true = cell(Nv, 1);
        
        sblock = cell(Nv, 1);
        eig_clp = zeros(n, Nv);
        for i = 1:Nv        
            sys_clp_true{i} = Avert{i} + sys.B*out.K{i};
            sblock{i} = Avert{i}*diag(out.y) + sys.B*out.S(:, :, i);
            
            eig_clp(:, i) = (eig(sys_clp_true{i})');
        end
        eig_clp
    else
        disp('infeasible')        
    end
end