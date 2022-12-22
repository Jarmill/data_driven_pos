rng(44, 'twister')

n= 3;
m =2;
L = 3;

Th_vert = [ 1  1 1 1;
           -1 -1 1 1;
           1  -1 1 -1];

Th_vert = diag([1; 1; 0.7])*Th_vert + [0; 0; 0.2];
Nv = size(Th_vert, 2);

PS = possim_lpv(n, m, 0.1, L);

PS.sampler.th = @() [1;  [1; 0.7].*(2*rand(PS.L-1,1)-1) + [0; 0.2]];


sys = PS.rand_sys(1.3);

T = 30;
traj = PS.sim(T, sys);

%% perform a stabilizing control using gains u = K x

v = sdpvar(n, 1);
Y = sdpvar(m, n, Nv, 'full');

eta = 1e-3;
cons = [sum(v)==1; v>=eta];
pos_cons = [];
stab_cons = [];
Avert= cell(Nv, 1);
for i = 1:Nv
    Avert{i} = 0;
    for k = 1:L
        Avert{i} = Avert{i} + sys.A{k}*Th_vert(k, i);
    end
    stab_cons = [stab_cons; v - Avert{i}*v - sys.B*Y(:, :, i)*ones(n, 1) - eta >= 0];
    pos_cons = [pos_cons; reshape(Avert{i}*diag(v) + sys.B*Y(:, :, i), [], 1)>=0];

    
end

cons = [cons; stab_cons; pos_cons];
   
% ST = posstab(traj);

% out = ST.run();

opts = sdpsettings;
opts.robust.lplp='duality';

sol = optimize(cons,0, opts);


%% recover
if sol.problem
    disp('Infeasible')
else

    v_rec = value(v);    

    Y_rec = value(Y);
    
    K_rec = cell(Nv, 1);
    sys_clp = cell(Nv, 1);
    e_clp = zeros(n, Nv);
    for i = 1:Nv
        K_rec{i} = Y_rec(:, :, i)*diag(1./v_rec);
        sys_clp{i} = Avert{i} + sys.B*K_rec{i};
        e_clp(:, i) = abs(eig(sys_clp{i})');
    end
    
    disp('Feasible')


    
    
    
    out_clean = struct('v', v_rec, 'Y', Y_rec, 'Th_vert', Th_vert);
    out_clean.K = K_rec;
    out_clean.A = sys.A;
    out_clean.B = sys.B;
    out_clean.Avert = Avert;
        
end