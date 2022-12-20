rng(44, 'twister')

n= 3;
m =2;
Nsys = 3;

PS = possim_switch(n, m, 0.1, Nsys);


sys = PS.rand_sys(1.2);

T = 30;
traj = PS.sim(T, sys);

%perform a stabilizing control using gains u = K x

v = sdpvar(n, 1);
Y = sdpvar(m, n, Nsys, 'full');

eta = 1e-3;
cons = [sum(v)==1; v>=eta];
pos_cons = [];
stab_cons = [];
for i = 1:Nsys
    stab_cons = [stab_cons; v - sys.A{i}*v - sys.B{i}*Y(:, :, i)*ones(n, 1) - eta >= 0];
    pos_cons = [pos_cons; reshape(sys.A{i}*diag(v) + sys.B{i}*Y(:, :, i), [], 1)>=0];

    
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
    
    K_rec = cell(Nsys, 1);
    sys_clp = cell(Nsys, 1);
    e_clp = zeros(n, Nsys);
    for i = 1:Nsys
        K_rec{i} = Y_rec(:, :, i)*diag(1./v_rec);
        sys_clp{i} = sys.A{i} + sys.B{i}*K_rec{i};
        e_clp(:, i) = abs(eig(sys_clp{i})');
    end
    


    
    
    
    out_clean = struct('v', v_rec, 'Y', Y_rec, 'K', K_rec,...
        'A', sys.A, 'B', sys.B);
end