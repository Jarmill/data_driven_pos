%discrete-time switching-based stabilization in which the control signal
%selects which subsystem should be active

%uses Lyapunov-Metzler Inequalities
rng(44, 'twister')

% n= 3;
% m =2;
% Nsys = 3;

n= 2;
m =1;
Nsys = 2;



PS = possim_switch(n, m, 0.1, Nsys);


sys.A = {[0.5, 0; 0, 1.3], [1.3, 0; 0, 0.5]};
sys.B = {[0; 1], [1; 0]};

% sys = PS.rand_sys(1.2);

T = 30;
traj = PS.sim(T, sys);

%% perform a stabilizing control using subsystem selection

Lam = ones(n, n)/n; %a nonnegative matrix satisfying Lam*1 = 1 (markov)

v = sdpvar(n, n, 'full');
Y = sdpvar(m, n, Nsys);


lm_expr = zeros(Nsys, n, 'like', sdpvar);
for i = 1:Nsys
    for j = 1:n
        if j ~= i
            lm_expr(:, i) = lm_expr(:, i) + Lam(i, j)*(v(:, j)- v(:, i));
        end
    end
end

% end

eta = 1e-3;
cons = [reshape(v, [], 1)>=eta];
pos_cons = [];
stab_cons = [];
M = metzler_indexer(n);
for i = 1:Nsys   

%     LM_curr = 
    stab_cons = [stab_cons; -lm_expr(:, i)- sys.A{i}*v(:, i)  - eta >= 0];    
    Acl_vec = reshape(sys.A{i}*diag(v(:, i)) + sys.B{i}*Y(:, :, i), [], 1);
    pos_cons = [pos_cons;Acl_vec(M)>=0];   
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

%     Y_rec = value(Y);
%     
%     K_rec = cell(Nsys, 1);
%     sys_clp = cell(Nsys, 1);
%     e_clp = zeros(n, Nsys);
%     for i = 1:Nsys
%         K_rec{i} = Y_rec(:, :, i)*diag(1./v_rec(:, i));
%         sys_clp{i} = sys.A{i} + sys.B{i}*K_rec{i};
%         e_clp(:, i) = abs(eig(sys_clp{i})');
%     end
%     
   
    out_clean = struct('v', v_rec, 'A', sys.A);
    
    
%     out_clean = struct('v', v_rec, 'Y', Y_rec, 'K', K_rec,...
%         'A', sys.A, 'B', sys.B);
end