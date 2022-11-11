SOLVE = 1;
SAMPLE = 1;
PLOT = 1;
n = 3;
m=3;

rng(40, 'twister')

PS = possim(n, m);


sys = PS.rand_sys(1.4);

T = 20;
traj = PS.sim(T, sys);

%% Solve system
if SOLVE
ST = posstab(traj);

out = ST.stab();

% recover and evaluate
sys_clp_true = sys.A + sys.B*out.K
eig_clp = abs(eig(sys_clp_true))'
end

%% Sample trajectories

if SAMPLE
    Nsys = 30;
    x0 = ones(n, 1);
    
    sys_vec = cprnd(Nsys, ST.poly.C, ST.poly.d);
    
    PS.epsilon = 0;
    PS.sampler.u = @(x) out.K*x;
    
    traj_smp = cell(Nsys, 1);
    for i = 1:Nsys
        sys_mat = reshape(sys_vec(i, :), n, n+m);
        sys_curr = struct('A', sys_mat(:, 1:n), 'B', sys_mat(:, (n+1):end));
        traj_smp{i} = PS.sim(T, sys_curr, x0);
    end
end

%% Plot results

if PLOT
    figure(1)
    clf
    hold on
    for i = 1:Nsys
        scatter3(traj_smp{i}.Xn(1, :), traj_smp{i}.Xn(2, :), traj_smp{i}.Xn(3, :), [], 'k', 'filled');
    end
    scatter3(x0(1), x0(2), x0(3), 300, 'b')
    hold off
    
    xlim([0, 1.2])
    ylim([0, 1.2])
    zlim([0, 1.2])
    
    xlabel('x_1')
    ylabel('x_2')
    zlabel('x_3')
    title(sprintf('Positive System Control (Nsys = %d)', Nsys), 'fontsize', 16)
        
end


