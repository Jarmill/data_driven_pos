SOLVE = 1;
SAMPLE = 0;
PLOT = 0;

rng(30, 'twister')
n = 3; m = 2;  T = 12; %works

% %other tests don't work. why?
% rng(60, 'twister')
% n = 5; m = 4; T = 90;


PS = possim_cont(n, m);


sys = PS.rand_sys(1.4);


% T = 50;
traj = PS.sample_slope(T, sys);

ST = posstab_cont(traj);

pall_true = reshape([sys.A, sys.B], [], 1);
check_poly = ST.poly.d - ST.poly.C*pall_true;

%% solve
if SOLVE
out = ST.run();

%% recover and evaluate
if ~out.sol.problem
sys_clp_true = sys.A + sys.B*out.K;
sblock = sys.A*diag(out.y) + sys.B*out.S
-ones(1, n)*sblock
%why does this have an off-diagonal element
eig_clp = eig(sys_clp_true)'
else
    disp('infeasible')

end
end


%% Sample trajectories

if SAMPLE&&~out.sol.problem
    Nsys = 30;
    x0 = ones(n, 1);
    
    sys_vec = cprnd(Nsys, ST.poly.C, ST.poly.d);
    
    PS.epsilon = 0;
    PS.sampler.u = @(x) out.K*x;
    
    traj_smp = cell(Nsys, 1);
    Tsim = 10;
    for i = 1:Nsys
        sys_mat = reshape(sys_vec(i, :), n, n+m);
        sys_curr = struct('A', sys_mat(:, 1:n), 'B', sys_mat(:, (n+1):end));
%         traj_smp{i} = PS.sim(T, sys_curr, x0);
        
        [tcurr, xcurr] = ode45(@(t, x) (sys_curr.A + sys_curr.B*out.K)*x, [0, Tsim], x0);

        traj_smp{i} = struct('t', tcurr, 'Xn', xcurr');
    end
end

%% Plot results

if PLOT&&~out.sol.problem
    figure(1)
    clf
    hold on
    for i = 1:Nsys
        plot3(traj_smp{i}.Xn(1, :), traj_smp{i}.Xn(2, :), traj_smp{i}.Xn(3, :), 'color', [0.9153    0.2816    0.2878]);
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