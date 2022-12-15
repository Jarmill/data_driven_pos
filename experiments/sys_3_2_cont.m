SOLVE = 1;
SAMPLE = 1;
PLOT = 1;

rng(30, 'twister')
n = 3; m = 2;  T = 5; %works

% %other tests don't work. why?
% rng(60, 'twister')
% n = 5; m = 4; T = 90;


PS = possim_cont(n, m);


% sys = PS.rand_sys(1.4);

A = [-0.55, 0.3, 0.65; 0.06, -1.35, 0.25; 0.1, 0.15, 0.4];

B = [0.18, 0.08; 0.47, 0.25; 0.07, 0.95];

sys = struct('A', A, 'B', B);


% T = 50;
traj = PS.sample_slope(T, sys);

ST = posstab_cont_f(traj);

pall_true = reshape([sys.A, sys.B], [], 1);
check_poly = ST.poly.d - ST.poly.C*pall_true;

%% solve
if SOLVE
out = ST.run();

%% recover and evaluate
if ~out.sol.problem
sys_clp_true = sys.A + sys.B*out.K;
sblock = sys.A*diag(out.y) + sys.B*out.S
% -ones(1, n)*sblock
%why does this have an off-diagonal element
eig_clp = eig(sys_clp_true)'
else
    disp('infeasible')

end
end


%% Sample trajectories

if SAMPLE&&~out.sol.problem
%     Nsys = 30;
    Nsys = 100;
    x0 = ones(n, 1);
    
    sys_vec = cprnd(Nsys, ST.poly.C, ST.poly.d);
    
    PS.epsilon = 0;
    PS.sampler.u = @(x) out.K*x;
    
    traj_smp = cell(Nsys, 1);
    Tsim = 20;
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
    scatter3(x0(1), x0(2), x0(3), 300, 'k')
    hold off
    

    
    xlabel('$x_1$', 'interpreter', 'latex')
    ylabel('$x_2$', 'interpreter', 'latex')
    zlabel('$x_3$', 'interpreter', 'latex')
    title(sprintf('Positive System Control (Nsys = %d)', Nsys), 'fontsize', 16)
        


%% analyze the Lyapunov function
%the CLF we get is max(x./y), not y'x. 
figure(5)
clf
vlin = cell(length(traj_smp), 1);
vmax = cell(length(traj_smp), 1);

for i = 1:length(traj_smp)
    tc = traj_smp{i};
    vlin{i} = out.y'*tc.Xn;
    vmax{i} = max(tc.Xn./out.y, [], 1);
end

subplot(1, 2, 1)
hold on
for i = 1:length(traj_smp)
    plot(traj_smp{i}.t, vlin{i});
end
title('Linear Lyap', 'Fontsize', 16)
xlabel('t')
ylabel('V(x)')

subplot(1, 2, 2)
hold on
for i = 1:length(traj_smp)
    plot(traj_smp{i}.t, vmax{i});
end
title('Max Lyap', 'Fontsize', 16)
xlabel('t')
ylabel('V(x)')
end