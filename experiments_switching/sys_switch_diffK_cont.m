%test the DDC algorithm on a discrete-time switched system
SOLVE = 1;
SYSSAMPLE = 1;
TRAJ = 0;
PLOT = 0;

rng(482, 'twister')

%works
n= 3;
m =2;
Nsys = 2;

A1 = [-0.55, 0.3, 0.65; 0.06, -1.35, 0.25; 0.1, 0.15, 0.4];

B1 = [0.18, 0.08; 0.47, 0.25; 0.07, 0.95];

% A2 = [-0.5, 0.1, 0.2; 0.3 0.2 1; 0 0.5 -2];


% A2 = diag([0.05; -1; 0.5]);
A2 = diag([0; -2; 0.5])+0.1*ones(n);
B2 = [1 0; 0 0; 0 1];
% B2 = B1;

% A3 = diag([-1; 2; 1]);
% A3 = [0.5, 0.1, 0.2; 0.3 0.2 1; 0 0.5 1];
% B3 = [0.5 -0.5; 0.5 -0.5; 0 0];


sys = struct;
% sys.A = {A1; A2; A3};
% sys.B = {B1; B2; B3};
sys.A = {A1; A2};
sys.B = {B1; B2};
T = 25;
% T = 55;
% T = 90;

% n= 3;
% m =1;
% Nsys = 3;
% 
% T = 90;


epsilon = 0.1;
% epsilon = 0;
PS = possim_switch_cont(n, m, epsilon, Nsys);


% sys = PS.rand_sys(1.1);


% T = 160;
traj = PS.sample_slope(T, sys);

dopts = data_opts;
dopts.nontrivial = 1;

ST = posstab_switch_diffK_cont_f(traj, dopts);
ST_same = posstab_switch_cont_f(traj, dopts);
sys_combined = [];
for i = 1:Nsys
    sys_combined = [sys_combined, sys.A{i}, sys.B{i}];
end
pall_true = reshape(sys_combined, [], 1);
% pall_true = reshape([sys.A, sys.B], [], 1);
check_poly = ST.poly.d - ST.poly.C*pall_true;

%% solve
if SOLVE
out = ST.run();
out_same = ST_same.run();

disp([out.sol.problem, out_same.sol.problem])
%% recover and evaluate
    if ~out.sol.problem

        sys_clp_true = cell(Nsys, 1);
        sblock = cell(Nsys, 1);
        eig_clp = zeros(n, Nsys);
        for i = 1:Nsys
            sys_clp_true{i} = sys.A{i} + sys.B{i}*out.K{i};
            sblock{i} = sys.A{i}*diag(out.y) + sys.B{i}*squeeze(out.S(:, :, i));
            % -ones(1, n)*sblock
            %why does this have an off-diagonal element
            eig_clp(:, i) = (eig(sys_clp_true{i})');
        end
        eig_clp
    else
        disp('infeasible')        
    end
end
if SYSSAMPLE
% acquire samples
Nsystem= 15;
% Nsys = 1;
sys_smp = ST.sample_sys(Nsystem);

% validate samples
% if ~isempty(out)
% [valid, eig_out] = LP.validate_stab_multi(sys_smp, Th_vert, out.K);
% end
end

if TRAJ

    x0 = [0.5; 1.5; 1];
    sys_all = [{traj.ground_truth}; sys_smp];

    %     sys_test = sys_smp{8};
    
    Ntraj = 30;
%     Ntraj = 1;
%     Ntraj = 10;

    traj_cont = cell(Nsystem+1, Ntraj);
    T_cont = 30;
%     mu_cont = 0.3;
    mu_cont = 0.5;
    for i = 1:Nsystem+1
        rng(40, 'twister');
        for j = 1:Ntraj        
            [traj_cont{i, j}] = PS.sim_closed_cont(sys_all{i}, out.K, T_cont, x0, mu_cont);
        end
    end
end

%% Plot the system
if PLOT
    figure(1)
    clf
    hold on
    c = linspecer(2);
    tiledlayout(1, 2);
    nexttile;
    hold on
    for i = Nsystem:-1:1
        for j = 1:1
            if i>1
                ccurr = c(1, :);
            else
                ccurr = c(2, :);
            end
            plot3(traj_cont{i, j}.X(1, :), traj_cont{i, j}.X(2, :), traj_cont{i, j}.X(3, :), 'color', ccurr, 'linewidth', 2)
        end
    end
    scatter3(x0(1), x0(2), x0(3), 300, 'ok')

xlim([0, 1])
% ylim([0, 2])
    xlabel('$x_1$', 'interpreter', 'latex')
    ylabel('$x_2$', 'interpreter', 'latex')
    zlabel('$x_3$', 'interpreter', 'latex')
    title('1 Switching Sequence', 'FontSize', 16)
    view(3)

    nexttile;
    hold on
    for i = 1:Nsystem
        for j = 1:Ntraj
            plot3(traj_cont{i, j}.X(1, :), traj_cont{i, j}.X(2, :), traj_cont{i, j}.X(3, :), 'color', c(1, :))
        end
    end
    scatter3(x0(1), x0(2), x0(3), 300, 'ok')
    title(sprintf('%d Switching Sequences', Ntraj), 'FontSize', 16)

xlim([0, 1])
% ylim([0, 2])
xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
zlabel('$x_3$', 'interpreter', 'latex')
view(3)
end
%% plot the Lyapunov function
if PLOT
figure(6)
clf
hold on
for i = 1:Nsystem
    for j = 1:Ntraj
        vcurr = max(traj_cont{i, j}.X./out.y, [], 1);
        plot(traj_cont{i, j}.t, vcurr);
    end
end
title('Lyapunov Function along Trajectories', 'Fontsize', 16)
xlabel('$t$', 'interpreter', 'Latex')
ylabel('max$(x./v)$', 'interpreter', 'Latex')
end
