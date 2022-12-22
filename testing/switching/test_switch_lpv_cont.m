%generate the trajectory
SOLVE = 1;
SYSSAMPLE = 1;
TRAJ = 1;
PLOT = 1;

rng(40, 'twister');
n = 2;
m = 2;
L = 2;
Tmax = 35;
% Tmax = 45;
% 

% n = 3;
% m = 2;
% L = 2;
% Tmax = 35;


epsilon = 0.1;
% epsilon = 0.15; %0.15 also works, but 0.2 does not.
LS = lpvsim(n, m, L, epsilon);
% traj = LS.sim(Tmax);

Th_vert = [-1 -1 1 1;
           1  -1 1 -1];

Th_shift = [1; 0];
Th_vert = Th_vert + Th_shift;
LS.sampler.th = @() (2*rand(2, 1) - 1) + Th_shift;
       
traj = LS.sim(Tmax);

       
% Th_vert = Th_vert + [0; 0.5];
%run QMI solver
if SOLVE
LP = lpvstab_cont(traj);
% LP.const_K = 1;
out = LP.stab(Th_vert);
disp(out.sol.problem)
end

%% sample plants
if SYSSAMPLE
%% acquire samples
Nsys = 15;
sys_smp = LS.sample_sys(traj, Nsys);

%% validate samples
if ~isempty(out)
[valid, eig_out] = LP.validate_stab_multi(sys_smp, Th_vert, out.K);
end
end

%% sample a trajectory
if TRAJ
    x0 = [-2; 1.5];
    sys_all = [{traj.ground_truth}; sys_smp];

%     sys_test = sys_smp{8};

Pth = get_vertex_interp(Th_vert);
Kth = @(th) K_interp(Pth, out.K, th); %interpolating controller

Ntraj = 30;
traj_cont = cell(Nsys+1, Ntraj);
T_cont = 4;
mu_cont = 0.05;
for i = 1:Nsys+1
    rng(40, 'twister');
    for j = 1:Ntraj        
        traj_cont{i, j} = LS.sim_closed_cont(sys_all{i}, Kth, x0, T_cont, mu_cont);    
    end
end
end
%% plot the trajectory
if PLOT
figure(1)
clf
hold on
c = linspecer(2);
tiledlayout(2, 1);
nexttile;
hold on
for i = Nsys:-1:1
    for j = 1:1
        if i>1
            ccurr = c(1, :);
        else
            ccurr = c(2, :);
        end
        plot(traj_cont{i, j}.X(1, :), traj_cont{i, j}.X(2, :), 'color', ccurr, 'linewidth', 2)
    end
end
scatter(x0(1), x0(2), 300, 'ok')
xlim([x0(1)-0.1, 0.1])
xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
title('1 Parameter Sequence', 'FontSize', 16)

nexttile;
hold on
for i = 1:Nsys
    for j = 1:Ntraj
        plot(traj_cont{i, j}.X(1, :), traj_cont{i, j}.X(2, :), 'color', c(1, :))
    end
end
scatter(x0(1), x0(2), 300, 'ok')
title('30 Parameter Sequences', 'FontSize', 16)
xlim([x0(1)-0.1, 0.1])
xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
end