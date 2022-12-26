rng(42, 'twister')

SOLVE = 1;
SYSSAMPLE = 1;
TRAJ= 1;
PLOT = 1;
%
n= 2;
m =2;
L = 3;

T = 10; %works

Th_vert = [ 1  1 1 1;
           -1 -1 1 1;
           1  -1 1 -1];

th_scale = [1; 1; 0.7];
th_trans = [0; 0 ;0.2];
Th_vert = diag(th_scale)*Th_vert + th_trans;


Nv = size(Th_vert, 2);

%% sample a trajectory 
PS = possim_lpv_cont(n, m, 0.1, L);

PS.sampler.th = @() [1;  (th_scale(2:3)).*(2*rand(PS.L-1,1)-1) + th_trans(2:3)];

sys = PS.rand_sys(1.3);


traj = PS.sim(T, sys);

dopts = data_opts;
ST = posstab_lpv_cont_f(traj, Th_vert, dopts);

%check if the plant is inside the polytope (bug-check)
sys_combined = [];
for i = 1:L
    sys_combined = [sys_combined, sys.A{i}];
end
sys_combined = [sys_combined, sys.B];
pall_true = reshape(sys_combined, [], 1);
% pall_true = reshape([sys.A, sys.B], [], 1);
check_poly = ST.poly.d - ST.poly.C*pall_true;

Avert = cell(Nv, 1);
eig_open = zeros(n, Nv);
for i = 1:Nv
    Avert{i} = 0;
    for k = 1:L
        Avert{i} = Avert{i} + sys.A{k}*Th_vert(k, i);
    end
    eig_open(:, i) = eig(Avert{i});
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

%% sample plants
if SYSSAMPLE
% acquire samples
Nsys = 15;
sys_smp = ST.sample_sys(Nsys);

% validate samples
% if ~isempty(out)
% [valid, eig_out] = LP.validate_stab_multi(sys_smp, Th_vert, out.K);
% end
end

%% sample a trajectory
if TRAJ

        x0 = [0.5; 1.5];
    sys_all = [{traj.ground_truth}; sys_smp];

    %     sys_test = sys_smp{8};
    
    Ntraj = 30;

    traj_cont = cell(Nsys+1, Ntraj);
    T_cont = 6;
    mu_cont = 0.3;
    for i = 1:Nsys+1
        rng(40, 'twister');
        for j = 1:Ntraj        
            [traj_cont{i, j}, Kth_tar] = PS.sim_closed_cont(Th_vert, out.K, sys_all{i}, x0, T_cont, mu_cont);
        end
    end
end

%% Plot the system
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

    xlim([0, 0.65])
ylim([0, 1.75])
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
title(sprintf('%d Parameter Sequences', Ntraj), 'FontSize', 16)

xlim([0, 0.65])
ylim([0, 1.75])
xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
end
%% plot the Lyapunov function
if PLOT
figure(6)
clf
hold on
for i = 1:Nsys
    for j = 1:Nsys
        vcurr = max(traj_cont{i, j}.X./out.y, [], 1);
        plot(traj_cont{i, j}.t, vcurr);
    end
end
title('Lyapunov Function along Trajectories', 'Fontsize', 16)
xlabel('$t$', 'interpreter', 'Latex')
ylabel('max$(x./v)$', 'interpreter', 'Latex')
end


