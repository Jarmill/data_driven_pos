%test the DDC algorithm on a discrete-time switched system
%the controller is aware of the current system state, and can select a
%subsystem-dependent control gain. However, there is a common Lyapunov
%function
SOLVE = 1;
rng(44, 'twister')

%works
% n= 3;
% m =2;
% Nsys = 2;
% T = 60;

%works
n= 3;
m =2;
Nsys = 3;
T = 80;

% 
% n= 3;
% m =2;
% Nsys = 4;
% T = 120;



% %works
% n = 4;
% m = 2;
% Nsys = 2;
% T = 140;
% 
% %works
% n= 2;
% m =2;
% Nsys = 4;
% T = 115;

epsilon = 0.1;
% epsilon = 0;
PS = possim_switch(n, m, epsilon, Nsys);


sys = PS.rand_sys(1.2);


% T = 160;
traj = PS.sim(T, sys);

ST = posstab_switch_diffK_f(traj);

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

%% recover and evaluate
    if ~out.sol.problem

        sys_clp_true = cell(Nsys, 1);
        sblock = cell(Nsys, 1);
        eig_clp = zeros(n, Nsys);
        for i = 1:Nsys
            sys_clp_true{i} = sys.A{i} + sys.B{i}*out.K{i};
            sblock{i} = sys.A{i}*diag(out.y) + sys.B{i}*out.S(:, :, i);
            % -ones(1, n)*sblock
            %why does this have an off-diagonal element
            eig_clp(:, i) = abs(eig(sys_clp_true{i})');
        end
        eig_clp
    else
        disp('infeasible')        
    end
end
