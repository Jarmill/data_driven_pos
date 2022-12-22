%test the DDC algorithm on a discrete-time switched system
SOLVE = 0;
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

ST = posstab_lpv_f(traj);

sys_combined = [];
for i = 1:Nsys
    sys_combined = [sys_combined, sys.A{i}];
end
sys_combined = [sys_combined, sys.B];
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
            sys_clp_true{i} = sys.A{i} + sys.B{i}*out.K;
            sblock{i} = sys.A{i}*diag(out.y) + sys.B{i}*out.S;
            % -ones(1, n)*sblock
            %why does this have an off-diagonal element
            eig_clp(:, i) = abs(eig(sys_clp_true{i})');
        end
        eig_clp
    else
        disp('infeasible')        
    end
end
