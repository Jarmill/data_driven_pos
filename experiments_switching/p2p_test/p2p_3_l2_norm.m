%discrete-time p2p


A0 = [0.7, 0, 0.2; 0.3, 0.3, 0.3; 0.2, 0.1, 0.6];
B0 = [0 -1; 0 -1; 3 1];

A1 = [0.4, 0.4, 0; 0, 0.4, 0.4; 0.6, 0.3, 0];
B1 = [1 -1; 1 -1; -2 0];

A2 = [0, 0, 1; 0, 0, 1; 0.4, 0.2, 0.2];
B2 = [0 0; 1 0; -1 1];

% A0 = [-0.2, 0.2, 0.2; 0.4, -0.7, 0.2; 0, 0.8, -3];
% B0 = [-0.4, 0.5; 0.2, -0.8; -1, 2];
% 
% A1 = [-1, 2, 0; 2, -3, 1; 1, 1, 3];
% B1 = [1 0; -2 0; -1 -1];
% 
% A2 = diag([0.8; -0.7; 1.2]) + 0.1*ones(3, 3);
% B2 = [0 1; 1 1; 1 0];
% B2 = [0 1; 0 -1; 3 -2];

Atrue = {A0, A1, A2};
Btrue = {B0, B1, B2};


n = size(A1, 1);
m = size(B1, 2);
N = length(Atrue);
M = metzler_indexer(n);


%p2p properties
C = [eye(n); zeros(2, 3)];
D = [zeros(3, 2); eye(2)];
E = [1 0; 0 0; 0 1];
F = zeros(n, 1);

SCENARIO = 2;
%0: same controller on all modes
%1: different controller on all modes
%2: periodic switching A0 -> A1 -> A2 -> A0

PRIOR_KNOWLEDGE = 1;
%0: no prior knowledge
%1: the A matrices are nonnegative
%2: the original system is used


%% sample a trajectory
epsilon = 0.25;
% epsilon = 0.5;
% epsilon = 1;
% epsilon = 0;
PS = possim_switch_cont(n, m, epsilon, N);
sys = struct;
sys.A = Atrue;
sys.B = Btrue;

%element-wise L2 bounded noise
PS.sampler.w = @() ball_sample(1, n)';

T = 80;
traj = PS.sample_slope(T, sys);


%% define uncertain system
% lam = [-10, 10; 10, -10];

Avar = sdpvar(n, n, N, 'full');
Bvar = sdpvar(n, m, N, 'full');
A = cell(3, 1);
B = cell(3, 1);

for i = 1:N
    A{i} = squeeze(Avar(:, :, i));
    B{i} = squeeze(Bvar(:, :, i));
end
rob_vars = [squeeze(reshape(Avar, [], 1, 1)); squeeze(reshape(Bvar, [], 1, 1))];
rob_con = [uncertain(rob_vars)];


traj_con = [];
res_list = [];
for t = 1:T
    Acurr = squeeze(A{ traj.S(t)});
    Bcurr = squeeze(B{traj.S(t)});
    res_curr = traj.Xdelta(:, t) - Acurr*traj.Xn(:, t) - Bcurr*traj.U(:, t);
    
    %check that the residuals are valid
    res_sq = replace(res_curr, [Acurr, Bcurr], [Atrue{traj.S(t)}, Btrue{traj.S(t)}]);
    res_list = [res_list; sqrt(sum(res_sq.^2))];
    % traj_con = [traj_con; sum(res_curr.^2) <= epsilon^2];
    
    traj_con = [traj_con; norm(res_curr, 2) <= epsilon];
end

%check the residuals



%the original system is positive (with zero input)
prior_con = [];
for i = 1:N
    if PRIOR_KNOWLEDGE==2
        prior_con = [prior_con; A{i} == Atrue{i}];
        prior_con = [prior_con; B{i}== Btrue{i}];
        traj_con = [];
    
    elseif PRIOR_KNOWLEDGE==1    
        Avec = reshape(A{i}, [], 1);        
        prior_con = [prior_con; Avec >= 0];
    end
end

%% optimization (lyapunov)

v = sdpvar(n, N, 'full');
S = sdpvar(m, n, N);

gam = sdpvar(1, 1);

delta = 1e-3;
stab_con = [];
pos_con = [];

for i = 1:N
    if SCENARIO == 2
        vnext = v(:, mod(i, N)+1);
    else
        vnext = v(:, i);
    end

    stab_con_curr = -vnext+ sum(E, 2) + A{i}*v(:, i) + B{i}*squeeze(S(:, :, i))*ones(n, 1);
    stab_con = [stab_con; stab_con_curr];
end



if SCENARIO < 2

    vcon = [];
    for i = 1:(N-1)
        vcon = [vcon; v(:, i) == v(:, i+1)];
    end
else
    vcon = [];
end

for i = 1:N
    pos_con_curr = reshape(A{i}*diag(v(:, i)) + B{i}*squeeze(S(:, :, i)), [], 1);

    pos_out_curr = gam - (C*diag(v(:, i)) + D*squeeze(S(:, :, i)))*ones(3, 1);
    pos_out_curr_x = reshape((C*diag(v(:, i)) + D*squeeze(S(:, :, i))), [], 1);
    pos_con = [pos_con; pos_con_curr; pos_out_curr; pos_out_curr_x];
    

end

K_con = [];
if SCENARIO == 0
    for i = 1:N-1
        K_con = [K_con; S(:, :, i) == S(:, :, i+1)];
    end
end


cons = [prior_con; vcon; stab_con <= (-delta); reshape(v, [], 1) >= delta; pos_con>=0; K_con; rob_con; traj_con];

sol = optimize(cons, gam);

% z = {value(z1), value(z2)};
v_rec = value(v);
S_rec = value(S);
K_rec = cell(N, 1);
Acl_rec = cell(N, 1);
gam_rec = value(gam);

for i = 1:N
    Xi_rec = diag(1./v_rec(:, i));
    K_rec{i} = S_rec(:, :, i)*Xi_rec;
    Acl_rec{i} = Atrue{i}+Btrue{i}*K_rec{i};
end


if sol.problem==0
    % abs(eig(Acl_rec{1}))'
    gam_rec
else
    disp('infeasible')
end
