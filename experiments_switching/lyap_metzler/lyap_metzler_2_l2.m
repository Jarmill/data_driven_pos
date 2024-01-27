%data-driven stabilization of choosing the switching sequence
rng(55, "twister");
A1 = [2, 1; 1, 2];
A2 = [1, 0; 2, 2];

B1 = [1; 0];
B2 = [0; 1];


Atrue = {A1, A2};
Btrue = {B1, B2};

n = size(A1, 1);
m = size(B1, 2);
N = length(Atrue);

M = metzler_indexer(n);

PLOT = 1;
%% sample a trajectory
% epsilon = 0.05;
% epsilon = 1e-5;
% epsilon = 0;
% epsilon = 0.25;
% epsilon = 0.5;
epsilon = 1;
% epsilon = 0.01;
PS = possim_switch_cont(n, m, epsilon, N);
sys = struct;
sys.A = Atrue;
sys.B = Btrue;

%element-wise L2 bounded noise

NORM = 2;
if NORM == 2
    PS.sampler.w = @() ball_sample(1, n)';
end

% T = 120;
T = 15;
traj = PS.sample_slope(T, sys);


%% set up uncertain variables

A = sdpvar(n, n, N, 'full');
B = sdpvar(n, m, N, 'full');
rob_vars = [squeeze(reshape(A, [], 1, 1)); squeeze(reshape(B, [], 1, 1))];
rob_con = [uncertain(rob_vars)];

traj_con = [];
res_sq_list = [];
for t = 1:T
    Acurr = squeeze(A(:, :, traj.S(t)));
    Bcurr = squeeze(B(:, :, traj.S(t)));
    res_curr = traj.Xdelta(:, t) - Acurr*traj.Xn(:, t) - Bcurr*traj.U(:, t);
    res_sq_list = [res_sq_list; sum(res_curr.^2)];
    % traj_con = [traj_con; sum(res_curr.^2) <= epsilon^2];
    traj_con = [traj_con; norm(res_curr, NORM) <= epsilon];
end

%the original system is positive
prior_con = [];
KNOWN_SYS = 2;
for i = 1:N
    if KNOWN_SYS==1
        prior_con = [prior_con; A(:, :, i) == Atrue{i}];
        prior_con = [prior_con; B(:, :, i) == Btrue{i}];
    elseif KNOWN_SYS == 2
        prior_con = [prior_con; B(:, :, i) == Btrue{i}];
        Avec = reshape(squeeze(A(:, :, i)), [], 1);
        prior_con = [prior_con; Avec(M) >= 0];
    else    
        Avec = reshape(squeeze(A(:, :, i)), [], 1);
        Bvec = reshape(squeeze(B(:, :, i)), [], 1);
        prior_con = [prior_con; Bvec >= 0; Avec(M) >= 0];
    end
end
% prior_con = [B>=0; A>=0; ]

%% perform optimization problem

% for i = 1:T
% end

v = sdpvar(n, N, 'full');
S = sdpvar(m, n, N);
Lam = ones(n, N)*10; %a nonnegative matrix satisfying Lam*1 = 1 (markov)
% Lam = zeros(n, N);

lm_expr = zeros(n, N, 'like', sdpvar);
for i = 1:N
    for j = 1:N
        if j ~= i
            lm_expr(:, i) = lm_expr(:, i) + Lam(i, j)*(v(:, j)- v(:, i));
        end
    end
end

stab_con = [];
pos_con = [];
delta = 1e-4;

for i = 1:N
    Acl_curr = A(:, :, i)*diag(v(:, i)) + B(:, :, i)*squeeze(S(:, :, i));
    stab_con_curr = Acl_curr*ones(n, 1)+ lm_expr(:, i);
    stab_con = [stab_con; stab_con_curr];

    
    pos_con_curr = reshape(Acl_curr, [], 1);

    pos_con = [pos_con; pos_con_curr(M)];
end


cons = [stab_con <= (-delta); reshape(v, [], 1) >= delta; pos_con>=0];

cons = [cons; rob_con; traj_con; prior_con];

sol = optimize(cons, []);

if sol.problem==0
% z = {value(z1), value(z2)};
v_rec = value(v);
S_rec = value(S);
K_rec = cell(N, 1);
Acl_rec = cell(N, 1);
e_rec = zeros(N, n);

for i = 1:N
    Xi_rec = diag(1./v_rec(:, i));
    K_rec{i} = S_rec(:, :, i)*Xi_rec;
    Acl_rec{i} = Atrue{i}+Btrue{i}*K_rec{i};
    e_rec(i, :) = eig(Acl_rec{i});
end

disp(e_rec)

    %% plot the system
if PLOT 
    Box = 1;
    Npts = 10;
    xx = linspace(0, Box, Npts);
    
    % x0 = [1;0];
    figure(1)
    clf
    hold on 
    T = 5;
        for j = 1:Npts
            for k = 1:Npts
                x0 = [ xx(j), xx(k)];
                [tc, xc] = ode23(@(t, x) switch_clp(t, x, Atrue, Btrue, v_rec, K_rec), [0, T], x0);
                plot(xc(:, 1), xc(:, 2))
            end
        end   
    
    xprev = xlim;
    yprev = xlim;
    xlim([0, xprev(2)]);
    ylim([0, yprev(2)]);
    end

else
    disp('Infeasible!')
end



function dx = switch_clp(t, x, A, B, v_rec, K_rec)
    
    Nsys = length(A);
    v = zeros(Nsys, 1);
    for i = 1:Nsys
        v(i) = max(x./v_rec(:, i));
    end
    
    [m, k] = max(v);
    
    dx = A{k}*x + B{k}*K_rec{k}*x;
    

end


