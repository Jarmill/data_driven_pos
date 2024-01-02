%discrete-time p2p
%means that the system should be stable already

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

A = {A0, A1, A2};
B = {B0, B1, B2};


n = size(A1, 1);
m = size(B1, 2);
N = length(A);

%p2p properties
C = [eye(n); zeros(2, 3)];
D = [zeros(3, 2); eye(2)];
E = [1 0; 0 0; 0 1];
F = zeros(n, 1);


SCENARIO = 1;
%0: same controller on all modes
%1: different controller on all modes
%2: periodic switching A0 -> A1 -> A2 -> A0


% lam = [-10, 10; 10, -10];

v = sdpvar(n, N, 'full');
S = sdpvar(m, n, N);
gam = sdpvar(1,1);
% Lam = ones(n, N)*50; %a nonnegative matrix satisfying Lam*1 = 1 (markov)



delta = 1e-3;
stab_con = [];
pos_con = [];
% M = metzler_indexer(n);
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


cons = [stab_con <= (-delta); reshape(v, [], 1) >= delta; pos_con>=0; K_con; vcon];

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
    Acl_rec{i} = A{i}+B{i}*K_rec{i};
end


if sol.problem==0
    % abs(eig(Acl_rec{1}))'
    gam_rec
else
    disp('infeasible')
end
