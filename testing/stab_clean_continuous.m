rng(30, 'twister')

n = 3;
m = 2; %feasible

% n = 5; %feasible
% m = 4;

% n = 7; %feasible
% m = 7;

% n = 4; %feasible
% m = 5;
PS = possim_cont(n, m);

sys = PS.rand_sys(1.4);
% sys = PS.rand_sys(0.5);

%this is unstable at the current settings
es = eig(sys.A);


%% Declare Variables

delta = 1e-3;

%design
y = sdpvar(n, 1);
S = sdpvar(m, n);

%uncertain
% A = sdpvar(n, n, 'full');
% B = sdpvar(n, m, 'full');

%% form constraints

%basic constraints
cons = [y >= delta; sum(y)==1];

% stabcon = -ones(1, n)*(sys.A*diag(y) + sys.B*S);
stabcon = -(sys.A*diag(y) + sys.B*S)*ones(n, 1);
poscon =  reshape((1-eye(n)).*(sys.A* diag(y) + sys.B*S), [], 1);
cons = [cons; stabcon >= delta; poscon >= 0];


%V = v'x
%Vdot = v'xdot = v'(Ax+BKx) = v'(A+BK)x <0
%Vdot = x'(A' + K' B')v
%Vdot = (A' + K' B')v < 0
%Vdot = diag(y) (A' + K' B') v diag(y) < 0
%Vdot = diag(y) (A' + K' B') 1 < 0

%stabilization constraints


%% solve program
opts = sdpsettings;
opts.robust.lplp='duality';

sol = optimize(cons,0, opts);


%% recover
if sol.problem
    disp('Infeasible')
else

    y_rec = value(y);
    v_rec = 1./y_rec;

    S_rec = value(S);
    K_rec = S_rec*diag(v_rec);


    sys_clp = sys.A + sys.B*K_rec
    e_clp = eig(sys_clp)'
    
    out_clean = struct('y', y_rec, 'v', v_rec, 'S', S_rec, 'K', K_rec,...
        'A', sys.A, 'B', sys.B);
end