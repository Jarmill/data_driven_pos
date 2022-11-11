rng(30, 'twister')

n = 3;
m = 2; %feasible

n = 5; %feasible
m = 4;

% n = 6; %infeasible
% m = 4;
PS = possim(n, m);

% sys = PS.rand_sys(1.4);
sys = PS.rand_sys(0.5);

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

stabcon = y' - ones(1, n)*(sys.A*diag(y) + sys.B*S);
poscon =  reshape(sys.A* diag(y) + sys.B*S, [], 1);
cons = [cons; stabcon >= delta; poscon >= 0];



%stabilization constraints


%% solve program
opts = sdpsettings;
opts.robust.lplp='duality';

sol = optimize(cons, norm(y,'inf'), opts);


%% recover
if sol.problem
    disp('Infeasible')
else

    y_rec = value(y);
    v_rec = 1./y_rec;

    S_rec = value(S);
    K_rec = S_rec*diag(v_rec);


    sys_clp = sys.A + sys.B*K_rec
    e_clp = abs(eig(sys_clp))'
end