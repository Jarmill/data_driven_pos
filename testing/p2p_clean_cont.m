rng(30, 'twister')
    % xdot = A x + B u + E w
    % z    = C x + D u + F w
    %peak-to-peak w->z
    
    
%P2P: controlled = 3.7416
%   uncontrolled = 32.1781
    
NOM = 1; %use the nominal controller
    
n = 3;
m = 2; %feasible

% sys = struct('A', [-0.562888446377643,0.304206397381707,0.656593915356461;0.062043810844475,-1.344996926420559,0.254522592715944;0.113187181708781,0.132354437133599,0.019971221719390],...
%              'B', [0.182453070500558,0.076396428548559;0.472833184316772,0.245344456124525;0.065942821511921,0.940264967809949]);

sys= struct('A', [-0.2, 0.2, 0.2; 0.4, -0.7, 0.2; 0, 0.8, -3],...
    'B', [-0.4, 0.5; 0.2, -0.8; -1, 2]);

%a nominal controller that makes A+BK stable
% K_nom = [0.086850860753820,-0.475392839969307,6.013859422926600;0.354285451225413,0.798178668811147,-8.034350355031790];         
K_nom = zeros(2, 3);
%this is unstable at the current settings
es = eig(sys.A);


%% Declare Variables

delta = 1e-3;

%design
y = sdpvar(n, 1);
S = sdpvar(m, n);
gamma = sdpvar(1, 1);
%parameters of p2p problem

p = 2; %number of external inputs
% q=n+m; %number of controlled outputs
q = n+m;
C = [eye(n);zeros(m, n)];
D = [zeros(n, m); eye(m)];
% D = 0.1*ones(n, m);
% C = [eye(n); zeros(m, n)]; %like an LQR setup for (C, D)
% C = rand(q, n);
% D = zeros(q, m);
% D = [zeros(n, m); eye(m)*0.5];
E = eye(n, p);
F = zeros(q, p);

%% form constraints

%basic constraints
% cons = [y >= delta; sum(y)==1];
cons = [y>=delta];
X = diag(y);

%condition (6.2) of https://ieeexplore.ieee.org/document/8618689
%equation (12) of https://arxiv.org/pdf/1204.3554.pdf
Abar = (sys.A*X + sys.B*S);
Cbar = (C*X + D*S);
term_x = -Abar *ones(n, 1) - E*ones(p, 1);
term_z = gamma*ones(q, 1) - Cbar*ones(n, 1) - F*ones(p, 1);

poscon_x =  reshape((1-eye(n)).*Abar, [], 1);
poscon_z = reshape(Cbar, [], 1);
% poscon_z = 0;
cons = [cons; [term_x; term_z] >= delta; [poscon_x; poscon_z] >= 0];

if NOM
    nomcon = (S == K_nom * X);
    cons = [cons; nomcon];
end

%% solve program
opts = sdpsettings;

sol = optimize(cons,gamma, opts);


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
    gamma_rec = value(gamma)
    
    out_clean = struct('y', y_rec, 'v', v_rec, 'S', S_rec, 'K', K_rec,...
        'A', sys.A, 'B', sys.B, 'gamma', gamma);
end