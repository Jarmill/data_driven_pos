rng(30, 'twister')
    % xdot = Acl x + E w
    % z    = C x  F w
    %peak-to-peak w->z
    
    
NOM = 1; %use the nominal controller
    
n = 3;
m = 2; %feasible

sys = struct('A', [-0.562888446377643,0.304206397381707,0.656593915356461;0.062043810844475,-1.344996926420559,0.254522592715944;0.113187181708781,0.132354437133599,0.019971221719390],...
             'B', [0.182453070500558,0.076396428548559;0.472833184316772,0.245344456124525;0.065942821511921,0.940264967809949]);

%a nominal controller that makes A+BK stable
K_nom = [0.086850860753820,-0.475392839969307,6.013859422926600;0.354285451225413,0.798178668811147,-8.034350355031790];         



%% Declare Variables

delta = 1e-3;

%design
y = sdpvar(n, 1);
% S = sdpvar(m, n);
gamma = sdpvar(1, 1);
%parameters of p2p problem

p = 3; %number of external inputs
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

%closed loop 
         
Acl = sys.A + sys.B*K_nom;
Ccl = C + D*K_nom;
%% form constraints

%basic constraints
% cons = [y >= delta; sum(y)==1];
cons = [y>=delta];
X = diag(y);

%condition (6.2) of https://ieeexplore.ieee.org/document/8618689
%equation (12) of https://arxiv.org/pdf/1204.3554.pdf
Abar = Acl*X;
Cbar = Ccl*X;
term_x = -Abar *ones(n, 1) - E*ones(p, 1);
term_z = gamma*ones(q, 1) - Cbar*ones(n, 1) - F*ones(p, 1);

poscon_x =  reshape((1-eye(n)).*Abar, [], 1);
poscon_z = reshape(Cbar, [], 1);
% poscon_z = 0;
cons = [cons; [term_x; term_z] >= delta; [poscon_x; poscon_z] >= 0];


%% solve program
opts = sdpsettings;

sol = optimize(cons,gamma, opts);


%% recover
if sol.problem
    disp('Infeasible')
else

    y_rec = value(y);
    v_rec = 1./y_rec;

    gamma_rec = value(gamma);

end