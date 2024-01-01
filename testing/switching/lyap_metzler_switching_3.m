
A1 = [-1, 2, 0; 2, -3, 1; -1, -1, 3];
A2 = [-4, 2.3, 0.5; 2.3, -1, -0.5; -2.5, -2.5, 2.5];
A3 = [1 2 0; 1 -1 0; 0 1 1];

B1 = [1; 2; -1];
B2 = [0; -0.5; -0.5];
B3 = [0; 0; 1];

A = {A1, A3};
B = {B1, B3};

n = size(A1, 1);
m = size(B1, 2);
N = length(A);

% lam = [-10, 10; 10, -10];

v = sdpvar(n, N, 'full');
S = sdpvar(m, n, N);
Lam = ones(n, N)*10; %a nonnegative matrix satisfying Lam*1 = 1 (markov)

lm_expr = zeros(n, N, 'like', sdpvar);
for i = 1:N
    for j = 1:N
        if j ~= i
            lm_expr(:, i) = lm_expr(:, i) + Lam(i, j)*(v(:, j)- v(:, i));
        end
    end
end

n = length(A1);
z1 = sdpvar(n, 1);
z2 = sdpvar(n, 1);

% z1 = [0.4483; 0.5001];
% z2 = [0.447; 0.5493];


stab_con = [];
pos_con = [];
M = metzler_indexer(n)
for i = 1:N
    stab_con_curr = A{i}*v(:, i) + B{i}*squeeze(S(:, :, i))*ones(n, 1)+ lm_expr(:, i);
    stab_con = [stab_con; stab_con_curr];

    pos_con_curr = reshape(A{i}*diag(v(:, i)) + B{i}*squeeze(S(:, :, i)), [], 1);

    pos_con = [pos_con; pos_con_curr(M)];
end


cons = [stab_con <= (-delta); reshape(v, [], 1) >= delta; pos_con>=0];

sol = optimize(cons, []);

% z = {value(z1), value(z2)};
v_rec = value(v);
S_rec = value(S);
K_rec = cell(N, 1);
Acl_rec = cell(N, 1);

for i = 1:N
    Xi_rec = diag(1./v_rec(:, i));
    K_rec{i} = S_rec(:, :, i)*Xi_rec;
    Acl_rec{i} = A{i}+B{i}*K_rec{i};
end



%% figure out the system

if sol.problem==0
Box = 1;
Npts = 3;
xx = linspace(0, Box, Npts);

% x0 = [1;0];
figure(1)
clf
hold on 
T = 5;
for i = 1:Npts
    for j = 1:Npts
        for k = 1:Npts
            x0 = [xx(i), xx(j), xx(k)];
            [tc, xc] = ode23(@(t, x) switch_clp(t, x, A, B, v_rec, K_rec), [0, T], x0);
            plot3(xc(:, 1), xc(:, 2), xc(:, 3))
        end
    end   
end

xprev = xlim;
yprev = xlim;
zprev = xlim;
xlim([0, xprev(2)]);
ylim([0, yprev(2)]);
zlim([0, zprev(2)]);
view(3);
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


