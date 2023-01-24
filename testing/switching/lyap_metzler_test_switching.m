A1 = [-1, 2; 2, -3];
A2 = [-4, 2.3; 2.3, -1];

A = {A1, A2};

lam = [-10, 10; 10, -10];
n = length(A1);
z1 = sdpvar(n, 1);
z2 = sdpvar(n, 1);

% z1 = [0.4483; 0.5001];
% z2 = [0.447; 0.5493];

delta = 1e-3;
q1 = A1*z1 + lam(2, 1)*(z2 - z1);
q2 = A2*z2 + lam(1, 2)*(z1 - z2);


cons = [[q1; q2] <= (-delta); z1>=delta; z2>=delta];

sol = optimize(cons, []);

z = {value(z1), value(z2)};

%% figure out the system

Box = 1;
Npts = 5;
xx = linspace(-Box, Box, Npts);

% x0 = [1;0];
figure(1)
clf
hold on 
T = 20;
for i = 1:Npts
    for j = 1:Npts
        x0 = [xx(i), xx(j)];
        [tc, xc] = ode23(@(t, x) switch_clp(t, x, A, z), [0, T], x0);
        plot(xc(:, 1), xc(:, 2))
    end
end






function dx = switch_clp(t, x, A, z)
    
    Nsys = length(A);
    v = zeros(Nsys, 1);
    for i = 1:Nsys
        v(i) = max(x./z{i});
    end
    
    [m, k] = max(v);
    
    dx = A{k}*x;
    

end


