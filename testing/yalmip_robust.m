rng(4,'twister');
n=10;
m=15;
L = 8;

A0=randn(m,n);
x0= normalize(randn(n,1));

b0= 2*A0*x0;


bL = randn(m, L);

x=sdpvar(n,1);
w=sdpvar(L,1);


C= randn(20, L);

W = C*w <= 0.1;
F = [uncertain(w); A0*x <= b0+bL*w; W];


opts = sdpsettings;
opts.robust.lplp='duality';
sol = optimize(F, [], opts)

value(x)