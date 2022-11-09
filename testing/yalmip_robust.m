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

mw=200;
C= randn(mw, L);
d = ones(mw,1);

[C_red, d_red] = nontrivial_constraints(C, d); 
W = C_red*w <= d_red;
F = [uncertain(w); A0*x <= b0+bL*w; W];


opts = sdpsettings;
opts.robust.lplp='duality';
sol = optimize(F, norm(x,'inf'), opts)

value(x)