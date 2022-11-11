rng(30, 'twister')

PS = possim(3, 2);


sys = PS.rand_sys(1.4);


T = 20;
traj = PS.sim(T, sys);

poly = data_cons(traj, 1);

pall = reshape([sys.A,sys.B], [], 1);

check = poly.b - poly.A*pall;
