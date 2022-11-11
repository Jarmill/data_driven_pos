function poly = data_cons(traj, nontrivial)
    %DATA_CONS generate the Linf-polytope constraints associated
    %with the data in traj
    %
    %Inputs:
    %   traj:   structure with fields
    %               epsilon:    intensity of Linf noise
    %               Xn:         current state
    %               Xdelta:     next state (discrete) or derivative
    %                           (continuous)
    %               U:          applied input

    %column-wise vectorization of [A; B];

    if nargin < 2
        nontrivial = 1;
    end
    
    %unpack the important inputs
    [n, T] = size(traj.Xdelta);
    m = size(traj.U, 1);
    epsilon = traj.epsilon;

    C = zeros(2*(n*T), n*(n+m));
    d = zeros(2*(n*T), 1);
    ind = 0;
    
    if nargin >2 && test
        A = sdpvar(n, n, 'full');
        B = sdpvar(n, m, 'full');
        pa = reshape([A,B], [], 1);
    end
    
    for t = 1:T
        xcurr = traj.Xn(:, t);
        xnext = traj.Xdelta(:, t);
        ucurr = traj.U(:, t);

        Amat = kron(xcurr', eye(n));
        Bmat = kron(ucurr', eye(n));
        
        dcurr = epsilon + kron([1;-1], xnext);
        blkcurr = [Amat, Bmat; -Amat, -Bmat];
        

%         dcurr = epsilon + xnext;
%         blkcurr = [Amat, Bmat];
        
        C(ind + (1:(2*n)), :) = blkcurr;
        
        d(ind + (1:(2*n))) = dcurr;

        ind = ind + (2*n);
    end
    
    if nontrivial    
        [C_red, d_red] = nontrivial_constraints(C, d);
        poly = struct('A', C_red, 'b', d_red);
    else
        poly = struct('A', C, 'b', d);
    end
    

end