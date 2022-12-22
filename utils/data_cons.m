function [C2, d2] = data_cons(traj, nontrivial)
    %DATA_CONS generate the Linf-polytope constraints associated
    %with the data in traj for a system xdelta = A x + B u
    %
    %Inputs:
    %   traj:   structure with fields
    %               epsilon:    intensity of Linf noise
    %               Xn:         current state
    %               Xdelta:     next state (discrete) or derivative
    %                           (continuous)
    %               U:          applied input
    %
    %Outputs:   
    %   [C, d]:     Polytope C x <= d

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
    
%     for t = 1:T
%         xcurr = traj.Xn(:, t);
%         xnext = traj.Xdelta(:, t);
%         ucurr = traj.U(:, t);
% 
%         Amat = kron(xcurr', eye(n));
%         Bmat = kron(ucurr', eye(n));
%         
%         dcurr = epsilon + kron([1;-1], xnext);
%         blkcurr = [Amat, Bmat; -Amat, -Bmat];
%         
% 
% %         dcurr = epsilon + xnext;
% %         blkcurr = [Amat, Bmat];
%         
%         C(ind + (1:(2*n)), :) = blkcurr;
%         
%         d(ind + (1:(2*n))) = dcurr;
% 
%         ind = ind + (2*n);
%     end
    
    Amat2 = kron(traj.Xn', eye(n));
    Bmat2 = kron(traj.U', eye(n));
    
    C2=  [Amat2, Bmat2; -Amat2, -Bmat2];
    
    d2 = epsilon + kron([1;-1],reshape(traj.Xdelta, [], 1));
    
    
    if nontrivial    
%         [C_red, d_red] = nontrivial_constraints(C, d);
        [C2, d2] = nontrivial_constraints(C2, d2);
        poly = struct('A', C2, 'b', d2, 'F_old', length(d), 'F_new', length(d_red));
    else
        poly = struct('A', C2, 'b', d2);
    end        
    

end