classdef pos_p2p_cont < pos_p2p
    %POS_P2P_CONT Peak-to-Peak (Linf) gain minimization
    % xdot = A x + B u + E w
    % z  = C x + D u + F w
    %
    % Choose a K: u = K x such that the Linf gain between w and z is
    % minimized.
    
    methods
        function obj = pos_p2p_cont(traj, param)
            %POS_P2P_CONT Construct an instance of this class
            %   Detailed explanation goes here
            obj@pos_p2p(traj, param);
        end
        
    function [Cpos, dpos] = pos_cons(obj, traj)
            %POS_CONS generate linear constraints such that the
            %ground-truth system is positive
            %
            %for continuous-time, A is Metzler.
            %all off-diagonal elements of A are positive
            %Inputs:
            %   traj:   trajectory
            %
            %Outputs:
            %   [Cpos, dpos]: Polytope Cpos x <= dpos for parameters x 
            %           (vectorized [A;B])     
            
            
            n = size(traj.Xdelta, 1);
            m = size(traj.U, 1);
            
            %only touch the off-diagonal elements of A
            mask_offdiag = logical(reshape(1-eye(n), [], 1));            
            jpos = [find(mask_offdiag); (n^2 + (1:(n*m)))'];
            
            ncon = sum(mask_offdiag);            
            ipos = (1:(ncon + (n*m)));
            vpos = -ones(1, ncon+(n*m));
            
            Cpos = sparse(ipos, jpos, vpos, ncon+(n*m), n*(n+m));
            dpos = sparse([], [], [], ncon+(n*m), 1);     
        end
        
        function [cons, stab_label]= cons_stab(obj, vars)
            %constraint to enforce stability
            n = size(vars.A, 1);
            e = size(obj.param.E, 2);
            p = size(obj.param.C, 1);
%             stab = vars.y' - ones(1, n)*(vars.A*diag(vars.y) + vars.B*vars.S);
            X = diag(vars.y);
            term_x = - (vars.A*X + vars.B*vars.S)*ones(n, 1)...
                - obj.param.E*ones(e, 1);
            
            term_z = vars.gamma*ones(p, 1) - (obj.param.C*X + obj.param.D*vars.S)*ones(n, 1) ...
                 - obj.param.F*ones(e, 1);
            cons = ([term_x; term_z] >= obj.delta);
            cons = [cons; vars.gamma >= 0];
            stab_label = 'Peak-to-Peak';
        end
        
        
        function cons = pos_cons_closed(obj, vars)
            %the closed-loop system should also be positive (Metzler)
            n = size(vars.A, 1);
            sblock = (vars.A*diag(vars.y) + vars.B*vars.S);
            pall = reshape((1-eye(n)).*sblock, [], 1);
            cons = (pall >= 0);
            
%             cons = [];

%             screen = 1- 2*eye(n);

        end
    end
end

