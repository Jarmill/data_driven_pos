classdef posstab_cont < posstab
    %POSSTAB_CONT Data-driven stabilization of a continuous-time linear 
    %positive system using full-state-feedback and a common linear 
    %copositive control lyapunov function
    
    
    methods
        function obj = posstab_cont(traj)
            %POSSTAB_CONT Construct an instance of this class
            %   Detailed explanation goes here
            obj@posstab(traj);
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
            

            mask_offdiag = logical(reshape(1-eye(n), [], 1));            
%             jpos = (1:n^2);
%             jpos = jpos(mask_offdiag);
            
            jpos = find(mask_offdiag);

            ncon = sum(mask_offdiag);            
            ipos = (1:ncon);
            vpos = -ones(1, ncon);
            
            Cpos = sparse(ipos, jpos, vpos, ncon, n*(n+m));
            dpos = sparse([], [], [], ncon, 1);
            
            
        end
    end
end

