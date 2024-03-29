classdef posstab_cont_f < posstab_f
    %POSSTAB_CONT Data-driven stabilization of a continuous-time linear 
    %positive system using full-state-feedback and a common linear 
    %copositive control lyapunov function
    
    
    methods
        function obj = posstab_cont_f(traj, dopts)
            %POSSTAB_CONT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 2
                dopts = data_opts;
            end
            obj@posstab_f(traj, dopts);
        end
        
        function [Cpos, dpos] = pos_cons(obj, traj, data_options)
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
            
            if data_options.pos_A
            
                %only touch the off-diagonal elements of A
                mask_offdiag = logical(reshape(1-eye(n), [], 1));            
%                 jpos = [find(mask_offdiag); (n^2 + (1:(n*m)))'];

                ncon = n^2 - n;
                jpos = find(mask_offdiag)';
                ipos = 1:ncon;
                vpos = -ones(1, ncon);
                ncon = sum(mask_offdiag);            
%                 ipos = (1:(ncon + (n*m)));
%                 vpos = -ones(1, ncon+(n*m));
                
                if data_options.pos_B
                    ipos = [ipos,  (ncon + (1:(n*m)))];
                    jpos = [jpos, (n^2 + (1:n*m))];
                    vpos = [vpos, -ones(1, n*m)];
                    ncon = ncon + n*m;
                end

                Cpos = sparse(ipos, jpos, vpos, ncon, n*(n+m));
                dpos = sparse([], [], [], ncon, 1); 
            else
                Cpos = [];
                dpos = [];
            end
            
            
            
        end
        
        function poly_out = poly_stab(obj, vars)
            %POLY_STAB generate the polytope of positive-stabilizable 
            %systems given the values in vars
            
            n = length(vars.y);
            %stable system
            Gc_stab = [kron(vars.y', eye(n)), kron((vars.S*ones(n, 1))', eye(n))];
            
            %positive system
            M = metzler_indexer(n);
            Gc_pos_all = -[kron(diag(vars.y), eye(n)), kron(vars.S', eye(n))];
                
            Gc_pos = Gc_pos_all(M, :);
            
            poly_out = struct;
            poly_out.C = [Gc_stab; Gc_pos];
            poly_out.d = [-obj.delta*ones(n,1); zeros(n^2-n, 1)];                        
        end

    end
end

