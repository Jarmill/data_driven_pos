classdef posstab_lpv_cont_f < posstab_lpv_f % & posstab_cont (how to do multiple inheritance in a non-clashing manner?)
    %POSSTAB_SWITCH_CONT Data-driven stabilization of a switched 
    %continuous-time linear positive system using full-state-feedback 
    % and a common linear copositive control lyapunov function
    
    
    methods
        function obj = posstab_lpv_cont_f(traj, Th_vert, dopts)
            %POSSTAB_CONT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                dopts = data_opts;
            end
            obj@posstab_lpv_f(traj, Th_vert, dopts);
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
            L = traj.L;
            
            if data_options.pos_A
            
                %only touch the off-diagonal elements of A
                mask_offdiag = logical(reshape(1-eye(n), [], 1));            
%                 jpos = [find(mask_offdiag); (n^2 + (1:(n*m)))'];

                
                mask_ind= find(mask_offdiag)';
                jpos = reshape(mask_ind'+ (n^2)*(0:(L-1)), 1, []);
                ncon = L*(n^2 - n);
                ipos = 1:ncon;
                vpos = -ones(1, ncon);
                
%                 ipos = (1:(ncon + (n*m)));
%                 vpos = -ones(1, ncon+(n*m));
                
                if data_options.pos_B
                    ipos = [ipos,  (ncon + (1:(n*m)))];
                    jpos = [jpos, ((n^2)*L + (1:n*m))];
                    vpos = [vpos, -ones(1, n*m)];
                    ncon = ncon + n*m;
                end

                Cpos = sparse(ipos, jpos, vpos, ncon, n*(n*L+m));
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
            Nv = size(obj.Th_vert, 2);
            %stable system
            M = metzler_indexer(n);
            Gc_stab = [];
            Gc_pos = [];

            for i = 1:Nv
                thcurr = obj.Th_vert(:, i);
                Scurr = vars.S(:, :, i);
                %stable system
                Gc_stab_curr = [kron(thcurr', kron(vars.y', eye(n))), kron((Scurr*ones(n, 1))', eye(n))];
                Gc_stab = [Gc_stab; Gc_stab_curr];
                
                %positive system
                Gc_pos_curr = -[kron(thcurr', kron(diag(vars.y), eye(n))), kron(Scurr', eye(n))];                
                Gc_pos = [Gc_pos; Gc_pos_curr(M, :)];
            end

            poly_out = struct;
            poly_out.C = [Gc_stab; Gc_pos];
            poly_out.d = [kron(ones(Nv, 1), -obj.delta*ones(n,1)); zeros(size(Gc_pos, 1), 1)];                        
%             poly_out = poly_stab@posstab_lpv_f(obj, vars);
% 
%             poly_out.d(1:(Nv*n)) = poly_out.d(1:(Nv*n)) - kron(ones(Nv, 1), vars.y);
% 			
        end

    end
end

