classdef posstab_switch_diffK_cont_f < posstab_switch_cont_f
    %posstab_switch_diffk_cont _f  Data-driven stabilization of a switched 
    % continuous-time linear positive system using a different full-state-
    % feedback (gain)  matrix in each subsystem but a common dual-linear 
    % copositive control lyapunov function
    %
    %Utilizes the Extended Farkas Lemma, does not require the Yalmip robust
    %optimization toolbox anymore
    
    methods
        function obj = posstab_switch_diffK_cont_f(traj, dopts)
            %POSSTAB_CONT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 2
                dopts = data_opts;
            end
            obj@posstab_switch_cont_f(traj, dopts);					
			
        end
		
		
		function [vars] = make_vars(obj)
            %generate uncertain and decision variables
           n = size(obj.traj{1}.Xdelta, 1);
           m = size(obj.traj{1}.U, 1);
		   Nsys = obj.Nsys;
           %decision variables
           y = sdpvar(n, 1, 'full');
           S = sdpvar(m, n, Nsys, 'full');
           
           vars = struct('y', y, 'S', S);                   
        end
		
				        		
        
        function poly_out = poly_stab(obj, vars)
            %POLY_STAB generate the polytope of positive-stabilizable 
            %subsystems given the values in vars
			%with a common Lyapunov function and controller K (switching-independent)
			
            n = length(vars.y);
			C_stab_cell = cell(obj.Nsys);
			d_stab = kron(ones(obj.Nsys, 1), [-obj.delta*ones(n,1); zeros(n^2-n, 1)]);
			

            M = metzler_indexer(n);
			for i = 1:obj.Nsys
				%stabilization
                Gc_stab = [kron(vars.y', eye(n)), kron((vars.S(:, :, i)*ones(n, 1))', eye(n))];
                
                Gc_pos_all = -[kron(diag(vars.y), eye(n)), kron(vars.S(:, :, i)', eye(n))];                
                Gc_pos = Gc_pos_all(M, :);
				
				C_stab_cell{i} = [Gc_stab; Gc_pos];
				
			end
			
			C_stab = blkdiag(C_stab_cell{:});
			
			poly_out = struct('C',C_stab, 'd', d_stab);
			                 
        end

    end
end

