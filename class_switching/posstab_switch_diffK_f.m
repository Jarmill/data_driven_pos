classdef posstab_switch_diffK_f < posstab_switch_f
    %posstab_switch_diffk_f  Data-driven stabilization of a switched discrete-time linear 
    % positive system using a different full-state-feedback (gain)  matrix
    % in each subsystem but a common dual-linear copositive control lyapunov function
    %
    %Utilizes the Extended Farkas Lemma, does not require the Yalmip robust
    %optimization toolbox anymore
    
    methods
        function obj = posstab_switch_diffK_f(traj, dopts)
            %POSSTAB_CONT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 2
                dopts = data_opts;
            end
            obj@posstab_switch_f(traj, dopts);					
			
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
		
		function cons_vars = var_cons(obj, vars)
            %constraints on the variables in the program
            cons_vars = [reshape(vars.y, [], 1) >= obj.delta; sum(vars.y, 1)==1]:'Lyap Nonneg';            
        end
				        		
        
        function poly_out = poly_stab(obj, vars)
            %POLY_STAB generate the polytope of positive-stabilizable 
            %subsystems given the values in vars
			%with a common Lyapunov function and controller K (switching-independent)
			
			C_stab_cell = cell(obj.Nsys);
			d_stab = [];
			
			for i = 1:obj.Nsys
				vars_curr = struct('y', vars.y, 'S', squeeze(vars.S(:, :, i)));
				
				poly_out_curr = poly_stab@posstab_f(obj, vars_curr);
				
				C_stab_cell{i} = poly_out_curr.C;
				d_stab = [d_stab; poly_out_curr.d];
			end
			
			C_stab = blkdiag(C_stab_cell{:});
			
			poly_out = struct('C',C_stab, 'd', d_stab);
			                 
        end

    end
end

