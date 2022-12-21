classdef posstab_switch_diffK_f < posstab_switch_f
    %posstab_switch_diffk_f  Data-driven stabilization of a switched discrete-time linear 
    % positive system using a different full-state-feedback (gain)  matrix
    % in each subsystem and also different dual-linear copositive control lyapunov functions
    %
    %Utilizes the Extended Farkas Lemma, does not require the Yalmip robust
    %optimization toolbox anymore       

    methods
        function obj = posstab_switch_diffK_f(traj, dopts)
            %POSSTAB_CONT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                dopts = data_opts;
            end
            obj@posstab_switch_f(traj, dopts);	

				
			if nargin < 2
				obj.graph = ones(obj.Nsys);
			end					
			
        end
		
		
		function [vars] = make_vars(obj)
            %generate uncertain and decision variables
           n = size(obj.traj.Xdelta, 1);
           m = size(obj.traj.U, 1);
		   Nsys = obj.Nsys;
           %decision variables
           y = sdpvar(n, Nsys, 'full');
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
			%with a common Lyapunov function and different controllers K per subsystem (switching-dependent)
			
			C_stab_cell = cell(obj.Nsys, 1);
			d_stab = [];
			
			%add a stability constraint for all allowable arcs (i -> j)					
			for i = 1:obj.Nsys
				for j = 1:obj.Nsys
				  	if obj.graph(i,j)
						vars_curr = struct('y', vars.y, 'S', vars.S(:, :, i);
						
						poly_out_curr = posstab_f@poly_stab(vars);
						
						C_stab_cell{i} = poly_out_curr.C;
						d_stab = [d_stab; poly_out_curr.d];
					
				end
			end
			
			%add a positivity constraint for all subsystems
			C_pos = cell(obj.Nsys, 1
			
			C_stab = blkdiag(C_stab_cell{:});
			
			poly_out = struct('C',C_stab, 'd', d_stab);
			                 
        end
		
		function out = recover(obj, vars, sol)
            %RECOVER get the controllers and parameters
            out = struct;
            %variables
            out.sol = sol;
            
            %copositive linear control lyapunov function
            out.y = value(vars.y);
%             out.v = 1./out.y;
            out.Z = value(vars.Z);
            
            %control action
            out.S = value(vars.S);
			out.K = cell(obj.Nsys, 1);
			for i = 1:obj.Nys
				out.K{i} = out.S(:, :, i)*diag(1./out.y(:, i));
			end
            
            
        end

    end
end

