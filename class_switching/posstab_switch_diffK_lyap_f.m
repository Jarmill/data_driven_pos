classdef posstab_switch_diffK_f < posstab_switch_f
    %posstab_switch_diffk_f  Data-driven stabilization of a switched discrete-time linear 
    % positive system using a different full-state-feedback (gain)  matrix
    % in each subsystem and also different dual-linear copositive control lyapunov functions
    %
    %Utilizes the Extended Farkas Lemma, does not require the Yalmip robust
    %optimization toolbox anymore       

	properties
		graph = []; %graph of allowable transitions when switching
	end

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
			d_stab = cell(obj.Nsys, 1);
			
			%add a stability constraint for all allowable arcs (i -> j)					
			%v(i) - A v(j) - B S(j) 1 > 0 
			for i = 1:obj.Nsys
				%origin system i
				d_curr = vars.y(:, i) - ones(n, 1);
				for j = 1:obj.Nsys
					%destination system j
				  	if obj.graph(i,j)2						
						
						C_stab_curr = [kron(vars.y(:, j)', eye(n)), kron((vars.S(:, :, j)*ones(n, 1))', eye(n))];
						
						C_stab_cell{j} = [C_stab_cell{j}; C_stab_curr];
						d_stab{j} = [d_stab{j}; d_curr];
					
				end
			end
			
			C_stab = blkdiag(C_stab_cell{:});
			d_stab = vertcat(d_stab{:});
			
			%add a positivity constraint for all subsystems
			C_pos_cell = cell(obj.Nsys, 1);
			for i = 1:obj.Nsys
				C_pos_cell{i} = -[kron(diag(vars.y(:, i)), eye(n)), kron(vars.S(:, :, i)', eye(n))];
			end
			
			C_pos = blkdiag(C_pos_cell);
			d_pos = zeros(size(C_pos, 1), 1);
			
			
			
			poly_out = struct('C',[C_stab; C_pos], 'd', [d_stab; d_pos]);
			                 
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

