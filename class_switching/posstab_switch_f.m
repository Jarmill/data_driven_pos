classdef posstab_switch_f < posstab_f
    %POSSTAB_F Data-driven stabilization of a switched discrete-time linear 
    % positive system using full-state-feedback (gain) and a common 
    %dual-linear copositive control lyapunov function
    %
    %Utilizes the Extended Farkas Lemma, does not require the Yalmip robust
    %optimization toolbox anymore
    
    properties
        Nsys;
    end

    methods
        function obj = posstab_switch_f(traj, dopts)
            %POSSTAB_CONT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 2
                dopts = data_opts;
            end
            obj@posstab_f(traj, dopts);
			
			%split up the trajectory based on the active subsystem
			obj.Nsys = traj.Nys;
			obj.traj = traj_split(traj);
			
        end
		
		function traj = traj_split(obj, traj_in)
			%TRAJ_SPLIT split up the trajectory based on the active subsystem
			%primarily for use in switched systems
			
			
			traj = cell(obj.Nsys, 1);
			
			for i = 1:obj.Nsys
				%isolate data occurring at subsystem i
				traj_curr = struct('n', traj_in.n, 'm', traj_in.m);
				mask_sys = traj_in.sys==i;
				traj_curr.X = traj_in.X(:, mask_sys);
				traj_curr.Xdelta = traj_in.Xdelta(:, mask_sys);
				traj_curr.U = traj_in.U(:, mask_sys);
				
			end
		
		end
        		
		function [C, d, F_old] = data_cons(obj, traj, data_options)
            %DATA_CONS generate the polytope constraint associated with the
            %given trajectory/data
            %Inputs:
            %   traj:   trajectory
            %   data_options:   struct: (nontrivial, pos_A, pos_B)
            %Outputs:
            %   [C, d]: Polytope C x <= d for parameters x 
            %           (vectorized [A;B])
            %   F_old:  Prior number of faces
            
			
			C_cell = cell(obj.Nsys, 1);
			d = [];
			F_old = 0;

			for i = 1:obj.Nsys
				[C_curr, D_curr, F_old_curr] = data_cons@posstab_f(obj, traj{i}, data_options);
				dpos = [dpos; dpos_curr];
				F_old = F_old + F_old_curr;				
			end
			
			%block diagonal matrix
			Cpos = blkdiag(Cpos_cell{:});
            
        end
		
        
        function poly_out = poly_stab(obj, vars)
            %POLY_STAB generate the polytope of positive-stabilizable 
            %subsystems given the values in vars
			%with a common Lyapunov function and controller K (switching-independent)
			
			poly_out_sys = poly_stab@posstab_f(obj, vars);
			
			poly_out = struct('C', kron(poly_out_sys.C, eye(obj.Nsys)), 'd', kron(poly_out_sys.d, ones(obj.Nsys, 1)));
			                 
        end

    end
end

