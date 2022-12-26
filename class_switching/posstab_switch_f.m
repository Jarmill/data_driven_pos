classdef posstab_switch_f < posstab_f
    %POSSTAB_F Data-driven stabilization of a switched discrete-time linear 
    % positive system using full-state-feedback (gain) and a common 
    %dual-linear copositive control lyapunov function
    %
    %Utilizes the Extended Farkas Lemma, does not require the Yalmip robust
    %optimization toolbox anymore
    
    properties
        Nsys; %number of subsystems
    end

    methods
        function obj = posstab_switch_f(traj, dopts)
            %POSSTAB_CONT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 2
                dopts = data_opts;
            end
            traj_split_sys = traj_split(traj);

            obj@posstab_f(traj_split_sys, dopts);
			
			%split up the trajectory based on the active subsystem			
			
			obj.traj = traj_split_sys;
            obj.Nsys = length(obj.traj);

            %report data-consistency polytope information
            %sys_cons:  number of linear constraints describing each
            %           subsystem
            %F_old:     number of linear constraints that were present
            %           before eliminating trivial faces

            obj.poly.sys_cons = obj.poly.F_old.sys_cons;
            obj.poly.F_old = obj.poly.F_old.F_old;

        end
		
        function [vars] = make_vars(obj)
            %generate uncertain and decision variables
           n = size(obj.traj{1}.Xdelta, 1);
           m = size(obj.traj{1}.U, 1);
           %decision variables
           y = sdpvar(n, 1);
           S = sdpvar(m, n, 'full');
           
           vars = struct('y', y, 'S', S);                   
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

            sys_cons = zeros(1, obj.Nsys);
			for i = 1:length(traj)
				[C_cell{i}, d_curr, F_old_curr] = data_cons@posstab_f(obj, traj{i}, data_options);
				d= [d; d_curr];
                sys_cons(i) = length(d_curr);
				F_old = F_old + F_old_curr;				
            end

            F_old = struct('F_old', F_old, 'sys_cons', sys_cons);
			
			%block diagonal matrix
			C= blkdiag(C_cell{:});
            
        end

        %possibly implement per-subsystem sign constraints later.
        %right now, the sign constraints are common among all subsystems
%        function cons_K = controller_cons(obj, vars, dopts)

		
        
        function poly_out = poly_stab(obj, vars)
            %POLY_STAB generate the polytope of positive-stabilizable 
            %subsystems given the values in vars
			%with a common Lyapunov function and controller K (switching-independent)
			
			poly_out_sys = poly_stab@posstab_f(obj, vars);
			
			poly_out = struct('C', kron(eye(obj.Nsys), poly_out_sys.C), 'd', kron(ones(obj.Nsys, 1), poly_out_sys.d));
			                 
        end

        function sys_out = sample_sys(obj, Nsystems)
            %SAMPLE_SYS: randomly sample systems inside the polytope
            %consistency set obj.poly
            %
            

            sys_raw = cprnd(Nsystems, obj.poly.C, obj.poly.d);

            n = obj.traj{1}.n;
            m = obj.traj{1}.m;
            

            sys_out = cell(Nsystems, 1);
            for i = 1:Nsystems
                sys_curr = struct;
                sys_curr.A = cell(obj.Nsys, 1);
                sys_curr.B = cell(obj.Nsys, 1);
                ind = 0;
                for j = 1:obj.Nsys
                    sys_curr.A{j} = reshape(sys_raw(i, ind + (1:n^2)), n, n);
                    ind = ind + n^2;
                    sys_curr.B{j} = reshape(sys_raw(i, ind + (1:(n*m))), n, m);
                    ind = ind + n*m;
                end
                
                

                sys_out{i} = sys_curr;


            end
        end

    end
end

