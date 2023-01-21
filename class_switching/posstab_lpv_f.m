classdef posstab_lpv_f < posstab_f
    %POSSTAB_LPV_F Data-driven stabilization of a switched discrete-time 
    % linear parameter-varying positive system using full-state-feedback 
    % (gain) and a common dual-linear copositive control lyapunov function
    %
    %Utilizes the Extended Farkas Lemma, does not require the Yalmip robust
    %optimization toolbox anymore
    
    properties        
        Th_vert; 
        L;
    end

    methods
        function obj = posstab_lpv_f(traj, Th_vert, dopts)
            %POSSTAB_CONT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                dopts = data_opts;
            end

            
            

            obj@posstab_f(traj, dopts);
			
			%load in the parameter data
			obj.Th_vert  = Th_vert;
            obj.L = traj.L;

        end
		
        function [vars] = make_vars(obj)
            %generate uncertain and decision variables
           n = size(obj.traj.Xdelta, 1);
           m = size(obj.traj.U, 1);
           %decision variables
           y = sdpvar(n, 1);
           S = sdpvar(m, n,size(obj.Th_vert, 2), 'full');
           
           vars = struct('y', y, 'S', S);                   
        end


        function [C2, d2] = data_cons_poly(obj, traj)
            %DATA_CONS generate the Linf-polytope constraints associated
            %with the data in traj for a system xdelta = A x + B u
            %
            %Inputs:
            %   traj:   structure with fields
            %               epsilon:    intensity of Linf noise
            %               Xn:         current state
            %               Xdelta:     next state (discrete) or derivative
            %                           (continuous)
            %               U:          applied input
            %
            %Outputs:   
            %   [C, d]:     Polytope C x <= d
        
            %           (vectorized [A1, A2, A3... AL, B])     
            %a modification of utils/data_cons.m
            
            %unpack the important inputs
            [n, T] = size(traj.Xdelta);
            m = size(traj.U, 1);
            L = traj.L;
            epsilon = traj.epsilon;
 
            %the column operations to get the x data description
            Amat2 = [];
            for i = 1:traj.L
                Amat2 = [Amat2, kron((traj.Xn.*traj.Th(i, :))', eye(n))];
            end
            
            Bmat2 = kron(traj.U', eye(n));
            
            C2=  [Amat2, Bmat2; -Amat2, -Bmat2];
            
            d2 = epsilon + kron([1;-1],reshape(traj.Xdelta, [], 1));
        
        end

        function [Cpos, dpos] = pos_cons(obj, traj, data_options)
            %POS_CONS generate linear constraints such that the
            %ground-truth system is positive
            %
            %for discrete-time, all elements of A must be nonnegative
            %also, all elements of B must be nonnegative
            %Inputs:
            %   traj:   trajectory
            %
            %Outputs:
            %   [Cpos, dpos]: Polytope Cpos x <= dpos for parameters x 
            %           (vectorized [A1, A2, A3... AL, B])     
            n = size(traj.Xdelta, 1);
            L = traj.L;
            m = size(traj.U, 1);
            Ncon = 0;            
            if data_options.pos_A
                ipos = (1:(L*n^2));
                jpos = (1:(L*n^2));
                vpos = -ones(1, (L*n^2));
                
                dpos = sparse([], [], [], (L*n^2), 1);
                Ncon = Ncon + (L*n^2);
                
                if data_options.pos_B
                    ipos = [ipos, (L*n^2) + (1:(n*m))];
                    jpos = [jpos, (L*n^2) + (1:(n*m))];
                    vpos = [vpos, -ones(1, n*m)];

                    dpos = [dpos; sparse([], [], [], n*m, 1)];
                    Ncon = Ncon + n*m;
                end
            else
                ipos = [];
                jpos = [];
                vpos = [];
                dpos = [];
            end
            
            Cpos = sparse(ipos, jpos, vpos, Ncon, n*(L*n+m));
        end

        %possibly implement per-subsystem sign constraints later.
        %right now, the sign constraints are common among all subsystems
%        function cons_K = controller_cons(obj, vars, dopts)

		
        function sys_out = sample_sys(obj, Nsystems)
            %SAMPLE_SYS: randomly sample systems inside the polytope
            %consistency set obj.poly

            opts_cprnd = struct('discard', 100);
            sys_raw = cprnd(Nsystems, obj.poly.C, obj.poly.d, opts_cprnd);

            n = obj.traj.n;
            m = obj.traj.m;
            L = obj.traj.L;

            sys_out = cell(Nsystems, 1);
            for i = 1:Nsystems
                sys_curr = struct;
                sys_curr.A = cell(L, 1);
                ind = 0;
                for j = 1:L
                    sys_curr.A{j} = reshape(sys_raw(i, ind + (1:n^2)), n, n);
                    ind = ind + n^2;
                end
                
                sys_curr.B = reshape(sys_raw(i, ind + (1:(n*m))), n, m);

                sys_out{i} = sys_curr;


            end
        end

        
        function poly_out = poly_stab(obj, vars)
            %POLY_STAB generate the polytope of positive-stabilizable 
            %subsystems given the values in vars
			%with a common Lyapunov function and controller K (switching-independent)
			
            Nv = size(obj.Th_vert, 2);
            n = length(vars.y);
            Gd_stab = [];
            Gd_pos = [];
            for i = 1:Nv
			
                thcurr = obj.Th_vert(:, i);
                Scurr = vars.S(:, :, i);
                %stable system
                Gd_stab_curr = [kron(thcurr', kron(vars.y', eye(n))), kron((Scurr*ones(n, 1))', eye(n))];
                Gd_stab = [Gd_stab; Gd_stab_curr];
                
                %positive system
                Gd_pos_curr = -[kron(thcurr', kron(diag(vars.y), eye(n))), kron(Scurr', eye(n))];
                Gd_pos = [Gd_pos; Gd_pos_curr];
            end
%             Gc_pos = Gc_pos_all(M, :);
            
            poly_out = struct;
            poly_out.C = [Gd_stab; Gd_pos];
            poly_out.d = [kron(ones(Nv, 1), vars.y-obj.delta*ones(n,1)); zeros((Nv*n^2), 1)];       
        end

        function out = recover_controller(obj, vars)
            %RECOVER get the controllers and parameters                       

            Nv = size(obj.Th_vert, 2);
            for i = 1:Nv
                out.K{i} = out.S(:, :, i)*diag(1./out.y);
            end
        end

    end
end

