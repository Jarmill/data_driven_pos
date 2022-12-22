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
            obj.L = length(traj.L);

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
           n = size(obj.traj.Xdelta, 1);
           m = size(obj.traj.U, 1);
           %decision variables
           y = sdpvar(n, 1);
           S = sdpvar(m, n,obj.L, 'full');
           
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

		
        
        function poly_out = poly_stab(obj, vars)
            %POLY_STAB generate the polytope of positive-stabilizable 
            %subsystems given the values in vars
			%with a common Lyapunov function and controller K (switching-independent)
			
            Nth = size(obj.Th_vert, 2);
            n = length(vars.y);
            for i = 1:Nth
			
                thcurr = obj.Th_vert(:, i);
                Scurr = vars.S(:, :, i);
                %stable system
                Gd_stab = [kron(thcurr', kron(vars.y', eye(n))), kron((Scurr*ones(n, 1))', eye(n))];
                
                %positive system
                Gd_pos = -[kron(thcurr', kron(diag(vars.y), eye(n))), kron(Scurr', eye(n))];
                
            end
%             Gc_pos = Gc_pos_all(M, :);
            
            poly_out = struct;
            poly_out.C = [Gd_stab; Gd_pos];
            poly_out.d = [kron(ones(Nth, 1), vars.y-obj.delta*ones(n,1)); zeros(n^2, 1)];       
        end

    end
end

