classdef posstab_f
    %POSSTAB_F Data-driven stabilization of a discrete-time linear positive 
    %system using full-state-feedback and a common linear copositive 
    %control lyapunov function
    %
    %Utilizes the Extended Farkas Lemma, does not require the Yalmip robust
    %optimization toolbox anymore
    
    properties
        %start traj with only one trajectory, then expand to multiple
        traj; %trajectory trace
        
        delta = 1e-3; %tolerance for a strict inequality to be positive       
        
        opts = sdpsettings('solver', 'mosek','robust.lplp', 'duality');

        poly;


    end
    
    methods
        function obj = posstab_f(traj, data_options)
            %POSSTAB Construct an instance of this class
            %   Detailed explanation goes here
            obj.traj = traj;
             
            if nargin < 2
                data_options = data_opts;
            end
            
            [C, d, F_old] = obj.data_cons(traj, data_options);
            obj.poly = struct('C', C, 'd', d, 'F_old', F_old);
        end
        
        function [out] = solve_program(obj, cons, objective, vars)
            %run the program
            sol = optimize(cons, objective, obj.opts);
            out.sol = sol;
            
            if sol.problem==0
                %successful trajectory execution
                out = obj.recover(vars, sol);
            end
        end   
        
        function [out] = run(obj)
            %STAB: main stabilization routine and execution
            [cons, vars] = obj.make_program();
            objective = obj.make_objective(vars);
            out = obj.solve_program(cons, objective, vars);
        end
        
        %% helper programs
        
        function [cons, vars] = make_program(obj)
           %MAKE_PROGRAM form the LMI program in YALMIP
            vars =  obj.make_vars();
            
            %valid and normalized (inverse) lyapunov weights
            cons_vars = [vars.y >= obj.delta; sum(vars.y)==1];
            
%             pall = reshape([vars.A, vars.B], [], 1);
%             cons_data = (obj.poly.d - obj.poly.C*pall) >= 0;

            poly_s = obj.poly_stab(vars);                       
            
            [cons_contain, Z]= obj.ext_farkas(poly_s, obj.poly);
            vars.Z = Z;
            
%             cons_pos_closed = obj.pos_cons_closed(vars);
            
            cons = [cons_vars:'Lyap Nonneg'; cons_contain];
        end       
        
        function objective = make_objective(obj, vars)
            objective = 0;
        end
        
        
        function [cons_contain, Z]= ext_farkas(obj, poly_s, poly_d)
            %EXT_FARKAS apply the extended farkas lemma
            %Input:
            %   poly_s:     Stabilized polytope     (outer)
            %   poly_d:     Data-consistent polytope(inner)
            %
            %Output:
            %   cons_contain:   constraints
            %   Z:              Nonnegative variable enforcing the
            %                   containment
            
            Z = sdpvar(length(poly_s.d), length(poly_d.d), 'full');
            
            cons_contain = [(Z*poly_d.C == poly_s.C):'Farkas eq';
            (Z*poly_d.d <= poly_s.d):'Farkas ineq';
            (Z>=0):'Farkas nonneg'];
            
            
        end
        
        function poly_out = poly_stab(obj, vars)
            %POLY_STAB generate the polytope of positive-stabilizable 
            %systems given the values in vars
            
            n = length(vars.y);
            %stable system
            Gd_stab = [kron(vars.y', eye(n)), kron((vars.S*ones(n, 1))', eye(n))];
            
            %positive system
            Gd_pos = -[kron(diag(vars.y), eye(n)), kron(vars.S', eye(n))];
                
%             Gc_pos = Gc_pos_all(M, :);
            
            poly_out = struct;
            poly_out.C = [Gd_stab; Gd_pos];
            poly_out.d = [vars.y-obj.delta*ones(n,1); zeros(n^2, 1)];                        
        end
        
        function [vars] = make_vars(obj)
            %generate uncertain and decision variables
           n = size(obj.traj.Xdelta, 1);
           m = size(obj.traj.U, 1);
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
            
            [Cdata, ddata] = data_cons(traj, 0); %data constraints from trajectory
            [Cpos, dpos] = obj.pos_cons(traj, data_options);   %must have been generated by a positive system
            
            F_old = length(ddata) + length(dpos);
            
            if data_options.nontrivial            
                [C, d] = nontrivial_constraints([Cdata; Cpos], [ddata; dpos]);
            else
                C = [Cdata; Cpos];
                d = [ddata; dpos];
            end
            
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
            %           (vectorized [A;B])     
            n = size(traj.Xdelta, 1);
            m = size(traj.U, 1);
            
%             ipos = (1:n^2);
%             jpos = (1:n^2);
%             vpos = -ones(1, n^2);
            
            Ncon = 0;
            if data_options.pos_A
                ipos = (1:n^2);
                jpos = (1:n^2);
                vpos = -ones(1, n^2);
                
                dpos = sparse([], [], [], n^2, 1);
                Ncon = Ncon + n^2;
                
                if data_options.pos_B
                    ipos = [ipos, n^2 + (1:(n*m))];
                    jpos = [jpos, n^2 + (1:(n*m))];
                    vpos = [vpos, -ones(1, n*m)];

                    dpos = [dpos; sparse([], [], [], n*m, 1)];
                    Ncon = Ncon + n*m;
                end
            else
                ipos = [];
                jpos = [];
                vpos = [];
            end
            

            
            Cpos = sparse(ipos, jpos, vpos, Ncon, n*(n+m));
                                  
        end

        
        %% recovery
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
            out.K = out.S*diag(1./out.y);
            
        end
    end
end

