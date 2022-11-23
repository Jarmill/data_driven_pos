classdef pos_p2p_cont_f < posstab_cont_f
    %POS_P2P Peak-to-Peak (Linf) gain minimization
    % x+ = A x + B u + E w
    % z  = C x + D u + F w
    %
    % Choose a K: u = K x such that the Linf gain between w and z is
    % minimized.
    
    properties
        param; %parameters (C, D, E, F) as a struct
    end
    
    methods
        function obj = pos_p2p_cont_f(traj, param,  data_options)
            %POS_P2P Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                data_options = data_opts;
            end
            obj@posstab_cont_f(traj, data_options);
            
%             obj.param = param_process(param);
            obj.param = param;
        end
        
        function vars = make_vars(obj)
            vars = make_vars@posstab_cont_f(obj);
            vars.gamma = sdpvar(1, 1);
        end
        
        function poly_out = poly_stab(obj, vars)
            %constraint to enforce stability with p2p performance
            e = size(obj.param.E, 2);
            n = length(vars.y);
%             stab = vars.y' - ones(1, n)*(vars.A*diag(vars.y) + vars.B*vars.S);
%             X = diag(vars.y);

            poly_out = poly_stab@posstab_cont_f(obj, vars);
            
            %stabilization
            poly_out.d(1:n) = poly_out.d(1:n) - obj.param.E*ones(e, 1);

%             term_x = vars.y - (vars.A*X + vars.B*vars.S)*ones(n, 1)...
%                 - obj.param.E*ones(e, 1);
            
%             term_z = gamma*ones(n, 1) - (obj.param.C*X + obj.param.D*vars.S)*ones(n, 1) ...
%                  - obj.param.F*ones(e, 1);
%             cons = ([term_x; term_z] >= obj.delta);
%             stab_label = 'Peak-to-Peak';
        end
        
        function cons_vars = var_cons(obj, vars)
            %constraints on the variables in the program            
            
            %it is no longer required that sum(vars.y) == 1.
            %such a normalization might hurt performance
            [q, n] = size(obj.param.C);
            e = size(obj.param.E, 2);
            
            
            zblk = obj.param.C*diag(vars.y) + obj.param.D*vars.S;
            
            term_z = vars.gamma*ones(q, 1) - zblk*ones(n, 1) - obj.param.F*ones(e, 1);
            
            
%             con_lyap = 
            cons_vars = [(vars.y >= obj.delta):'Lyap Nonneg'; (term_z >= 0):'P2P Output'; ...
                (zblk >= 0):'P2P Mat Nonneg'];            
        end
        
        
        function objective = make_objective(obj, vars)
            objective = vars.gamma;
        end
        
        function out = recover(obj, vars, sol)
            out = recover@posstab_cont_f(obj, vars, sol);
            
            out.gamma = value(vars.gamma);
        end
    end
end

