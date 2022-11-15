classdef pos_p2p < posstab
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
        function obj = pos_p2p(traj, param)
            %POS_P2P Construct an instance of this class
            %   Detailed explanation goes here
            obj@posstab(traj);
            
%             obj.param = param_process(param);
            obj.param = param;
        end
        
        function vars = make_vars(obj)
            vars = make_vars@posstab(obj);
            vars.gamma = sdpvar(1, 1);
        end
        
        function [cons, stab_label]= cons_stab(obj, vars)
            %constraint to enforce stability
            n = size(vars.A, 1);
            e = size(obj.param.E, 2);
%             stab = vars.y' - ones(1, n)*(vars.A*diag(vars.y) + vars.B*vars.S);
            X = diag(vars.y);
            term_x = vars.y - (vars.A*X + vars.B*vars.S)*ones(n, 1)...
                - obj.param.E*ones(e, 1);
            
            term_z = gamma*ones(n, 1) - (obj.param.C*X + obj.param.D*vars.S)*ones(n, 1) ...
                 - obj.param.F*ones(e, 1);
            cons = ([term_x; term_z] >= obj.delta);
            stab_label = 'Peak-to-Peak';
        end
        
        function objective = make_objective(obj, vars)
            objective = vars.gamma;
        end
        
        function out = recover(obj, vars, sol)
            out = recover@posstab(obj, vars, sol);
            
            out.gamma = value(gamma);
        end
    end
end

