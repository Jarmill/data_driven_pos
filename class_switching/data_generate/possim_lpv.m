classdef possim_lpv
    %possim_lpv Simulation of a LPVA discrete-time positive system    
    
    properties
        n = 3;      %number of states
        m = 2;      %number of inputs
        L = 2;      %number of parameters
                
        epsilon = 0.1;

        sampler = [];

        A_scale = 1.2;      %perturb randomly generated stable ss systems to make them possibly unstable
    end

    methods

        function obj = possim_lpv(n, m, epsilon, L)
            %POSSIM_SWITCH create the data simulator
            %   Detailed explanation goes here

            if nargin >= 3
                obj.n = n;
                obj.m = m;
                obj.L = L;
            end
            
            if nargin >=4
                obj.epsilon = epsilon;
            end
                
            

            obj.sampler = struct('u', @(x, th) 2*rand(obj.m, 1)-1, ...
                                 'w', @() 2*rand(obj.n, 1)-1, ...
                                 'th', @() 2*rand(obj.L,1)-1);
        end

        function out = sim(obj, T, sys, x0)
            %simulate a switched discrete-time positive system trajectory 
            % with Linf bounded noise
            out = struct;

            if nargin < 2
                T = 15;
            end

            if nargin < 3
                sys = obj.rand_sys();
            end
            if nargin < 4
                x0 = ones(obj.n, 1);                
            end           

            
            Xn = [x0, zeros(obj.n, T-1)];      %state
            Xdelta = zeros(obj.n, T);      %state
            U = zeros(obj.m, T);            %input
            W_true = zeros(obj.n, T);       %noise term
            Th = zeros(obj.L, T);           %parameter
            %main simulation loop
            xcurr = x0;
            for i = 1:T
                Xn(:, i) = xcurr;
                %inputs
                
                wcurr = obj.sampler.w()*obj.epsilon;
                thcurr = obj.sampler.th();
                ucurr = obj.sampler.u(xcurr, thcurr);                
                
                %propagation 
                %this is where the noise enters (process ?)
                xnext = sys.B*ucurr + wcurr;
                for k = 1:obj.L
                    xnext = xnext + sys.A{k}*xcurr*thcurr(k);
                end   
                
                %storage
                
                Xdelta(:, i) = xnext;
                U(:, i) = ucurr;
                W_true(:, i) = wcurr;
                Th(:, i) = thcurr;
%                 Th(:, i) = thcurr;
                xcurr = obj.transition_x(xnext);                
            end
            

            ground_truth = struct;
            ground_truth.A = sys.A;
            ground_truth.B = sys.B;
            ground_truth.W = W_true;

            %package up the output
            
            out.Xn = Xn;
            out.Xdelta = Xdelta;
            out.U = U;
            out.epsilon = obj.epsilon;            
            out.ground_truth = ground_truth;
            out.n = obj.n;
            out.m = obj.m;
            out.Th = Th;
            out.L = obj.L;
        end
        


        function xnext = transition_x(obj, xnext_in)
            xnext = xnext_in;
        end

        %% generate sample plants inside the consistency set
        
        function sys_pos = rand_sys(obj, A_scale, bneg)

            %randomly generate the positive system            
            %all entries must be nonnegative
            
            
            
            if nargin < 2
                A_scale = obj.A_scale;
            end
            
            if nargin < 3
                bneg = 0;
            end

            sys_pos = struct;
            sys_pos.A = cell(obj.L, 1);
            sys_pos.B = cell(obj.L, 1);

            for s = 1:obj.L                
                A = abs(randn(obj.n, obj.n));                           
                %package the output                    
                sys_pos.A{s} = A/norm(A)*A_scale;                
            end
            B = randn(obj.n, obj.m);
            if ~bneg
                B = abs(B); %B should be nonnegative if requested
            end
            sys_pos.B = B/norm(B);
        end
        
        function sys_smp = sample_sys(obj, traj, Nsys)
            %sample a set of systems consistent with the data
            %optima of linear objective over the quadratic region
            if nargin <3
                Nsys = 1;
            end
            n = size(traj.X,1);
            [m, T] = size(traj.U);
            
            Xp = traj.X(:, 2:end);
            Xn = traj.X(:, 1:end-1);
            
            %declare variables and error object
            B = sdpvar(n, m, 'full');
            A = cell(L, 1);
            W = Xp -B*traj.U;           
            for k = 1:obj.L
                W = W - A{k}*(traj.Th(k, :) .* Xn);
            end
            
            cons = [];
            for t = 1:T
                cons = [cons; norm(W(:, t), 'inf') <= traj.epsilon]; 
            end
%             W2 = sum(W.^2, 1);

%             cons = (W2' <= traj.epsilon^2);
            
            sys_smp = cell(Nsys, 1);
            %I'm not sure why the yalmip optimizer isn't working here.
%             opts = sdpsettings('solver','mosek', 'verbose', 2);
            opts = sdpsettings('solver','mosek', 'verbose', 0);

            for i = 1:Nsys
%                 sys_curr = struct('A', [], 'B', []);
                CB = randn(n, m);
                
                objective = sum(CB.*B, 'all');
                objective = objective + sum(CA.*A, 'all');
%                 for k = 1:obj.L
%                     CA = randn(n, n);
%                     objective = objective + sum(CA.*A{k}, 'all');
%                 end
                %solve the program                
                sol = optimize(cons, objective, opts);
                
                %recover the solution
                sys_curr = struct;
                sys_curr.B = value(B);
                sys_curr.A = cellfun(@value, A, 'UniformOutput', false);
%                 sys_curr.W2 = value(W2);
                sys_curr.W = value(W);
                sys_curr.Wnorm =max(abs(W), [], 1);
                sys_smp{i} = sys_curr;
            end           
        end
              


    
    end
%     function sa

end

