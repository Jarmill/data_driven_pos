classdef possim
    %possim Simulation of a discrete-time positive system    
    
    properties
        n = 3;      %number of states
        m = 2;      %number of inputs
        
        epsilon = 0.1;

        sampler = [];

        A_scale = 1.2;      %perturb randomly generated stable ss systems to make them possibly unstable
    end

    methods

        function obj = possim(n, m, epsilon)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here

            if nargin >= 2
                obj.n = n;
                obj.m = m;
            end
            
            if nargin ==3
                obj.epsilon = epsilon;
            end
                

            obj.sampler = struct('u', @(x) 2*rand(obj.m, 1)-1, ...
                                 'w', @() 2*rand(obj.n, 1)-1);
%                                  'w', @() normalize(randn(obj.n,1), 1, 'norm', 2));

            %             obj.Property1 = inputArg1 + inputArg2;
        end

        function out = sim(obj, T, sys, x0)
            %simulate a discrete-time positive system trajectory with Linf 
            %bounded noise
            out = struct;

            if nargin < 2
                T = 15;
            end

            if nargin < 3
                sys = obj.rand_sys();
            end
            if nargin < 4
                x0 = ones(obj.n, 1);
%                 x0(1) = 1;
            end           

            
            X = [x0, zeros(obj.n, T)];
            U = zeros(obj.m, T);
            W_true = zeros(obj.n, T);
            %main simulation loop
            xcurr = x0;
            for i = 1:T
                %inputs
                
                wcurr = obj.sampler.w()*obj.epsilon;
                ucurr = obj.sampler.u(xcurr);
                
                %propagation 
                %this is where the noise enters (process ?)
                xnext = sys.A*xcurr + sys.B*ucurr + wcurr;
%                 for k = 1:obj.L
%                     xnext = xnext + sys.A{k}*xcurr*thcurr(k);
%                 end                
                
                %storage
                X(:, i+1) = xnext;
                U(:, i) = ucurr;
                W_true(:, i) = wcurr;
%                 Th(:, i) = thcurr;
                xcurr = xnext;
            end
            

            ground_truth = struct;
            ground_truth.A = sys.A;
            ground_truth.B = sys.B;
            ground_truth.W = W_true;

            %package up the output
            out.X = X;
            out.Xn = X(:, 1:end-1);
            out.Xdelta = X(:, 2:end);
            out.U = U;
            out.epsilon = obj.epsilon;            
            out.ground_truth = ground_truth;
            out.n = obj.n;
            out.m = obj.m;
%             out.L = obj.L;
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
            
            A = abs(randn(obj.n, obj.n));
            B = randn(obj.n, obj.m);
            if ~bneg
                B = abs(B); %B should be nonnegative if requested
            end
            %package the output

            sys_pos = struct;
            sys_pos.A = A/norm(A)*A_scale;
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
            W = Xp -B*traj.U - A*traj.X;           
            
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

