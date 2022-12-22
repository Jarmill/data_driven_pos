classdef possim_lpv_cont < possim_lpv
    %possim_lpv Simulation of a LPVA continuous-time positive system    
    
    properties
            Rmax = 5; %randomly sample in cube [0, Rmax]^n
        p_pos = 0.5; %diagonal entries of system matrix will have a positive sign with this probability
    end

    methods

        function obj = possim_lpv_cont(n, m, epsilon, L)
            %POSSIM_SWITCH create the data simulator
            %   Detailed explanation goes here

            if nargin < 2
                n = 3;
                m=2 ;
            end
            
            if nargin <3
                epsilon = 0.1;
            end
                
            if nargin <4
                L = 1;
            end


                obj@possim_lpv(n, m, epsilon, L)                      
                                                        
            obj.sampler.x=@() rand(obj.n, 1)*obj.Rmax;
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
                
                ssign = 2*sign(rand(obj.n, 1) - obj.p_pos)-1;   
                %Metzler matrix: off-diagonal signs are allowed to be
                %nonnegative
                for i= 1:obj.n                
                    A(i, i) = ssign(i)*A(i,i);
                end                    

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

         function out = sim_closed_cont(obj, sys, Kth, x0, T, mu)
            %simulate an continuous-time controlled LPV trajectory 
            %there is no sample noise here.
            
            %Kth(th): gain-scheduled controller at parameter value theta
            out = struct;

            if nargin < 5
                T = 15;
            end
            if nargin < 6
                mu = 0.3;
            end

            if nargin < 2
                sys = obj.rand_sys();
            end
            
            if nargin < 3
                Kth = @(th) 0;
            end
            
            if nargin < 4
                x0 = zeros(obj.n, 1);
                x0(1) = 1;
            end           

            X = [];
            Th = [];
            U = [];
            
%             X = [x0, zeros(obj.n, T)];
%             U = zeros(obj.m, T);
%             W_true = zeros(obj.n, T);
%             Th = zeros(obj.L, T);
%             %main simulation loop
%             traj = {};
            xprev = x0;
%             trajcurr = struct;
            tmax_curr = exprnd(mu);          
            t_all = 0;
            switch_times = 0;
            while t_all < T
                tmax_curr = min(t_all + tmax_curr, T);
                %inputs
                
%                 wcurr = obj.sampler.w()*obj.epsilon;
                thcurr = obj.sampler.th();
                Kcurr = Kth(thcurr);
                
                %propagation 
                Ath = 0;
                for ell = 1:length(thcurr)
                    Ath = sys.A{ell}*thcurr(ell);
                end
                
                Acl = Ath + sys.B*Kcurr;
                [tcurr, xcurr] = ode45(@(t, x) Acl*x, [0, tmax_curr], xprev);
                

     
                X = [X, xcurr'];
                U = [U, Kcurr*xcurr'];
                Th = [Th, thcurr];
                %storage
%                 X(:, i+1) = xnext;
%                 U(:, i) = ucurr;
%                 W_true(:, i) = wcurr;
%                 Th(:, i) = thcurr;
                xprev = xcurr(end, :)';
                t_all = t_all + tmax_curr;
                switch_times = [switch_times; tmax_curr];
            end
            

            ground_truth = struct;
            ground_truth.A = sys.A;
            ground_truth.B = sys.B;
            
%             struct('W', W_true, 'A', sys.A, 'B', sys.B)

            %package up the output
            out.X = X;
            out.U = U;
            out.Th = Th;    
            out.ground_truth = ground_truth;
            out.n = obj.n;
            out.m = obj.m;
            out.L = obj.L;
            out.switch_times  = switch_times ;

        end
              

        function xnext = transition_x(obj, xnext_in)
            xnext = obj.sampler.x();
        end


    
    end
%     function sa

end

