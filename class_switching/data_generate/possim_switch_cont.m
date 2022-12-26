classdef possim_switch_cont < possim_switch
    %possim_switch_cont Simulation of a switched continuous-time positive system    
    
    properties        
        Rmax = 5; %randomly sample in cube [0, Rmax]^n
        p_pos = 0.5; %diagonal entries of system matrix will have a positive sign with this probability
    end

    methods

        function obj = possim_switch_cont(n, m, epsilon, Nsys)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here

%             if nargin ==0
%                 obj@possim()
%             end
            
%             if nargin ==3

            if nargin < 2
                n = 3;
                m=2 ;
            end
            
            if nargin <3
                epsilon = 0.1;
            end
                
            if nargin <4
                Nsys = 1;
            end


                obj@possim_switch(n, m, epsilon, Nsys)                      
                

            obj.sampler.x=@() rand(obj.n, 1)*obj.Rmax;

            obj.Nsys = Nsys;
        end

        function out = sample_slope(obj, T, sys)
            %randomly sample x points in [0, Rmax] to obtain data
            out = struct;

            if nargin < 2
                T = 15;
            end

            if nargin < 3
                sys = obj.rand_sys();
            end        

            
            Xn = zeros(obj.n, T);
            Xdelta = zeros(obj.n, T);
            U = zeros(obj.m, T);
            W_true = zeros(obj.n, T);
            S = zeros(1, T);
            %main simulation loop
%             xcurr = x0;
            for i = 1:T
                %inputs
                
                wcurr = obj.sampler.w()*obj.epsilon;
                xcurr = obj.sampler.x();
                scurr = obj.sampler.sys(xcurr);
                ucurr = obj.sampler.u(xcurr, scurr);
                
                
                %propagation 
                %this is where the noise enters (process ?)
                xdelta = sys.A{scurr}*xcurr + sys.B{scurr}*ucurr + wcurr;
%                 for k = 1:obj.L
%                     xnext = xnext + sys.A{k}*xcurr*thcurr(k);
%                 end                
                
                %storage
                Xn(:, i) = xcurr;
                Xdelta(:, i) = xdelta;
                U(:, i) = ucurr;
                W_true(:, i) = wcurr;
                S(1, i) = scurr;
%                 Th(:, i) = thcurr;
%                 xcurr = xnext;
            end
            

            ground_truth = struct;
            ground_truth.A = sys.A;
            ground_truth.B = sys.B;            
            ground_truth.W = W_true;

            %package up the output
%             out.X = X;
            out.Xn = Xn;
            out.Xdelta = Xdelta;
            out.U = U;
            out.epsilon = obj.epsilon;            
            out.ground_truth = ground_truth;
            out.n = obj.n;
            out.m = obj.m;
            out.Nsys = obj.Nsys;
            out.S = S;
%             out.L = obj.L;
        end
        
        
        function traj = sim_closed_cont(obj, sys, K, Tsim, x0, mu0)
            %SIM_CONT simulate a trajectory in continuous-time

            if nargin < 2
                sys = obj.rand_sys();
            end

            if nargin < 3
                Tsim = 15;
            end

            if nargin < 4
                x0 = ones(obj.n, 1);
            end

            if nargin < 6
                mu0 = 0.3;
            end

            if ~iscell(K)
                Nsys = length(sys.A);
                K0 = K;
                K = cell(Nsys, 1);
                for i = 1:Nsys
                    K{i} = K0;
                end
            end


            T = [];
            X = [];
            U = [];
            S = [];

            xprev = x0;
            
            t_all = 0;
            switch_times = 0;
            while t_all < Tsim
                tmax_curr = exprnd(mu0);          
                tmax_curr = min(t_all + tmax_curr, Tsim) - t_all;
                scurr = obj.sampler.sys(xprev);

                Kcurr = K{scurr};

                Acl = sys.A{scurr} + sys.B{scurr}*Kcurr;

                [tcurr, xcurr] = ode45(@(t, x) Acl*x, [0, tmax_curr], xprev);                    
                X = [X, xcurr'];
                U = [U, Kcurr*xcurr'];
                
                xprev = xcurr(end, :)';
                
                T = [T, t_all + tcurr'];
                switch_times = [switch_times; t_all + tmax_curr];
                t_all = t_all + tmax_curr;
                S = [S; scurr];
            end
            
            traj = struct('t', T, 'X', X, 'U', U, 'switch_times', switch_times, 'S', S);

        end

        %% generate sample plants inside the consistency set
        
        function sys_pos = rand_sys(obj, A_scale, bneg)

            %randomly generate the positive system            
            %must be a Metzler matrix                        
            
            if nargin < 2
                A_scale = obj.A_scale;
            end
            
            
            if nargin < 3
                bneg = 0;
            end
            
            
            sys_pos = struct;
            sys_pos.A = cell(obj.Nsys, 1);
            sys_pos.B = cell(obj.Nsys, 1);

            for s = 1:obj.Nsys
                
                %State matrix A
                A = abs(randn(obj.n, obj.n));

                ssign = 2*sign(rand(obj.n, 1) - obj.p_pos)-1;   
                %Metzler matrix: off-diagonal signs are allowed to be
                %nonnegative
                for i= 1:obj.n                
                    A(i, i) = ssign(i)*A(i,i);
                end
                
                %input matrix B
                B = randn(obj.n, obj.m);
                if ~bneg
                    B = abs(B); %B should be nonnegative if requested
                end
                %package the output
    
                
                sys_pos.A{s} = A/norm(A)*A_scale;
                sys_pos.B{s} = B/norm(B);
            end
    
        end               

    
    end

end

