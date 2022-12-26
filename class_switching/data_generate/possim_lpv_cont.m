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
               

         function [out, Kth] = sim_closed_cont(obj, Th_vert, K, sys, x0, T, mu)
            %simulate an continuous-time controlled LPV trajectory 
            %there is no sample noise here.
            
            %Kth(th): gain-scheduled controller at parameter value theta
            out = struct;

            Pth = get_vertex_interp(Th_vert);
            Kth = @(th) K_interp(Pth, K, th); %interpolating controller


            if nargin < 6
                T = 15;
            end
            if nargin < 7
                mu = 0.3;
            end

            if nargin < 4
                sys = obj.rand_sys();
            end
            
            if nargin < 1
                Th_vert = [];
                Kth = @(th) 0;
            end
            
            if nargin < 5
                x0 = zeros(obj.n, 1);
                x0(1) = 1;
            end           

            X = [];
            Th = [];
            U = [];
            Tlog = [];
            
%             X = [x0, zeros(obj.n, T)];
%             U = zeros(obj.m, T);
%             W_true = zeros(obj.n, T);
%             Th = zeros(obj.L, T);
%             %main simulation loop
%             traj = {};
            xprev = x0;
%             trajcurr = struct;
            
            t_all = 0;
            switch_times = 0;
            while t_all < T
                tmax_curr = exprnd(mu);
                tmax_curr = min(t_all + tmax_curr, T) - t_all;
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
                Tlog = [Tlog; t_all + tcurr];
                %storage
%                 X(:, i+1) = xnext;
%                 U(:, i) = ucurr;
%                 W_true(:, i) = wcurr;
%                 Th(:, i) = thcurr;
                xprev = xcurr(end, :)';
                switch_times = [switch_times; t_all + tmax_curr];
                t_all = t_all + tmax_curr;
                
            end
            

            ground_truth = struct;
            ground_truth.A = sys.A;
            ground_truth.B = sys.B;
            
%             struct('W', W_true, 'A', sys.A, 'B', sys.B)

            %package up the output
            out.X = X;
            out.U = U;
            out.Th = Th;    
            out.t = Tlog;
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

