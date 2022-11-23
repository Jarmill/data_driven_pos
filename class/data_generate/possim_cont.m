classdef possim_cont < possim
    %possim_cont Simulation of a continuous-time positive system    
    
    properties        
        Rmax = 5; %randomly sample in cube [0, Rmax]^n
        p_pos = 0.5; %diagonal entries of system matrix will have a positive sign with this probability
    end

    methods

        function obj = possim_cont(n, m, epsilon)
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
            
            if nargin < 3
                epsilon = 0.1;
            end

                obj@possim(n, m, epsilon)
%             end
            
            
                

            obj.sampler.x=@() rand(obj.n, 1)*obj.Rmax;
%                                  'w', @() normalize(randn(obj.n,1), 1, 'norm', 2));

            %             obj.Property1 = inputArg1 + inputArg2;
        end

        function out = sample_slope(obj, T, sys)
            %randomly sample points in [0, Rmax] to obtain data
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
            %main simulation loop
%             xcurr = x0;
            for i = 1:T
                %inputs
                
                wcurr = obj.sampler.w()*obj.epsilon;
                ucurr = obj.sampler.u();
                xcurr = obj.sampler.x();
                
                %propagation 
                %this is where the noise enters (process ?)
                xdelta = sys.A*xcurr + sys.B*ucurr + wcurr;
%                 for k = 1:obj.L
%                     xnext = xnext + sys.A{k}*xcurr*thcurr(k);
%                 end                
                
                %storage
                Xn(:, i) = xcurr;
                Xdelta(:, i) = xdelta;
                U(:, i) = ucurr;
                W_true(:, i) = wcurr;
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
%             out.L = obj.L;
        end
        
        
        function xout = sim_cont(obj)
            %SIM_CONT simulate a trajectory in continuous-time
            xout = 0;
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
            
            
            A = abs(randn(obj.n, obj.n));
            
            B = randn(obj.n, obj.m);
            if ~bneg
                B = abs(B); %B should be nonnegative
            end
            
            
            s = 2*sign(rand(obj.n, 1) - obj.p_pos)-1;
            %Metzler matrix: off-diagonal signs are allowed to be
            %nonnegative
            for i= 1:obj.n                
                A(i, i) = s(i)*A(i,i);
            end
            
            %package the output

            sys_pos = struct;
            sys_pos.A = A/norm(A)*A_scale;
            sys_pos.B = B/norm(B);
        end               

    
    end

end

