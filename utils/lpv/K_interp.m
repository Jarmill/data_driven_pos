function [K_th] = K_interp(Pth, K, th_in)
% K_INTERP return the self-scheduled controller given vertex controllers K
% evaluated at the parameter value th_in
%
% Input:
%   Pth:    optimizer object returning the map th->c
%   K:      cell array of controllers formed by lpvstab
%   th_in:  current value of the parameters theta
%
% Output:
%   K_th:   gain-scheduled controller formed by a convex combination of
%           controllers K weighted by Pth(th) = c

%input processing
if ~iscell(K)
    K = {K};
end

%solve for weights
c = Pth(th_in);

%form the controller
K_th= zeros(size(K{1}));
for v = 1:length(c)
    K_th = K_th + K{v}*c(v);
end

end

