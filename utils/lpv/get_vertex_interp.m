function [Pth] = get_vertex_interp(Th_vert)
%GET_VERTEX_INTERP return a Yalmip optimizer object that finds a
%representation of the point th (input th) based on a convex combination of
%points in Th_vert (output c)
%
%Input:
%   Th_vert: vertices of the parameter polytope (Th)
%
%Output:
%   Pth:    Yalmip optimizer object that performs th->c


[L, Nv] = size(Th_vert);
th = sdpvar(L, 1);
c = sdpvar(Nv, 1);

%requires solution of a linear program
cons = [c>=0; sum(c)==1; Th_vert*c==th];
Pth = optimizer(cons, 0, sdpsettings('solver', 'linprog'), th, c);

end

