classdef data_opts
    %DATA_OPTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nontrivial = 1;
        pos_A = 1;
        pos_B = 1;
        
        %sign conditions on the controller
        gez = []; %greater than or equal to zero
        lez = []; %less than or equal to zero
        
        %gez and lez are either empty or are 0/1 arrays of size m by n.
        %when gez(i,j) =lez(i,j)=1, then that element of K is 0.
    end
    

end

