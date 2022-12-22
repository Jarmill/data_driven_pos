function traj = traj_split(traj_in)
	%TRAJ_SPLIT split up the trajectory based on the active subsystem
	%primarily for use in switched systems
	
	
	traj = cell(traj_in.Nsys, 1);
	
	for i = 1:traj_in.Nsys
		%isolate data occurring at subsystem i
		traj_curr = struct('n', traj_in.n, 'm', traj_in.m, 'epsilon', traj_in.epsilon);
		mask_sys = traj_in.S==i;
		traj_curr.Xn = traj_in.Xn(:, mask_sys);
		traj_curr.Xdelta = traj_in.Xdelta(:, mask_sys);
		traj_curr.U = traj_in.U(:, mask_sys);
		
        traj{i} = traj_curr;
	end

end