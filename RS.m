function [Z, step_table] = RS(n_sam, min_z, max_z, min_l, max_l)
% generates a vertical vector of random steps 
	% number of samples to generate [n_sam]
	% minimum and maximum magnitude of input values [min_z] and [max_z]
	% minimum and maximum step lengths [min_l] and [max_l]
	
Z = zeros(n_sam,1); 
range_z = abs(max_z - min_z);	% step value range
range_l = abs(max_l - min_l);	% step length range

round_lim = 1 - 1e-15;	% rounding limit

bkmk = 0;	% check value
i_step = 1;	% step counter

% step generation
	% calculation would generate a matrix [step_table] of 2 vertical vectors [l, z]
	% i_step is the i-th step, z is the signal value at the i-th step, and l is the step length

n_step = ceil(n_sam/min_l);		% calculate the maximum number of steps (defines size of step_table)
step_table = zeros(n_step,2);	% defines step_table

bkmk = 0;
for i_step = 1:n_step;

	% roll value for l
	step_table(i_step,1) = min_l + floor(rand(1)*(range_l + round_lim));	
	
	% roll value for z
	if i_step == 1;
		step_table(i_step,2) = min_z + floor(rand(1)*(range_z + round_lim));	
	else
		step_table(i_step,2) = min_z + floor(rand(1)*(range_z + round_lim));
		while 	step_table(i_step,2) == step_table(i_step-1,2);	% ensure next step has a new z
				step_table(i_step,2) = min_z + floor(rand(1)*(range_z + round_lim));
		end
	end
	
	% populate vector Z
	if bkmk + step_table(i_step,1) >= n_sam	
		Z(bkmk+1:end) = step_table(i_step,2);
	else
		Z(bkmk+1:bkmk+step_table(i_step,1)) = step_table(i_step,2);	
	end
	bkmk = bkmk+step_table(i_step,1);
end
