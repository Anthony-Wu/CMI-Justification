function [con_A con_b] = CMI_con(nZ, K_sign, min_phase)
% sets up the inequality constraint matrix for constrained model identification
% nZ = number of non-zero coefficients
% K_sign = direction of static gain (1 +ve, 0 do not set, -1 -ve)
% min_phase = minimum phase? (0 = no, 1 = yes)
% assumes sum(first row) == number of parameters
	n_paras = sum(nZ);
	n_var = size(nZ,2);	% total number of variables 

	
	% Define con_A
	con_A = zeros(1,n_paras);	% define row length of con_A
	bkmk_c = 0;	% condition bookmark
	bkmk_z = 0;	% variable bookmark
	
	for i_var = 1:n_var;
		if K_sign(i_var) == 1;
			if 	min_phase(i_var) == 0;	% not minimum phase system (only 1 condition)
				con_A(bkmk_c+1,bkmk_z+1:bkmk_z+nZ(i_var)) = -ones(1,nZ(i_var));
			elseif min_phase(i_var) == 1;
				con_A(bkmk_c+1:bkmk_c+nZ(i_var),bkmk_z+1:bkmk_z+nZ(i_var)) = -diag(ones(1,nZ(i_var)));
			end
		end
		
		if K_sign(i_var) == -1;
			if 	min_phase(i_var) == 0;	% not minimum phase system (only 1 condition)
				con_A(bkmk_c+1,bkmk_z+1:bkmk_z+nZ(i_var)) = ones(1,nZ(i_var));
			elseif min_phase(i_var) == 1;
				con_A(bkmk_c+1:bkmk_c+nZ(i_var),bkmk_z+1:bkmk_z+nZ(i_var)) = diag(ones(1,nZ(i_var)));
			end
		end
		bkmk_z = bkmk_z + nZ(i_var);
		bkmk_c = size(con_A,1);
	end
	
	con_b = zeros(bkmk_c,1);