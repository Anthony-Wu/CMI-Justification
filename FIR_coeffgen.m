function [B] = FIR_coeffgen (n_coeff, K, tau, dt)
% generates FIR coefficients for a first-order system with [n_coeff] number of non-zero coefficients, [K] gain, [tau] time constant and [dt] dead time

n_var = size(n_coeff,2); 	% define number of input parameters

for i_var = 1:n_var
	B_int{i_var} = lsim(tf(K(i_var),[tau(i_var) 1]),[zeros(dt(i_var),1);1;zeros(n_coeff(i_var),1)],0:1:dt(i_var)+n_coeff(i_var));
	B{i_var} = B_int{i_var}(2:end);
end