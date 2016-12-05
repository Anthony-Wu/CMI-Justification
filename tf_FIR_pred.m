function [Y, Y_cont] = tf_FIR_pred (U, B, dt)
% estimates the output of a system using an FIR model structure

[n_sam n_var] = size(U);

Y_cont_int = zeros(n_sam+1, n_var);
Y_cont = zeros(n_sam, n_var);
Y = zeros(n_sam, 1);

Ts = 1;
for i_var = 1:n_var;
	G(i_var) = tf( B{i_var}', 1, Ts,'InputDelay', dt(i_var), 'variable', 'q^-1');
	Y_cont_int(:,i_var) = lsim(G(i_var),[0;U(:,i_var)],[0:1:n_sam]);
	Y_cont(:,i_var) = Y_cont_int(2:end,i_var);
end

Y = sum(Y_cont,2);
