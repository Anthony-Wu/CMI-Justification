function [Y_p, Y_FIR_i, Y_AR_i] = tf_ARX_pred_inc2 (U_i, A, B, dt, Y)

[nB NA] = cellfun(@size, B);	% nB is a row vector of the number of coefficients for each input
Bmax = max(nB);

[n_sam n_var] = size(U_i);	% number of samples & number of input variable

%for i_var = 1:n_var;
	%iB{1} = B{i_var};	% intermediate term
	%Y_i_cont(:,i_var) = tf_ARX_pred(U_i(:,i_var), A, iB, dt(i_var)); % calculate increment vector
%end
%Y_i = sum(Y_i_cont,2);	% calculate increment at each step

%FIR part
Y_FIR_i = tf_FIR_pred (U_i, B, dt);

% AR part
Y_i = [0; Y(2:end) - Y(1:end-1)];
Y_AR_i = tf_FIR_pred (Y_i, A, 1);

Y_p = zeros(n_sam,1);
Y_p_i = zeros(n_sam,1);

Y_p_i(Bmax+1:end) = Y_p(Bmax:end-1) + Y_FIR_i(Bmax+1:end) - Y_AR_i(Bmax+1:end);

Y_p(1:Bmax) = Y(1:Bmax);
for i_sam = Bmax+1:n_sam;
	Y_p(i_sam) = Y_p(i_sam-1) + Y_p_i(i_sam);
end