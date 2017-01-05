function [Y_p, Y_i, Y_i_cont] = tf_ARX_pred_inc (U_i, A, B, dt, Y)

[nB NA] = cellfun(@size, B);	% nB is a row vector of the number of coefficients for each input
Bmax = max(nB);

[n_sam n_var] = size(U_i);	% number of samples & number of input variable

Y_p = zeros(n_sam,1);
Y_i = zeros(n_sam,1);		
Y_i_cont = zeros(n_sam, n_var);	

for i_var = 1:n_var;
	iB{1} = B{i_var};	% intermediate term
	Y_i_cont(:,i_var) = tf_ARX_pred(U_i(:,i_var), A, iB, dt(i_var)); % calculate increment vector
end
Y_i = sum(Y_i_cont,2);	% calculate increment at each step

	Y_p(Bmax,:) = Y(Bmax,:);
for i_sam = Bmax+1:n_sam
	Y_p(i_sam,1) = Y_p(i_sam-1,1) + Y_i(i_sam,1); % calculate the predicted output associated with each input
end	
