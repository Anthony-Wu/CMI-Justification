function [B, A] = QP_ARX_MISO (U, y, nB, dt, nA, con_A, con_b)
% uses ordinary least squares to identify the FIR coefficients
% U is the input signal
% y is the (single) output signal
% nB is a row vector for FIR coefficients for each input signal
% dt is a row vector for the dead time for each input signal
% nA is the number of autoregression coefficients
% con_A and con_b are the inequality constraints (set as [] if not used, i.e. unconstrained)

	[n_sam, n_var] = size(U);	% number of samples & number of variables
	
	% Define matrix sizes
	X = zeros(n_sam,sum(nB)+nA); 	% model input
	X_u_int = zeros(n_sam,n_var); 	% intermediate matrix
	X_u = zeros(n_sam,sum(nB));		% model input FIR part
	X_y = zeros(n_sam,nA); 			% model input ARX part
	
	% FIR part
	bkmk = 0;
	for i_var = 1:n_var
		X_u_int(:,i_var) = shiftdwn(U(:,i_var),dt(i_var));	% shift for dead time
		X_u(:,bkmk+1:bkmk+nB(i_var)) = spreaddwn(X_u_int(:,i_var),nB(i_var)); 	% spread to populate FIR
		bkmk = bkmk + nB(i_var);	% update bookmark for next variable
	end
	
	% ARX part
	X_y = spreaddwn(shiftdwn(y,1),nA); 
	
	% Combined X matrix
	X = [X_u, X_y]; 
	
	% Quadratic Programming
	paras = zeros(sum(nB)+nA,1);
	paras0 = zeros(sum(nB)+nA,1);
	
	H = 2*X'*X;
	f = -2*X'*y;
	
	if 	(isempty(con_A) + isempty(con_b)) > 0;	% if either part of the inequality constraints is not specified
		paras = quadprog(H,f,[],[],[],[],[],[],paras0);
	else
		paras = quadprog(H,f,con_A,con_b,[],[],[],[],paras0);
	end
	
	% reformat results to cell array
	bkmk = 0;
	for i_var = 1:n_var
		B{i_var} = [zeros(dt(i_var),1);[paras(bkmk+1:bkmk+nB(i_var))]];
		bkmk = bkmk + nB(i_var);
	end
		A{1} = -paras(bkmk+1:end);
end
