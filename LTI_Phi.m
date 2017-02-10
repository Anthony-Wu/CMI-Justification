function [Phi] = LTI_Phi(U, y_hat, nB, nA, dt)

[n_sam, n_varU] = size(U);	% number of samples & number of variables

Phi = zeros(n_sam,sum(nB)+nA);
Phi_u = zeros(n_sam,sum(nB));		% model input FIR part
Phi_u_int = zeros(n_sam,n_varU);	% intermediate matrix

if nA == 0;	% FIR
	bkmk = 0;
	for i_var = 1:n_varU
		Phi_u_int(:,i_var) = shiftdwn(U(:,i_var),dt(i_var));	% shift for dead time
		Phi(:,bkmk+1:bkmk+nB(i_var)) = spreaddwn(Phi_u_int(:,i_var),nB(i_var)); 	% spread to populate FIR
		bkmk = bkmk + nB(i_var);	% update bookmark for next variable
	end
elseif nA > 0; % ARX
	bkmk = 0;
	for i_var = 1:n_varU
		Phi_u_int(:,i_var) = shiftdwn(U(:,i_var),dt(i_var));	% shift for dead time
		Phi_u(:,bkmk+1:bkmk+nB(i_var)) = spreaddwn(Phi_u_int(:,i_var),nB(i_var)); 	% spread to populate FIR
		bkmk = bkmk + nB(i_var);	% update bookmark for next variable
	end
	Phi_y_int = spreaddwn(shiftdwn(y_hat,1),nA);
end

end