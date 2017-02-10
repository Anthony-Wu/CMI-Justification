function [M_norm_hat] = FIM_LTI(U, theta_hat, nB, nA, dt_hat)
% Calculate FIM 
% Reference: Cai et al 2016 - Optimizing Zone Temperature Setpoint Excitation to Minimize Training Data for Data-Driven Dynamic Building Models
% eqn1: M_norm_hat(theta) = diag(theta) * PSI(theta)' * PSI(theta) * diag(theta)
% eqn2: PSI(theta) = [psi(1,theta)  psi(2,theta) ... psi(N,theta)]
% eqn3: psi(t,theta) = d/d(theta) {y_hat(t,theta)}

% For a linear, time-invariant (LTI) model
% eqn4: y_hat(t, theta) = Phi(t) * theta
% eqn5: d/d(theta) {y_hat(t), theta} = sum(Phi(t))

% Prerequisite variables: theta_hat, U, y
[n_sam, n_var] = size(U);

% convert theta_hat to vectors of B (and A) for model prediction
bkmk = 0;
for i_var = 1:n_var
	B_hat{i_var} = [zeros(dt_hat(i_var),1);[theta_hat(bkmk+1:bkmk+nB(i_var))]];
	bkmk = bkmk + nB(i_var);
end
A_hat{1} = theta_hat(bkmk+1:end);

% estimate output
y_hat = zeros(n_sam,1);
if nA == 0;
	[y_hat] = tf_FIR_pred (U, B_hat, dt_hat);
elseif nA > 0;
	[y_hat] = tf_ARX_pred (U, A_hat, B_hat, dt_hat);
end

% eqn2: determine the vector PSI
[Phi] = LTI_Phi(U, y_hat, nB, nA, dt_hat);
PSI = sum(Phi,2);	% eqn5

% eqn1: determine FIM, M_norm_hat(theta)
M_norm_hat = diag(theta_hat) * (PSI' * PSI) * diag(theta_hat);

end
