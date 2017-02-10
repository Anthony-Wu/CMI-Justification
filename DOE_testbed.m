% Estimate output
n_var = 2;
B{1} = [4 3 1 1 1]';
B{2} = [2 1 1]';
nB = [5 3];
nA = 0;
dt = [0 0];

% Define input steps
n_sam_t = 400; % number of samples in training dataset

U1 = zeros(n_sam_t, n_var);	
for i_var = 1:n_var;
	U1(:,i_var) = [RS(n_sam_t, 0, 10, 1, 5)];
end

% Estimate output
y1=tf_FIR_pred (U1, B, dt);

dt_hat = [0 0];

% test 1
theta_hat_1 = [4 3 1 1 1 2 1 1]';
[M_1] = FIM_LTI(U1, theta_hat_1, nB, nA, dt_hat);


% test 2
theta_hat_2 = [3 2 1 1 1 1 1 1]';
[M_2] = FIM_LTI(U1, theta_hat_2, nB, nA, dt_hat);