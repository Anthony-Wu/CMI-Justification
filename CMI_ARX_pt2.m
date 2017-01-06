% CMI_ARX_pt2.m
% Justification for constrained model identification - ARX models
% part 2: Data Collection

% variables from earlier parts: n_var, n_coeff, K, tau, dt

n_sam_t = 200; % number of samples in training dataset
n_sam_v = 300; % number of samples in validation dataset

Ut = zeros(n_sam_t, n_var);
Uv = zeros(n_sam_v, n_var);	

toggle_Y_noise = 0; % 1 = add noise to output, 0 = do not add noise
SNR = 10;	% specify level of noise to add (only applies if toggle is set)

%Define input steps
for i_var = 1:n_var
	Ut(:,i_var) = RS(n_sam_t, 0, 1, floor(1*nB(i_var)), ceil(2*nB(i_var)));
	Uv(:,i_var) = RS(n_sam_v, 0, 1, nB(i_var), nB(i_var));
end

% Calculate output
[Yt, Yt_cont] = tf_ARX_pred (Ut, A, B, zeros(n_var,1));
[Yv] = tf_ARX_pred (Uv, A, B, zeros(n_var,1));

% Add noise
if toggle_Y_noise == 1;
	[Yt_n] = awgn(Yt,SNR,'measured');
	[Yt_measured] = Yt_n;
else
	[Yt_measured] = Yt;
end