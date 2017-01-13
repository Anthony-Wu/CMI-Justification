% CMI_FIR_pt2.m
% Justification for constrained model identification - FIR models
% part 2: Data Collection

% variables from earlier parts: n_var, n_coeff, K, tau, dt

% Define input steps
n_sam_t = 300; % number of samples in training dataset
n_sam_v = 500; % number of samples in training dataset

% Ut = ind_step(nB, dt);
Ut = zeros(n_sam_t, n_var);
Uv = ind_step(nB, dt);
% Uv = zeros(n_sam_v, n_var);

for i_var = 1:n_var
	Ut(:,i_var) = [RS(n_sam_t, 0, 1, 30, 60)];
	% Ut(size(Ut,1)+1:n_sam_t,i_var) = [RS(n_sam_t-size(Ut,1), 0, 1, 30, 60)];
	%Uv(:,i_var) = [RS(n_sam_v, 0, 1, 30, 60)];
	Uv(size(Uv,1)+1:n_sam_v,i_var) = [RS(n_sam_v-size(Uv,1), 0, 1, 30, 60)]; 
end

n_sam_v = size(Uv,1);

% Calculate output
[Yt, Yt_cont] = tf_FIR_pred (Ut, B, zeros(n_var,1));
[Yv] = tf_FIR_pred (Uv, B, zeros(n_var,1));

% Add noise
SNR_y = 10; % specify amount of noise to add to output signal (0 to switch off)
if SNR_y ~= 0;
	Yt_measured = awgn(Yt,SNR_y,'measured');
else
	Yt_measured = Yt;
end