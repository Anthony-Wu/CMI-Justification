% CMI_FIR_pt1
% Justification for constrained model identification - FIR models
% part 1: Specify model

nB = [30 30 30];
K = [20 5 -10];
tau = [10 20 30];
dt = [30 10 20];

%nB = [60];
%K = [10];
%tau = [10];
%dt = [20];

n_var = size(nB,2);

[B] = FIR_coeffgen (nB, K, tau, dt);	% generate FIR coefficients