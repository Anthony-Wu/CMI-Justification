% CMI_FIR_pt1
% Justification for constrained model identification - FIR models
% part 1: Specify model

nB = [60 60 60];
K = [20 5 -10];
tau = [10 10 10];
dt = [30 10 20];

%nB = [60];
%K = [10];
%tau = [10];
%dt = [40];

n_var = size(nB,2);

[B] = FIR_coeffgen (nB, K, tau, dt);	% generate FIR coefficients