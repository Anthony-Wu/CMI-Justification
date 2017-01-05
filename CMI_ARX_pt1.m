% CMI_ARX_pt1.m

nB = [10 10 10];
K = [20 5 -10];
tau = [10 10 10];
dt = [30 10 20];

%nB = [60];
%K = [10];
%tau = [10];
%dt = [40];

n_var = size(nB,2);

[B] = FIR_coeffgen (nB, K, tau, dt);	% generate FIR coefficients

A{1} = [-0.2 -0.1]';