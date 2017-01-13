% CMI_ARX_pt1.m

nB = [60 80 70];
K = [20 5 -15];
tau = [20 10 40];
dt = [25 5 15];

%nB = [60];
%K = [10];
%tau = [10];
%dt = [20];

n_var = size(nB,2);

[B] = FIR_coeffgen (nB, K, tau, dt);	% generate FIR coefficients
%B{1} = [7 5 3 2 1]';

A{1} = [-0.2 -0.1 -0.1]';