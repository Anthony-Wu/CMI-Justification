% 00 Estimated model parameters

% Currently assuming SISO FIR

B_hat{1} = [3,2,1,1]'; 	% used to estimate output
nB = length(B_hat{1});	% nB = number of FIR coefficients 
dt = 0; 				% dt = dead time

% 01 Design Input Trajectory
% ===

% Currently just randomised
n_sam = 100; % number of samples

U = RS(n_sam, 0, 1, 5, 10); % generates random step tests 
% RS(number of samples, min input magnitude, max input magnitude, min step legnth, max step legnth)

% 02 Estimate Output
% ===

% Y_hat = Phi * theta; Phi is the data matrix, theta is the model parameters

Phi = spreaddwn(U,nB);	% SISO FIR
theta = B_hat{1};		% SISO FIR

Y_hat = Phi * theta;

% Calculate FIM
% ===

% F = (delta y_hat / delta theta)' * Q * (delta y_hat / delta theta)

% For SISO FIR (delta y_hat / delta theta) = Phi
% Assume Q (weighting matrix) = Identity matric for now
% F = Phi' * Phi

F = Phi' * Phi;

% 03 Rate F
% ===
tog_F = 1; % 1 = A-criterion, 2 = E-criterion, 3 = ME-criterion (to be used later?)

F_A = min(trace(F)); % just using A-criterion for now

% Thoughts
% ===

% for SISO FIR...the coefficient values only affect the estimated output (more specifically whether or not it violates contraints)

% Alternative approach
% ===

% fun = @(U)trace((spreaddwn(U,nB))'*(spreaddwn(U,nB))); % basically compiling steage 02 and 03 as a function, and solve as a fmincon problem. 
% Input contraints (magntiude) can be easily done with inequlaity constraints A & b; 
% Input contraints (step length)...unless if you have fixed step lengths, this will be harder to do
% Output constraints...you would presumably set it as some diagonalised form for B{1} for A and b would be the max output magnitude

%...I think it may be possible to do this...but A is going be a huge matrix....

% U = fmincon(fun,x0,A,b)

