function [Y_p, Y_i] = tf_OE_pred_inc3 (U_i, F, B, dt)
% OE predictor (infinite step)
% MIMO compatible

% Inputs
% =======
% U_i = incremental input [nsam x nU]
% F = cell array of the AR coefficients [1 x nY]
% B = cell array of the FIR coefficients [nY x nU]
% dt = dead time of the input-output path [nY x nU]

[nB NA] = cellfun(@size, B); % nB is a matrix of the number of FIR coefficients for each input [nY x nU]
[nF NA] = cellfun(@size, F); % nF is a vector of the number of AR coefficients for each input [1 x nY]
nC = nB + nF; % total number of coefficients
Cmax = max(nC)';

[nsam nU] = size(U_i); % number of samples & number of input variable
[nsam nY] = size(Y); % number of samples & number of input variable

Y_p = zeros(nsam, nY); % predicted output
Y_i = zeros(nsam, nY); % predicted output increment	
Y_i_cont = zeros(nsam, nU*nY);	% predicted output increment associated with each input

for iY = 1:nY
	for iU = 1:nU;
		ivar = (iY - 1) * nU + iU;
		iB{1} = B{iY,iU};	% intermediate term
		iF{1} = F{iY};		% intermediate term
		Y_i_cont(:,ivar) = tf_OE_pred(U_i(:,iU), iF, iB, dt(iY,iU)); % calculate increment vector
	end
	
	Y_i(:,iY) = sum(Y_i_cont(:,((iY-1)*nU):ivar),2); % calculate increment at each step
	
	Y_p(1,iY) = Y_i(1,iY);
	for isam = 2:nsam;
		Y_p(i_sam,iY) = Y_p(i_sam-1,iY) + Y_i(i_sam,iY); % calculate the predicted output associated with each input
	end
end