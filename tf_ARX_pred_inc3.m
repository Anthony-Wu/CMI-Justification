function [Yi_p Yi_FIR Yi_FIRc] = tf_ARX_pred_inc3 (Ui, A, B, dt, Yi, nk)
% n-th sample ahead predictor
% MIMO compatible

% Inputs
% =======
% Ui = incremental input [nsam x nU]
% A = cell array of the AR coefficients [1 x nY]
% B = cell array of the FIR coefficients [nY x nU]
% dt = dead time of the input-output path [nY x nU]
% Y = measured output (absolute) [nsam x nY]
% nk = n-th sample ahead [scalar]

[nB NA] = cellfun(@size, B); % nB is a matrix of the number of FIR coefficients for each input [nY x nU]
[nA NA] = cellfun(@size, A); % nA is a vector of the number of AR coefficients for each input [1 x nY]
nC = nB + nA; % total number of coefficients
Cmax = max(nC)';

[nsam nU] = size(Ui); % number of samples & number of input variable
[nsam nY] = size(Yi); % number of samples & number of input variable

Yi_p = zeros(nsam, nY); % predicted output increment
Yi_FIR = zeros(nsam, nY); % predicted output increment
Yi_FIRc = zeros(nsam, nU*nY); % predicted output increment associated with each input

% FIR part (this part will stay the same since only inputs are taken)
for iY = 1:nY
	for iU = 1:nU;
		ivar = (iY-1)*nU+iU;
		
		iB{1} = B{iY,iU}; % intermediate term
		Yi_FIRc(:,ivar) = tf_FIR_pred (Ui(:,iU), iB, dt(iY,iU));
	end
	Yi_FIR(:,iY) = sum(Yi_FIRc(:,((iY-1)*nU)+1:ivar),2); % calculate increment at each step
end

% AR part (this will change a bit depending on the starting point

% nk-step ahead prediction
% you have y(t), u(t) and prior, you want to estimate y(t+nk)
for iY = 1:nY
	inA = nA(iY);
	Yi_p(1:Cmax(iY),iY) = Yi(1:Cmax(iY),iY); % not enough samples to fully populate model, do not estimate
	for isam = Cmax+1:nsam-nk
		Ypi_nk = zeros(nk+inA,1);
		Ypi_nk = [Yi(isam-inA+1:isam)' Yi_FIR(isam+1:isam+nk)']';
		for i = 1:nk;
			Ypi_nk(inA+i) = -A{iY}' * flipud(Ypi_nk(i:i+inA-1)) + Ypi_nk(inA+i);
		end
		
		if isam == Cmax+1;
			Yi_p(isam:isam+nk-1,iY) = Ypi_nk(nA:end-1);
		end
		
		Yi_p(isam+nk,iY) = Ypi_nk(end);
	end
end