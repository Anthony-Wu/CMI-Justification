function [F_A, Phi] = DOE_FA(paras, nS, B, A, dt)
% paras = input trajectory: [Ui_1(1),...,Ui_1(nS),...,Ui_nU(1),...,Ui_nU(nS)]'
% nS = number of samples in dataset

nU = length(paras) / nS;	% number of input signals

Ui = zeros(nS, nU);

% assign paras to input matrix
for iU = 1:nU;
	U(:,iU) = paras((iU-1)*nU+1:(iU-1)*nU+nS);
end

Ui = abs2inc(U);

% calculate the output
[Y, Yi] = tf_OE_pred_inc3 (Ui, A, B, dt);
nY = size(Y, 2);

% calculate the FIM (A criterion)
% FIM = Phi' * Q * Phi
% Q assumed as identity matrix for now

% check total nA and nB

tot_nA = 0;
for iY = 1:nY;
	tot_nA = tot_nA + length(A{iY});
end

tot_nB = 0;
for iY = 1:nY;
	for iU = 1:nU;
		tot_nB = tot_nB + length(B{iY,iU});
	end
end

Phi = zeros(nS,tot_nB + tot_nA); % model input
Phi_U_int = zeros(nS,nU); % intermediate matrix
Phi_U = zeros(nS,tot_nB); % model input FIR part
Phi_Y = zeros(nS,tot_nA); % model input ARX part

% Calcualte Phi
bkmk_B = 0;
for iY = 1:nY;
	for iU = 1:nU;
		nB = length(B{iY,iU});
		Phi_U_int(:,iU) = shiftdwn(U(:,iU),dt(iY,iU));	% shift for dead time
		Phi_U(:,bkmk_B+1:bkmk_B+nB) = spreaddwn(Phi_U_int(:,iU),nB); % spread to populate FIR
		bkmk_B = bkmk_B + nB;	% update bookmark for next variable
	end
end

bkmk_A = 0;
for iY = 1:nY;
	nA = length(A{iY});
	Phi_Y(:,bkmk_A+1:bkmk_A+nA) = -spreaddwn(shiftdwn(Y(:,iY),1),nA); 
	bkmk_A = bkmk_A + nA;
end

Phi = [Phi_U Phi_Y];

F = Phi' * Phi; % FIM

% A-criterion

F_A = trace (F^(-1));

end