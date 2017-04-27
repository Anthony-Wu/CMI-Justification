% CMI MC_Main

% 01 Randomise process characteristics
% =====================================
n_exp = 100;	% number of experiments
nU = 2; nY = 1;	% number of inputs and outputs

EXP_summary = zeros(n_exp,9); % summary table

for i_exp = 1:n_exp;

	% characteristic ranges
	rng_nB = [5 10]; rng_K = [-10 10]; rng_tau = [1 10]; rng_dt = [0 10];	rng_SNR = [10 30]; sLf = [1 1.5]; rng_mv = [1 10]; rng_AR = [-0.9 0.5]; rng_nSamt_iU = [50 150];% nB = number of FIR coefficients, K = gain, tau = time constant, dt = dead time, SNR = signal to noise ratio added to the output sLf = step length factor, mv = moving average window, AR = 1st autoregression coefficient, number of samples in training data
	
	nB = ones(nY, nU); K = zeros(nY, nU); tau = ones(nY, nU); dt = zeros(nY, nU);
	
	nB = round(rng_nB(1)+(rng_nB(2)-rng_nB(1)).*rand(nY,nU));
	tau = round(rng_tau(1)+(rng_tau(2)-rng_tau(1)).*rand(nY,nU));
	dt = round(rng_dt(1)+(rng_dt(2)-rng_dt(1)).*rand(nY,nU));	
	c_A{1} = round(rng_AR(1)+(rng_AR(2)-rng_AR(1)).*rand(1),1);	
	SNR_y = round(rng_SNR(1)+(rng_SNR(2)-rng_SNR(1)).*rand(1));	
	mov_avg = round(rng_mv(1)+(rng_mv(2)-rng_mv(1)).*rand(1)); 
		
	% process gain (process gain = 0 messes up the current algorithm, so if the min is -ve but max is +ve it can have problems)
	if 	rng_K(1) * rng_K(2) < 0; 
		for iU = 1:nU;
			rng = round(rand(1));
			if rng == 0;
				K(iU) = floor(rng_K(1).*rand(1));
			elseif rng == 1;
				K(iU) = ceil(rng_K(2).*rand(1));
			end
		end
	else
		K = round(rng_K(1)+(rng_K(2)-rng_K(1)).*rand(nY,nU));
	end
	
	tog_as = 1; tog_inc = 1; 
	
	%Dataset sizes
	nSamt_iU = round(rng_nSamt_iU(1)+(rng_nSamt_iU(2)-rng_nSamt_iU(1)).*rand(1)); % number of samples for step test per signal
	nSamt = nU * nSamt_iU; 	
	
	% model identifiication
	nB_in = [nB+dt]; dt_in = zeros(nY,nU); nA_in = [size(c_A{1},1)];
	
	CMI_EXP;
	
	EXP_summary(i_exp,:) = [i_exp, nSamt, data_para_ratio, SNR_y, c_A{1}, RMSE_Yt_fi_M1, RMSE_Yt_fi_M3, RMSE_Yv_fi_M1, RMSE_Yv_fi_M3];
	
	file_name_full = sprintf('fullworkspaceH_%d.mat',i_exp);
	save(file_name_full);
	
	file_sam_short = sprintf('summaryH_%d.mat',i_exp);
	save(file_sam_short,'c_B_sum','Yt_f_sum','Yv_f_sum','Yt_fi_sum','Yv_fi_sum','RMSE_sum');
	
	
	
	
end

% save(filename)
%for i_save = 1:10
%    file_name = sprintf('case study_#%d.mat',i_save);
%    save(file_name);
%end