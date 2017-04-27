% 00 Input Summary
% ----

% 01: Set up Synthetic Process
% ----

nA = zeros(nY,1);
for iY = 1:nY;
	nA(iY) = length(c_A{iY})';
end

% shortcut to convert single-output to multiple-output if parameters are the same for all outputs
[nB] = sct_SO2MO(nY, nB); [K] = sct_SO2MO(nY, K); [tau] = sct_SO2MO(nY, tau); [dt] = sct_SO2MO(nY, dt);

% rearrange FIR components to a single row vector
[c_B, ref_table] = FIR_coeffgen2 (nU, nY, nB, K, tau, dt); % generate FIR coefficients

% 02: Design Input trajectory and obtain data
% ----
n_paras = sum(sum(nB)) + sum(sum(nA)) + sum(sum(dt)); % number of parameters to identify

data_para_ratio = nSamt/n_paras; % calculate the data to parameter ratio

sLength = max(nB) + max(nA)*ones(1,nU); % step length needed to capture all coefficients

Ut = zeros(nSamt,nU);	% generate inputs (training)
nSamt_iU = nSamt/nU;
bkmk = 0;
for iU = 1:nU
	Ut(bkmk+1:bkmk+nSamt_iU,iU) = [RS3(nSamt_iU, 0, 1, floor(min(sLf)*sLength(iU)), ceil(max(sLf)*sLength(iU)))];
	bkmk = bkmk + nSamt_iU;
end

Uv = ind_step2(nB, dt);
nSamv = size(Uv,1);

%Uv = zeros(nSamv,nU);	% generate inputs (training)
%for iU = 1:nU
	%Uv(:,iU) = [RS(nSamv, 0, 1, floor(min(sLf)*sLength(iU)), ceil(max(sLf)*sLength(iU)))];
%end

% TODO: rewrite this as a single function later
Yt = zeros(nSamt,nY);	% generate outputs
Yv = zeros(nSamv,nY);
for iY = 1:nY;	% generate output
	c_A_int{1} = c_A{iY};
	for iU = 1:nU;
		c_B_int{iU} = c_B{nU*(iY-1)+iU};
	end
	Yt(:,iY) = tf_OE_pred(Ut, c_A_int, c_B_int, dt(iY,:));
	Yv(:,iY) = tf_OE_pred(Uv, c_A_int, c_B_int, dt(iY,:));
end

% add noise to output
Yv_meas = Yv;
if 	SNR_y == 0;
	Yt_meas = Yt;
else
	Yt_meas = awgn(Yt,SNR_y,'measured');
end

% 02: Data pre-treatment
% TODO: merge as single function
Yt_f = filter(ones(1,mov_avg)/mov_avg, 1, Yt_meas);
Yv_f = filter(ones(1,mov_avg)/mov_avg, 1, Yv_meas);

% autoscaled, absolute
[Ut_as, Ut_m, Ut_std] = ascale(Ut); 
[Yt_fas, Yt_m, Yt_std] = ascale(Yt_f);
Uv_as = (Uv - Ut_m) ./ Ut_std;
Yv_fas = (Yv_f - Yt_m) ./ Yt_std;

% autoscaled, incremental
Ut_asi 	= abs2inc(Ut_as); Yt_fasi = abs2inc(Yt_fas);
Uv_asi 	= abs2inc(Uv_as); Yv_fasi = abs2inc(Yv_fas);

% measured, incremental
Ut_i  = abs2inc(Ut); Yt_fi = abs2inc(Yt_f);
Uv_i  = abs2inc(Uv); Yv_fi = abs2inc(Yv_f);

% export to WaterMV
if tog_as == 1;
	MV_export = [Ut_as, shiftdwn(Yt_fas,1), Yt_fas];
else
	MV_export = [Ut, shiftdwn(Yt_f,1), Yt_f];
end

% 03: Model Identification
% ----

K_sign_in = (K./abs(K)); % actual steady state gain direction
A_sign_in = zeros(1,nY);
for iY = 1:nY;
	A_sign_in(iY) = sum(c_A{1})./abs(sum(c_A{1}));
end

% M1: UNC2: Unconstrained (Matlab ARX)
% ----
[B_M1 A_M1] = QP_ARX_MISO (Ut_asi, Yt_fasi, nB_in, dt_in, nA_in, [], []); % model identification

Yt_fas_M1 = tf_ARX_pred_inc2 (Ut_asi, A_M1, B_M1, dt_in, Yt_fas); % est. output (training, absolute)
Yv_fas_M1 = tf_ARX_pred_inc2 (Uv_asi, A_M1, B_M1, dt_in, Yv_fas); % est. output (validation, absolute)

Yt_fasi_M1 = abs2inc(Yt_fas_M1); % est. output (training, incremental)
Yv_fasi_M1 = abs2inc(Yv_fas_M1); % est. output (validation, incremental)

Yt_f_M1 = zeros(nSamt,nY); Yv_f_M1 = zeros(nSamv,nY); Yt_fi_M1 = zeros(nSamt,nY); Yv_fi_M1 = zeros(nSamv,nY);
for iY = 1:nY
	Yt_f_M1(:,iY) = Yt_fas_M1(:,iY) .* Yt_std(iY) + Yt_m(iY);	
	Yv_f_M1(:,iY) = Yv_fas_M1(:,iY) .* Yt_std(iY) + Yt_m(iY);
	Yt_fi_M1(:,iY) = abs2inc(Yt_f_M1(:,iY));
	Yv_fi_M1(:,iY) = abs2inc(Yv_f_M1(:,iY));
end

% calculate the RMSE
RMSE_Yt_fas_M1 = stat_RMSE2(Yt_fas_M1,Yt_fas,sLength);
RMSE_Yv_fas_M1 = stat_RMSE2(Yv_fas_M1,Yv_fas,sLength);
RMSE_Yt_fasi_M1 = stat_RMSE2(Yt_fasi_M1,Yt_fasi,sLength);
RMSE_Yv_fasi_M1 = stat_RMSE2(Yv_fasi_M1,Yv_fasi,sLength);
RMSE_Yt_f_M1 = stat_RMSE2(Yt_f_M1,Yt_f,sLength);
RMSE_Yv_f_M1 = stat_RMSE2(Yv_f_M1,Yv_f,sLength);
RMSE_Yt_fi_M1 = stat_RMSE2(Yt_fi_M1,Yt_fi,sLength);
RMSE_Yv_fi_M1 = stat_RMSE2(Yv_fi_M1,Yv_fi,sLength);

% M2 : CON01 Direction specified, not min phase
% ----
K_sign_M2 = [K_sign_in A_sign_in]; 		% specify sign direction vector
min_phase_M2 = [zeros(1,nU) zeros(1,nY)];	% 0 = not minimum phase

% constraints
nZ_M2 = [nB_in nA_in]; 
[con_A_M2 con_b_M2] = CMI_con(nZ_M2, K_sign_M2, min_phase_M2);	% set up inequality constraint matrices


[B_M2 A_M2] = QP_ARX_MISO (Ut_asi, Yt_fasi, nB_in, dt_in, nA_in, con_A_M2, con_b_M2); % model identification

Yt_fas_M2 = tf_ARX_pred_inc2 (Ut_asi, A_M2, B_M2, dt_in, Yt_fas); % est. output (training, absolute)
Yv_fas_M2 = tf_ARX_pred_inc2 (Uv_asi, A_M2, B_M2, dt_in, Yv_fas); % est. output (validation, absolute)

Yt_fasi_M2 = abs2inc(Yt_fas_M2); % est. output (training, incremental)
Yv_fasi_M2 = abs2inc(Yv_fas_M2); % est. output (validation, incremental)

Yt_f_M2 = zeros(nSamt,nY); Yv_f_M2 = zeros(nSamv,nY); Yt_fi_M2 = zeros(nSamt,nY); Yv_fi_M2 = zeros(nSamv,nY);
for iY = 1:nY
	Yt_f_M2(:,iY) = Yt_fas_M2(:,iY) .* Yt_std(iY) + Yt_m(iY);	
	Yv_f_M2(:,iY) = Yv_fas_M2(:,iY) .* Yt_std(iY) + Yt_m(iY);
	Yt_fi_M2(:,iY) = abs2inc(Yt_f_M2(:,iY));
	Yv_fi_M2(:,iY) = abs2inc(Yv_f_M2(:,iY));
end

% calculate the RMSE
RMSE_Yt_fas_M2 = stat_RMSE2(Yt_fas_M2,Yt_fas,sLength);
RMSE_Yv_fas_M2 = stat_RMSE2(Yv_fas_M2,Yv_fas,sLength);
RMSE_Yt_fasi_M2 = stat_RMSE2(Yt_fasi_M2,Yt_fasi,sLength);
RMSE_Yv_fasi_M2 = stat_RMSE2(Yv_fasi_M2,Yv_fasi,sLength);
RMSE_Yt_f_M2 = stat_RMSE2(Yt_f_M2,Yt_f,sLength);
RMSE_Yv_f_M2 = stat_RMSE2(Yv_f_M2,Yv_f,sLength);
RMSE_Yt_fi_M2 = stat_RMSE2(Yt_fi_M2,Yt_fi,sLength);
RMSE_Yv_fi_M2 = stat_RMSE2(Yv_fi_M2,Yv_fi,sLength);

% M3 : CON02 Direction specified, min phase (FIR ONLY)
% ----
K_sign_M3 = [K_sign_in A_sign_in]; 		% specify sign direction vector
min_phase_M3 = [ones(1,nU) zeros(1,nY)];	% 0 = not minimum phase

% constraints
nZ_M3 = [nB_in nA_in]; 
[con_A_M3 con_b_M3] = CMI_con(nZ_M3, K_sign_M3, min_phase_M3);	% set up inequality constraint matrices


[B_M3 A_M3] = QP_ARX_MISO (Ut_asi, Yt_fasi, nB_in, dt_in, nA_in, con_A_M3, con_b_M3); % model identification

Yt_fas_M3 = tf_ARX_pred_inc2 (Ut_asi, A_M3, B_M3, dt_in, Yt_fas); % est. output (training, absolute)
Yv_fas_M3 = tf_ARX_pred_inc2 (Uv_asi, A_M3, B_M3, dt_in, Yv_fas); % est. output (validation, absolute)

Yt_fasi_M3 = abs2inc(Yt_fas_M3); % est. output (training, incremental)
Yv_fasi_M3 = abs2inc(Yv_fas_M3); % est. output (validation, incremental)

Yt_f_M3 = zeros(nSamt,nY); Yv_f_M3 = zeros(nSamv,nY); Yt_fi_M3 = zeros(nSamt,nY); Yv_fi_M3 = zeros(nSamv,nY);
for iY = 1:nY
	Yt_f_M3(:,iY) = Yt_fas_M3(:,iY) .* Yt_std(iY) + Yt_m(iY);	
	Yv_f_M3(:,iY) = Yv_fas_M3(:,iY) .* Yt_std(iY) + Yt_m(iY);
	Yt_fi_M3(:,iY) = abs2inc(Yt_f_M3(:,iY));
	Yv_fi_M3(:,iY) = abs2inc(Yv_f_M3(:,iY));
end

% calculate the RMSE
RMSE_Yt_fas_M3 = stat_RMSE2(Yt_fas_M3,Yt_fas,sLength);
RMSE_Yv_fas_M3 = stat_RMSE2(Yv_fas_M3,Yv_fas,sLength);
RMSE_Yt_fasi_M3 = stat_RMSE2(Yt_fasi_M3,Yt_fasi,sLength);
RMSE_Yv_fasi_M3 = stat_RMSE2(Yv_fasi_M3,Yv_fasi,sLength);
RMSE_Yt_f_M3 = stat_RMSE2(Yt_f_M3,Yt_f,sLength);
RMSE_Yv_f_M3 = stat_RMSE2(Yv_f_M3,Yv_f,sLength);
RMSE_Yt_fi_M3 = stat_RMSE2(Yt_fi_M3,Yt_fi,sLength);
RMSE_Yv_fi_M3 = stat_RMSE2(Yv_fi_M3,Yv_fi,sLength);

% M4 : CON03 Direction specified, min phase (FIR & AR)
% ----
K_sign_M4 = [K_sign_in A_sign_in]; 		% specify sign direction vector
min_phase_M4 = [ones(1,nU) ones(1,nY)];	% 0 = not minimum phase

% constraints
nZ_M4 = [nB_in nA_in]; 
[con_A_M4 con_b_M4] = CMI_con(nZ_M4, K_sign_M4, min_phase_M4);	% set up inequality constraint matrices


[B_M4 A_M4] = QP_ARX_MISO (Ut_asi, Yt_fasi, nB_in, dt_in, nA_in, con_A_M4, con_b_M4); % model identification

Yt_fas_M4 = tf_ARX_pred_inc2 (Ut_asi, A_M4, B_M4, dt_in, Yt_fas); % est. output (training, absolute)
Yv_fas_M4 = tf_ARX_pred_inc2 (Uv_asi, A_M4, B_M4, dt_in, Yv_fas); % est. output (validation, absolute)

Yt_fasi_M4 = abs2inc(Yt_fas_M4); % est. output (training, incremental)
Yv_fasi_M4 = abs2inc(Yv_fas_M4); % est. output (validation, incremental)



Yt_f_M4 = zeros(nSamt,nY); Yv_f_M4 = zeros(nSamv,nY); Yt_fi_M4 = zeros(nSamt,nY); Yv_fi_M4 = zeros(nSamv,nY);
for iY = 1:nY
	Yt_f_M4(:,iY) = Yt_fas_M4(:,iY) .* Yt_std(iY) + Yt_m(iY);	
	Yv_f_M4(:,iY) = Yv_fas_M4(:,iY) .* Yt_std(iY) + Yt_m(iY);
	Yt_fi_M4(:,iY) = abs2inc(Yt_f_M4(:,iY));
	Yv_fi_M4(:,iY) = abs2inc(Yv_f_M4(:,iY));
end

% calculate the RMSE
RMSE_Yt_fas_M4 = stat_RMSE2(Yt_fas_M4,Yt_fas,sLength);
RMSE_Yv_fas_M4 = stat_RMSE2(Yv_fas_M4,Yv_fas,sLength);
RMSE_Yt_fasi_M4 = stat_RMSE2(Yt_fasi_M4,Yt_fasi,sLength);
RMSE_Yv_fasi_M4 = stat_RMSE2(Yv_fasi_M4,Yv_fasi,sLength);
RMSE_Yt_f_M4 = stat_RMSE2(Yt_f_M4,Yt_f,sLength);
RMSE_Yv_f_M4 = stat_RMSE2(Yv_f_M4,Yv_f,sLength);
RMSE_Yt_fi_M4 = stat_RMSE2(Yt_fi_M4,Yt_fi,sLength);
RMSE_Yv_fi_M4 = stat_RMSE2(Yv_fi_M4,Yv_fi,sLength);



% align B coefficient vectors to the same size (for later comparison)
for i = 1:nU*nY;
	iY = ref_table(i,2); iU = ref_table(i,3); 
	c_B_plot{i} = [c_B{iY,iU}; zeros(dt_in(iY,iU)+nB_in(iY,iU)-length(c_B{iY,iU}),1)];
	if 	tog_as == 1;
		c_B_mv_plot{iY,iU} = filter(ones(1,mov_avg)/mov_avg,1,c_B_plot{iY,iU}*(Ut_std(iU)/Yt_std));
	else
		c_B_mv_plot{iY,iU} = filter(ones(1,mov_avg)/mov_avg,1,c_B_plot{iY,iU});
	end
end

% convert FIR coefficients from as scale to eng scale
for i = 1:nU*nY;
	iY = ref_table(i,2); iU = ref_table(i,3);
	c_B_M1{iY,iU} = B_M1{iY,iU}.*Yt_std(iY)./Ut_std(iU); 
	c_B_M2{iY,iU} = B_M2{iY,iU}.*Yt_std(iY)./Ut_std(iU); 
	c_B_M3{iY,iU} = B_M3{iY,iU}.*Yt_std(iY)./Ut_std(iU); 	
	c_B_M4{iY,iU} = B_M4{iY,iU}.*Yt_std(iY)./Ut_std(iU);
	c_B_sum{iY,iU} = [c_B_mv_plot{iY,iU}, c_B_M1{iY,iU}, c_B_M2{iY,iU}, c_B_M3{iY,iU}, c_B_M4{iY,iU}];
end

Yt_f_sum = [Yt_f, Yt_f_M1, Yt_f_M2, Yt_f_M3, Yt_f_M4];
Yv_f_sum = [Yv_f, Yv_f_M1, Yv_f_M2, Yv_f_M3, Yv_f_M4];
Yt_fi_sum = [Yt_fi, Yt_fi_M1, Yt_fi_M2, Yt_fi_M3, Yt_fi_M4];
Yv_fi_sum = [Yv_fi, Yv_fi_M1, Yv_fi_M2, Yv_fi_M3, Yv_fi_M4];

RMSE_sum = [RMSE_Yt_f_M1 RMSE_Yt_f_M2 RMSE_Yt_f_M3 RMSE_Yt_f_M4;...
			RMSE_Yv_f_M1 RMSE_Yv_f_M2 RMSE_Yv_f_M3 RMSE_Yv_f_M4;...
			RMSE_Yt_fi_M1 RMSE_Yt_fi_M2 RMSE_Yt_fi_M3 RMSE_Yt_fi_M4;...
			RMSE_Yv_fi_M1 RMSE_Yv_fi_M2 RMSE_Yv_fi_M3 RMSE_Yv_fi_M4;...
			RMSE_Yt_fas_M1 RMSE_Yt_fas_M2 RMSE_Yt_fas_M3 RMSE_Yt_fas_M4;...
			RMSE_Yv_fas_M1 RMSE_Yv_fas_M2 RMSE_Yv_fas_M3 RMSE_Yv_fas_M4;...
			RMSE_Yt_fasi_M1 RMSE_Yt_fasi_M2 RMSE_Yt_fasi_M3 RMSE_Yt_fasi_M4;...
			RMSE_Yv_fasi_M1 RMSE_Yv_fasi_M2 RMSE_Yv_fasi_M3 RMSE_Yv_fasi_M4];

	gain_B = K./abs(K);
for i = 1:nU*nY;
	iY = ref_table(i,2); iU = ref_table(i,3);
	gain_B_M1{iY,iU} = sum(B_M1{iY,iU})/abs(sum(B_M1{iY,iU})); 
	gain_B_M2{iY,iU} = sum(B_M1{iY,iU})/abs(sum(B_M1{iY,iU}));  
	gain_B_M3{iY,iU} = sum(B_M1{iY,iU})/abs(sum(B_M1{iY,iU}));  	
	gain_B_M4{iY,iU} = sum(B_M1{iY,iU})/abs(sum(B_M1{iY,iU})); 
	gain_B_sum{iY,iU} = [gain_B_M1{iY,iU}, gain_B_M2{iY,iU}, gain_B_M3{iY,iU}, gain_B_M4{iY,iU}];
end
			
% save(filename)
%for i_save = 1:10
%    file_name = sprintf('case study_#%d.mat',i_save);
%    save(file_name);
%end
