% CMI_FIR_pt4.m
% Justification for constrained model identification - FIR models
% part 4: Model identification

nB_in = [50 50 50];
for i_var = 1:n_var;
	B_true{i_var} = [B{i_var}; zeros(nB_in(i_var) - length(B{i_var}),1)]; % erase this is not longer needed, replaced by mv
	if toggle_ascale == 1;
		if 	toggle_filter == 1; 	% moving average and no auto-scaling used
			B_true_mv_as{i_var} = [B_mv_as{i_var}; zeros(nB_in(i_var) - length(B{i_var}),1)];
			B_true_as{i_var} = [B_as{i_var}; zeros(nB_in(i_var) - length(B{i_var}),1)];			
		elseif toggle_filter == 0; 	% auto-scaling used, but no moving average
			B_true_as{i_var} = [B_as{i_var}; zeros(nB_in(i_var) - length(B{i_var}),1)];
		end
	elseif 	toggle_filter == 1;		% moving average used, but no auto-scaling
		B_true_mv{i_var} = [B_mv{i_var}; zeros(nB_in(i_var) - length(B{i_var}),1)];	
	end
end

nA_in = 2;
A_true{1} = [A{1}; zeros(nA_in(1) - length(A{1}),1)];

toggle_inc = 1; 	% 0 = use absolute model, 1 = use incremental model
if toggle_inc == 1;
	Ut_model_i = abs2inc(Ut_model);
	Yt_model_i = abs2inc(Yt_model);
	Uv_model_i = abs2inc(Uv_model);
	Yv_model_i = abs2inc(Yv_model);	
end

if toggle_ascale == 1;	% Y_val is the output used for validation in Eng units (if moving average is used, it can skew the results)
	Yt_val = Yt_model.*Yt_std+Yt_m;
	Yv_val = Yv_model.*Yt_std+Yt_m;
else
	Yt_val = Yt_model;
	Yv_val = Yv_model;
end		
		
% M1: UNC1: Unconstrained using WaterMV
% Parameter estimation carried out in WaterMV, coefficients c/p here
B_M1nz{1} = [-3.522005E-7	-2.9537256E-7	-2.4728087E-7	-7.4914915E-8	-1.4716596E-7	-1.9787508E-7	5.6318083E-8	1.7204042E-7	4.968644E-7	3.6311468E-7	1.5447766E-7	-1.9190325E-7	1.4244865E-9	5.193836E-7	7.258231E-7	5.0774446E-7	8.058068E-7	8.006897E-7	6.6031913E-7	8.0789937E-7	6.6456624E-7	5.913621E-7	7.5268383E-7	1.3702745E-6	1.316995E-6	1.3792649E-6	1.5323054E-6	1.3996022E-6	1.5824663E-6	0.17074963	0.12974167	0.09194646	0.07563061	0.06437585	0.056681916	0.050568826	0.045456342	0.040998664	0.03704098	-0.029322889	-0.017434401	-0.006417558	-0.0030257558	-0.0012454562	-5.502521E-4	-2.334995E-4	-1.00897385E-4	-4.2616368E-5	-1.781542E-5	-6.9715506E-6]';	% c/p coefficients from WaterMV as appropriate
B_M1nz{2} = [1.5623965E-7	2.5196988E-7	1.7826625E-7	-1.6691494E-7	-2.0464809E-7	4.3818753E-9	1.9285628E-7	3.3627487E-7	3.076991E-8	0.042688817	0.03243659	0.02298721	0.01890826	0.01609411	0.01417069	0.012642338	0.011363721	0.010249302	0.009259719	-0.007331802	-0.0043593883	-0.0016054588	-7.5696135E-4	-3.120041E-4	-1.3816524E-4	-5.8753823E-5	-2.5624844E-5	-1.1152956E-5	-5.4303027E-6	-2.593879E-6	-1.5233545E-6	-7.01108E-7	-5.470957E-7	-4.015664E-7	-3.896105E-8	-2.7073952E-7	-1.2310782E-7	8.3230496E-8	4.538303E-8	-9.683473E-8	-4.875315E-7	-4.0134125E-7	-4.011392E-7	-4.324527E-8	1.1676137E-7	1.6064482E-7	3.9354418E-9	5.7996676E-8	1.2635216E-7	-2.5209104E-7]';
B_M1nz{3} = [8.8003425E-8	3.7003707E-7	3.3904678E-7	5.818518E-7	6.3974807E-7	6.9289314E-7	7.3613154E-7	7.089389E-7	6.790159E-7	4.6268366E-7	5.766591E-7	7.054822E-7	8.504472E-7	7.633454E-7	1.1323505E-6	1.2282349E-6	1.0426762E-6	1.0752252E-6	1.2536007E-6	-0.085270636	-0.064790696	-0.045915812	-0.037768144	-0.032147665	-0.028305406	-0.025252193	-0.022699224	-0.02047302	-0.01849636	0.014645749	0.008708139	0.0032063436	0.0015119048	6.227476E-4	2.756567E-4	1.1743641E-4	5.1134928E-5	2.1625283E-5	8.757394E-6	3.6722572E-6	1.3149735E-6	2.3460905E-7	3.3929854E-8	-3.591528E-8	-4.3598916E-7	-4.5711968E-7	-4.82455E-7	-5.4154685E-7	-2.5697423E-7	2.6283823E-7]';

A_M1nz{1} = [-0.124965094	-0.22004129]';
B_M1{1} = [[B_M1nz{1}]' zeros(1,nB_in(1) - length(B_M1nz{1}))]';
B_M1{2} = [[B_M1nz{2}]' zeros(1,nB_in(2) - length(B_M1nz{2}))]';
B_M1{3} = [[B_M1nz{3}]' zeros(1,nB_in(3) - length(B_M1nz{3}))]';
A_M1{1} = [[A_M1nz{1}]' zeros(1,nA_in(1) - length(A_M1nz{1}))]';

dt_M1 = zeros(1,n_var);			% specify a priori dead-time

% Output Prediction
if toggle_inc == 1;		% incremental model
	% WaterMV's notation is slightly different and needs a shift in coefficeints by 1 to translate over
	for i_var = 1:n_var
		B_M1{i_var} = [0; [B_M1{i_var}(1:end-1)]];
	end
	[Yt_model_M1] = tf_ARX_pred_inc (Ut_model_i, A_M1, B_M1, dt_M1, Yt_model);
	[Yv_model_M1] = tf_ARX_pred_inc (Uv_model_i, A_M1, B_M1, dt_M1, Yv_model);
elseif toggle_inc == 0;	% absolute model
	[Yt_model_M1] = tf_ARX_pred (Ut_model, A_M1, B_M1, dt_M1);
	[Yv_model_M1] = tf_ARX_pred (Uv_model, A_M1, B_M1, dt_M1);
end

% Model accuracy
if toggle_ascale == 1;	% if auto-scaling is used, convert to Eng units
	Yt_M1 = Yt_model_M1 .* Yt_std + Yt_m;
	Yv_M1 = Yv_model_M1 .* Yt_std + Yt_m;
	[RMSE_Yt_M1] = stat_RMSE(Yt_val(max(nB_in)+1:end), Yt_M1(max(nB_in)+1:end));
	[RMSE_Yv_M1] = stat_RMSE(Yv_val(max(nB_in)+1:end), Yv_M1(max(nB_in)+1:end));	
else
	Yt_M1 = Yt_model_M1;
	Yv_M1 = Yv_model_M1;
	[RMSE_Yt_M1] = stat_RMSE(Yt_model(max(nB_in)+1:end), Yt_model_M1(max(nB_in)+1:end));
	[RMSE_Yv_M1] = stat_RMSE(Yv_model(max(nB_in)+1:end), Yv_model_M1(max(nB_in)+1:end));
end

% M2: CON1: Constrained (quadprog, specify: K_sign)
% set up constraints [B A]
K_sign_M2 = [1 1 -1 0];				% specify sign direction vector
min_phase_M2 = zeros(1,n_var+1);	% 0 = not minimum phase
dt_M2 = zeros(1,n_var);				% specify a priori dead-time
nB_M2 = nB_in - dt_M2;

nZ_M2 = [nB_M2 nA_in];
[con_A_M2 con_b_M2] = CMI_con(nZ_M2, K_sign_M2, min_phase_M2);

% Model Identification and Output Prediction
if toggle_inc == 1;		% incremental model
	[B_M2 A_M2] = QP_ARX_MISO (Ut_model_i, Yt_model_i, nB_M2, dt_M2, nA_in, con_A_M2, con_b_M2);	% parameter estimation
	[Yt_model_M2] = tf_ARX_pred_inc (Ut_model_i, A_M2, B_M2, dt_M2, Yt_model);
	[Yv_model_M2] = tf_ARX_pred_inc (Uv_model_i, A_M2, B_M2, dt_M2, Yv_model);
elseif toggle_inc == 0;	% absolute model
	[B_M2 A_M2] = QP_ARX_MISO (Ut_model, Yt_model, nB_M2, dt_M2, con_A_M2, con_b_M2);	% parameter estimation
	[Yt_model_M2] = tf_ARX_pred (Ut_model, A_M2, B_M2, dt_M2);
	[Yv_model_M2] = tf_ARX_pred (Uv_model, A_M2, B_M2, dt_M2);
end

% Model accuracy
if toggle_ascale == 1;	% if auto-scaling is used, convert to Eng units
	Yt_M2 = Yt_model_M2 .* Yt_std + Yt_m;
	Yv_M2 = Yv_model_M2 .* Yt_std + Yt_m;
	[RMSE_Yt_M2] = stat_RMSE(Yt_val(max(nB_in)+1:end), Yt_M2(max(nB_in)+1:end));
	[RMSE_Yv_M2] = stat_RMSE(Yv_val(max(nB_in)+1:end), Yv_M2(max(nB_in)+1:end));	
else
	Yt_M2 = Yt_model_M2;
	Yv_M2 = Yv_model_M2;
	[RMSE_Yt_M2] = stat_RMSE(Yt_model(max(nB_in)+1:end), Yt_model_M2(max(nB_in)+1:end));
	[RMSE_Yv_M2] = stat_RMSE(Yv_model(max(nB_in)+1:end), Yv_model_M2(max(nB_in)+1:end));
end

% M3: CON2: Constrained (quadprog, specify: K_sign min_phase)
% set up constraints
K_sign_M3 = [1 1 -1 0];			% specify sign direction vector
min_phase_M3 = ones(1,n_var+1);	% 0 = not minimum phase
dt_M3 = zeros(1,n_var);			% specify a priori dead-time
nB_M3 = nB_in - dt_M3;

nZ_M3 = [nB_M3 nA_in];
[con_A_M3 con_b_M3] = CMI_con(nZ_M3, K_sign_M3, min_phase_M3);

%%%%%%%%%% CONTINUE FROM HERE!!!!

% Model Identification and Output Prediction
if toggle_inc == 1;		% incremental model
	[B_M3 A_M3] = QP_ARX_MISO (Ut_model_i, Yt_model_i, nB_M3, dt_M3,  nA_in, con_A_M3, con_b_M3);	% parameter estimation
	[Yt_model_M3] = tf_ARX_pred_inc (Ut_model_i, A_M3, B_M3, dt_M3, Yt_model);
	[Yv_model_M3] = tf_ARX_pred_inc (Uv_model_i, A_M3, B_M3, dt_M3, Yv_model);
elseif toggle_inc == 0;	% absolute model
	[B_M3] = QP_ARX_MISO (Ut_model, Yt_model, nB_M3, dt_M3, nA_in, con_A_M3, con_b_M3);	% parameter estimation
	[Yt_model_M3] = tf_ARX_pred (Ut_model, A_M3, B_M3, dt_M3);
	[Yv_model_M3] = tf_ARX_pred (Uv_model, A_M3, B_M3, dt_M3);

end

% Model accuracy
if toggle_ascale == 1;	% if auto-scaling is used, convert to Eng units
	Yt_M3 = Yt_model_M3 .* Yt_std + Yt_m;
	Yv_M3 = Yv_model_M3 .* Yt_std + Yt_m;
	[RMSE_Yt_M3] = stat_RMSE(Yt_val(max(nB_in)+1:end), Yt_M3(max(nB_in)+1:end));
	[RMSE_Yv_M3] = stat_RMSE(Yv_val(max(nB_in)+1:end), Yv_M3(max(nB_in)+1:end));	
else
	Yt_M3 = Yt_model_M3;
	Yv_M3 = Yv_model_M3;
	[RMSE_Yt_M3] = stat_RMSE(Yt_model(max(nB_in)+1:end), Yt_model_M3(max(nB_in)+1:end));
	[RMSE_Yv_M3] = stat_RMSE(Yv_model(max(nB_in)+1:end), Yv_model_M3(max(nB_in)+1:end));
end