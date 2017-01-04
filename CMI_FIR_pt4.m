% CMI_FIR_pt4.m
% Justification for constrained model identification - FIR models
% part 4: Model identification

nB_in = [100 100 100];
for i_var = 1:n_var;
	B_true{i_var} = [B{i_var}; zeros(nB_in(i_var) - length(B{i_var}),1)];
	if toggle_ascale == 1;
		if 	toggle_filter == 1; 	% moving average and no auto-scaling used
			B_true_mv_as{i_var} = [B_mv_as{i_var}; zeros(nB_in(i_var) - length(B{i_var}),1)];				
		elseif toggle_filter == 0; 	% auto-scaling used, but no moving average
			B_true_as{i_var} = [B_as{i_var}; zeros(nB_in(i_var) - length(B{i_var}),1)];
		end
	elseif 	toggle_filter == 1;		% moving average used, but no auto-scaling
		B_true_mv{i_var} = [B_mv{i_var}; zeros(nB_in(i_var) - length(B{i_var}),1)];	
	end
end

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
B_M1nz{1} = [0.1 0 0 0 0]';	% c/p coefficients from WaterMV as appropriate
B_M1nz{2} = [0.1 0 0 0 0]';
B_M1nz{3} = [0.1 0 0 0 0]';
B_M1{1} = [[B_M1nz{1}]' zeros(1,nB_in(1) - length(B_M1nz{1}))]';
B_M1{2} = [[B_M1nz{2}]' zeros(1,nB_in(2) - length(B_M1nz{2}))]';
B_M1{3} = [[B_M1nz{3}]' zeros(1,nB_in(3) - length(B_M1nz{3}))]';

dt_M1 = zeros(1,n_var);			% specify a priori dead-time

% Model Identification and Output Prediction
if toggle_inc == 1;		% incremental model
	% WaterMV's notation is slightly different and needs a shift in coefficeints by 1 to translate over
	for i_var = 1:n_var
		B_M1{i_var} = [0; [B_M1{i_var}(1:end-1)]];
	end
	[Yt_model_M1] = tf_FIR_pred_inc (Ut_model_i, B_M1, dt_M1, Yt_model);
	[Yv_model_M1] = tf_FIR_pred_inc (Uv_model_i, B_M1, dt_M1, Yv_model);
elseif toggle_inc == 0;	% absolute model
	[Yt_model_M1] = tf_FIR_pred (Ut_model, B_M1, dt_M1);
	[Yv_model_M1] = tf_FIR_pred (Uv_model, B_M1, dt_M1);
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
% set up constraints
K_sign_M2 = [1 1 -1];			% specify sign direction vector
min_phase_M2 = zeros(1,n_var);	% 0 = not minimum phase
dt_M2 = zeros(1,n_var);			% specify a priori dead-time
nB_M2 = nB_in - dt_M2;

nZ_M2 = nB_M2;
[con_A_M2 con_b_M2] = CMI_con(nZ_M2, K_sign_M2, min_phase_M2);

% Model Identification and Output Prediction
if toggle_inc == 1;		% incremental model
	[B_M2] = QP_FIR_MISO (Ut_model_i, Yt_model_i, nB_M2, dt_M2, con_A_M2, con_b_M2);	% parameter estimation
	[Yt_model_M2] = tf_FIR_pred_inc (Ut_model_i, B_M2, dt_M2, Yt_model);
	[Yv_model_M2] = tf_FIR_pred_inc (Uv_model_i, B_M2, dt_M2, Yv_model);
elseif toggle_inc == 0;	% absolute model
	[B_M2] = QP_FIR_MISO (Ut_model, Yt_model, nB_M2, dt_M2, con_A_M2, con_b_M2);	% parameter estimation
	[Yt_model_M2] = tf_FIR_pred (Ut_model, B_M2, dt_M2);
	[Yv_model_M2] = tf_FIR_pred (Uv_model, B_M2, dt_M2);
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
K_sign_M3 = [1 1 -1];			% specify sign direction vector
min_phase_M3 = ones(1,n_var);	% 0 = not minimum phase
dt_M3 = zeros(1,n_var);			% specify a priori dead-time
nB_M3 = nB_in - dt_M3;

nZ_M3 = nB_M3;
[con_A_M3 con_b_M3] = CMI_con(nZ_M3, K_sign_M3, min_phase_M3);

% Model Identification and Output Prediction
if toggle_inc == 1;		% incremental model
	[B_M3] = QP_FIR_MISO (Ut_model_i, Yt_model_i, nB_M3, dt_M3, con_A_M3, con_b_M3);	% parameter estimation
	[Yt_model_M3] = tf_FIR_pred_inc (Ut_model_i, B_M3, dt_M3, Yt_model);
	[Yv_model_M3] = tf_FIR_pred_inc (Uv_model_i, B_M3, dt_M3, Yv_model);
elseif toggle_inc == 0;	% absolute model
	[B_M3] = QP_FIR_MISO (Ut_model, Yt_model, nB_M3, dt_M3, con_A_M3, con_b_M3);	% parameter estimation
	[Yt_model_M3] = tf_FIR_pred (Ut_model, B_M3, dt_M3);
	[Yv_model_M3] = tf_FIR_pred (Uv_model, B_M3, dt_M3);

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