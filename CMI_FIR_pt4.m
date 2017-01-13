% CMI_FIR_pt4.m
% Justification for constrained model identification - FIR models
% part 4: Model identification

nB_in = [60 60 60];	% total number of FIR coefficients (incl. dead time) for each input
dt_in = [0 0 0];	% dead time for each input

% align B coefficient vectors to the same size (for later comparison)
for i_var = 1:n_var;
	B_true{i_var} = [B{i_var}; zeros(nB_in(i_var) - length(B{i_var}),1)];	
	if  toggle_ascale == 1; % autoscaled
		B_true_mv_as{i_var} = filter(ones(1,mov_avg)/mov_avg, 1, B_true{i_var} * (Ut_std(i_var)/Yt_std)); 
	elseif toggle_ascale == 0; % not autoscaled
		B_true_mv{i_var} = filter(ones(1,mov_avg)/mov_avg, 1, B_true{i_var}); 
	end	
end
	
%%%	
	% M1: UNC1: Unconstrained using WaterMV
%%%
dt_M1 = dt_in;			% specify a priori dead-time

% Parameter estimation carried out in WaterMV, coefficients c/p here
B_M1nz{1} = []';	% c/p coefficients from WaterMV as appropriate
B_M1nz{2} = []';
B_M1nz{3} = []';

% align B coefficient vectors to the same size (for later comparison)
for i_var = 1:n_var
	B_M1{i_var} = [zeros(1,dt_M1(i_var)) [B_M1nz{i_var}]' zeros(1,nB_in(i_var) - length([zeros(1,dt_M1(i_var)) [B_M1nz{i_var}]']))]';
end

% Output Prediction
if toggle_inc == 1;		% incremental model
	[Yt_model_M1] = tf_FIR_pred_inc (Ut_model_i, B_M1, dt_M1, Yt_model);
	[Yv_model_M1] = tf_FIR_pred_inc (Uv_model_i, B_M1, dt_M1, Yv_model);
elseif toggle_inc == 0;	% absolute model
	[Yt_model_M1] = tf_FIR_pred (Ut_model, B_M1, dt_M1);
	[Yv_model_M1] = tf_FIR_pred (Uv_model, B_M1, dt_M1);
end

if toggle_ascale == 1;	% if auto-scaling is used, convert to Eng units
	Yt_M1 = Yt_model_M1 .* Yt_std + Yt_m;
	Yv_M1 = Yv_model_M1 .* Yt_std + Yt_m;
else	% if auto-scaling not used, then model is already in Eng units
	Yt_M1 = Yt_model_M1;
	Yv_M1 = Yv_model_M1;
end
[RMSE_Yt_M1] = stat_RMSE(Yt_val(max(nB_in)+1:end), Yt_M1(max(nB_in)+1:end));
[RMSE_Yv_M1] = stat_RMSE(Yv_val(max(nB_in)+1:end), Yv_M1(max(nB_in)+1:end));

%%%	
	% M2: CON1: Constrained (quadprog, specify: K_sign)
%%%
% set up constraints
K_sign_M2 = [1 1 -1];			% specify sign direction vector
min_phase_M2 = zeros(1,n_var);	% 0 = not minimum phase
dt_M2 = dt_in;					% specify a priori dead-time
nB_M2 = nB_in - dt_M2;

nZ_M2 = nB_M2;
[con_A_M2 con_b_M2] = CMI_con(nZ_M2, K_sign_M2, min_phase_M2);	% set up inequality constraint matrices

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

if toggle_ascale == 1;	% if auto-scaling is used, convert to Eng units
	Yt_M2 = Yt_model_M2 .* Yt_std + Yt_m;
	Yv_M2 = Yv_model_M2 .* Yt_std + Yt_m;
else	% if auto-scaling not used, then model is already in Eng units
	Yt_M2 = Yt_model_M2;
	Yv_M2 = Yv_model_M2;
end
[RMSE_Yt_M2] = stat_RMSE(Yt_val(max(nB_in)+1:end), Yt_M2(max(nB_in)+1:end));
[RMSE_Yv_M2] = stat_RMSE(Yv_val(max(nB_in)+1:end), Yv_M2(max(nB_in)+1:end));

%%%
	% M3: CON2: Constrained (quadprog, specify: K_sign min_phase)
%%%
% set up constraints
K_sign_M3 = [1 1 -1];			% specify sign direction vector
min_phase_M3 = ones(1,n_var);	% 0 = not minimum phase
dt_M3 = dt_in;					% specify a priori dead-time
nB_M3 = nB_in - dt_M3;

nZ_M3 = nB_M3;
[con_A_M3 con_b_M3] = CMI_con(nZ_M3, K_sign_M3, min_phase_M3);	% set up inequality constraint matrices

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

if toggle_ascale == 1;	% if auto-scaling is used, convert to Eng units
	Yt_M3 = Yt_model_M3 .* Yt_std + Yt_m;
	Yv_M3 = Yv_model_M3 .* Yt_std + Yt_m;	
else	% if auto-scaling not used, then model is already in Eng units
	Yt_M3 = Yt_model_M3;
	Yv_M3 = Yv_model_M3;
end
[RMSE_Yt_M3] = stat_RMSE(Yt_val(max(nB_in)+1:end), Yt_M3(max(nB_in)+1:end));
[RMSE_Yv_M3] = stat_RMSE(Yv_val(max(nB_in)+1:end), Yv_M3(max(nB_in)+1:end));

if toggle_inc == 1;
	Yt_i_M1 = abs2inc(Yt_M1);
	Yt_i_M2 = abs2inc(Yt_M2);
	Yt_i_M3 = abs2inc(Yt_M3);
	Yt_i_val = abs2inc(Yt_val);
	
	Yv_i_M1 = abs2inc(Yv_M1);
	Yv_i_M2 = abs2inc(Yv_M2);
	Yv_i_M3 = abs2inc(Yv_M3);
	Yv_i_val = abs2inc(Yv_val);
	
	[RMSE_Yt_i_M1] = stat_RMSE(Yt_i_val(max(nB_in)+1:end), Yt_i_M1(max(nB_in)+1:end));
	[RMSE_Yv_i_M1] = stat_RMSE(Yv_i_val(max(nB_in)+1:end), Yv_i_M1(max(nB_in)+1:end));	
	[RMSE_Yt_i_M2] = stat_RMSE(Yt_i_val(max(nB_in)+1:end), Yt_i_M2(max(nB_in)+1:end));
	[RMSE_Yv_i_M2] = stat_RMSE(Yv_i_val(max(nB_in)+1:end), Yv_i_M2(max(nB_in)+1:end));	
	[RMSE_Yt_i_M3] = stat_RMSE(Yt_i_val(max(nB_in)+1:end), Yt_i_M3(max(nB_in)+1:end));
	[RMSE_Yv_i_M3] = stat_RMSE(Yv_i_val(max(nB_in)+1:end), Yv_i_M3(max(nB_in)+1:end));
end