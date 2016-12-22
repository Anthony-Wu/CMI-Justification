% CMI_FIR_pt4.m
% Justification for constrained model identification - FIR models
% part 4: Model identification

n_coeff = [100];
	for i_var = 1:n_var;
		B_true{i_var} = [B{i_var}; zeros(n_coeff - size(B{i_var},1),1)];
%CHECK THIS BIT		
		if toggle_ascale == 1;
			if toggle_filter == 1;
				B_true_mv_as{i_var} = [B_mv_as{i_var}; zeros(n_coeff - size(B_mv_as{i_var},1),1)];				
			else
				B_true_as{i_var} = [B_as{i_var}; zeros(n_coeff - size(B_as{i_var},1),1)];
			end
		else
			B_true_mv{i_var} = [B_mv{i_var}; zeros(n_coeff - size(B_mv{i_var},1),1)];	
		end
			
		
		
	end
dt = [0];

nB_total = n_coeff + dt; %total number of FIR coefficients fed to the model

toggle_inc = 1; 	% 0 = use absolute model, 1 = use incremental model
if toggle_inc == 1;
	Ut_model_i = abs2inc(Ut_model);
	Yt_model_i = abs2inc(Yt_model);
	Uv_model_i = abs2inc(Uv_model);
	Yv_model_i = abs2inc(Yv_model);	
end

if toggle_ascale == 1;	% Y_val is the output used for validation in Eng units (if moving average is used, it can skew the results)
	Yt_val = Yt_model .* Yt_std + Yt_m;
	Yv_val = Yv_model .* Yt_std + Yt_m;
end		
		
% M1: UNC1: Unconstrained using WaterMV
	% Parameter estimation carried out in WaterMV, coefficients c/p here
	B_M1{1} = [0.1 0 0 0 0 zeros(1,95)]'; % export and c/p into here
	B_M1{2} = [0.1 0 0 0 0 zeros(1,95)]';
	B_M1{3} = [0.1 0 0 0 0 zeros(1,95)]';
	
	dt_M1 = zeros(1,n_var);			% specify a priori dead-time
	
	% Output Prediction
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
	[RMSE_Yt_model_M1] = stat_RMSE (Yt_model(n_coeff+1:end), Yt_model_M1(n_coeff+1:end));
	[RMSE_Yv_model_M1] = stat_RMSE (Yv_model(n_coeff+1:end), Yv_model_M1(n_coeff+1:end));
	
	if toggle_ascale == 1;	% if auto-scaling is used, convert to Eng units
		Yt_M1 = Yt_model_M1 .* Yt_std + Yt_m;
		Yv_M1 = Yv_model_M1 .* Yt_std + Yt_m;
		[RMSE_Yt_M1] = stat_RMSE (Yt_val(n_coeff+1:end), Yt_M1(n_coeff+1:end));
		[RMSE_Yv_M1] = stat_RMSE (Yv_val(n_coeff+1:end), Yv_M1(n_coeff+1:end));	
	end

% M2: CON1: Constrained (quadprog, specify: K_sign)
	% set up constraints
	K_sign_M2 = [1 1 -1];				% specify sign direction vector
	min_phase_M2 = zeros(1,n_var);	% 0 = not minimum phase
	dt_M2 = zeros(1,n_var);			% specify a priori dead-time
	nB_M2 = nB_total - dt_M2;

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
	[RMSE_Yt_model_M2] = stat_RMSE (Yt_model(n_coeff+1:end), Yt_model_M2(n_coeff+1:end));
	[RMSE_Yv_model_M2] = stat_RMSE (Yv_model(n_coeff+1:end), Yv_model_M2(n_coeff+1:end));
	
	if toggle_ascale == 1;	% if auto-scaling is used, convert to Eng units
		Yt_M2 = Yt_model_M2 .* Yt_std + Yt_m;
		Yv_M2 = Yv_model_M2 .* Yt_std + Yt_m;
		[RMSE_Yt_M2] = stat_RMSE (Yt_val(n_coeff+1:end), Yt_M2(n_coeff+1:end));
		[RMSE_Yv_M2] = stat_RMSE (Yv_val(n_coeff+1:end), Yv_M2(n_coeff+1:end));	
	end

% M3: CON2: Constrained (quadprog, specify: K_sign min_phase)
	% set up constraints
	K_sign_M3 = [1 1 -1];				% specify sign direction vector
	min_phase_M3 = ones(1,n_var);	% 0 = not minimum phase
	dt_M3 = zeros(1,n_var);			% specify a priori dead-time
	nB_M3 = nB_total - dt_M3;

	nZ_M3 = nB_M3;
	[con_A_M3 con_b_M3] = CMI_con(nZ_M3, K_sign_M3, min_phase_M3);

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
	[RMSE_Yt_model_M3] = stat_RMSE (Yt_model(n_coeff+1:end), Yt_model_M3(n_coeff+1:end));
	[RMSE_Yv_model_M3] = stat_RMSE (Yv_model(n_coeff+1:end), Yv_model_M3(n_coeff+1:end));
	
	if toggle_ascale == 1;	% if auto-scaling is used, convert to Eng units
		Yt_M3 = Yt_model_M3 .* Yt_std + Yt_m;
		Yv_M3 = Yv_model_M3 .* Yt_std + Yt_m;
		[RMSE_Yt_M3] = stat_RMSE (Yt_val(n_coeff+1:end), Yt_M3(n_coeff+1:end));
		[RMSE_Yv_M3] = stat_RMSE (Yv_val(n_coeff+1:end), Yv_M3(n_coeff+1:end));	
	end	
	