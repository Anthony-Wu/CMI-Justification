% CMI_FIR_pt4.m
% Justification for constrained model identification - FIR models
% part 4: Model identification

nB_in = [30 30];
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
B_M1nz{1} = [4.0602956E-8	9.544616E-8	1.9004692E-7	1.9277677E-7	2.545831E-7	3.0050276E-7	2.715955E-7	3.0153683E-7	3.2304686E-7	0.0845655	0.08004317	0.073393695	0.07000753	0.067113966	0.06469593	0.062471442	0.060380884	0.058382478	0.05646055	-0.005987558	-0.0045387954	-0.0015063493	-7.549232E-4	-3.0142855E-4	-1.355284E-4	-5.698552E-5	-2.468315E-5	-1.0351333E-5	-4.254277E-6	-1.6373056E-6]';	% c/p coefficients from WaterMV as appropriate
B_M1nz{2} = [-6.690077E-8	-5.6318093E-8	-6.324389E-9	7.3003754E-9	0.061370354	0.054260135	0.046175323	0.041069977	0.03672721	0.03307412	0.029851656	0.026980093	0.024398932	0.022071294	-0.0026085938	-0.0018940879	-6.3970604E-4	-3.1737742E-4	-1.2748322E-4	-5.725073E-5	-2.419282E-5	-1.054802E-5	-4.51672E-6	-2.0088546E-6	-9.5111284E-7	-4.3632437E-7	-2.6927898E-7	-1.5663855E-7	-9.331918E-8	-6.301523E-8]';
%B_M1nz{3} = [-1.7058079E-6	-3.5796213E-6	-1.9182355E-6	-2.2785678E-6	0.10811276	0.09417671	0.07973014	0.07068096	0.06311351	0.056792792	0.05124112	0.046306983	0.041870568	0.037874434	-0.005506495	-0.0036425316	-0.001281566	-6.225813E-4	-2.5240966E-4	-1.15536444E-4	-5.194282E-5	-2.3937884E-5	-9.318589E-6	-5.3212316E-6	-3.3944075E-6	-3.0280764E-6	-1.74369E-6	-1.4976952E-6	6.255219E-7	2.1069493E-6]';

A_M1nz{1} = [-0.077233404	-0.14346129]';
B_M1{1} = [[B_M1nz{1}]' zeros(1,nB_in(1) - length(B_M1nz{1}))]';
B_M1{2} = [[B_M1nz{2}]' zeros(1,nB_in(2) - length(B_M1nz{2}))]';
%B_M1{3} = [[B_M1nz{3}]' zeros(1,nB_in(3) - length(B_M1nz{3}))]';
A_M1{1} = [[A_M1nz{1}]' zeros(1,nA_in(1) - length(A_M1nz{1}))]';

dt_M1 = zeros(1,n_var);			% specify a priori dead-time

% Output Prediction
% WaterMV's notation is slightly different and needs a shift in coefficeints by 1 to translate over
for i_var = 1:n_var
	B_M1{i_var} = [0; [B_M1{i_var}(1:end-1)]];
end
[Yt_model_M1] = tf_ARX_pred_inc (Ut_model_i, A_M1, B_M1, dt_M1, Yt_model);
[Yv_model_M1] = tf_ARX_pred_inc (Uv_model_i, A_M1, B_M1, dt_M1, Yv_model);

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
K_sign_M2 = [1 1 -1];				% specify sign direction vector
min_phase_M2 = zeros(1,n_var+1);	% 0 = not minimum phase
% dt_M2 = zeros(1,n_var);				% specify a priori dead-time
dt_M2 = [10 5];				% specify a priori dead-time
nB_M2 = nB_in - dt_M2;

nZ_M2 = [nB_M2 nA_in];
[con_A_M2 con_b_M2] = CMI_con(nZ_M2, K_sign_M2, min_phase_M2);

% Model Identification and Output Prediction
[B_M2 A_M2] = QP_ARX_MISO (Ut_model_i, Yt_model_i, nB_M2, dt_M2, nA_in, con_A_M2, con_b_M2);	% parameter estimation
[Yt_model_M2] = tf_ARX_pred_inc (Ut_model_i, A_M2, B_M2, zeros(1,n_var), Yt_model);
[Yv_model_M2] = tf_ARX_pred_inc (Uv_model_i, A_M2, B_M2, zeros(1,n_var), Yv_model);

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
min_phase_M3 = ones(1,n_var+1);	% 0 = not minimum phase
% dt_M3 = zeros(1,n_var);			% specify a priori dead-time
dt_M3 = [10 5];
nB_M3 = nB_in - dt_M3;

nZ_M3 = [nB_M3 nA_in];
[con_A_M3 con_b_M3] = CMI_con(nZ_M3, K_sign_M3, min_phase_M3);

% Model Identification and Output Prediction
[B_M3 A_M3] = QP_ARX_MISO (Ut_model_i, Yt_model_i, nB_M3, dt_M3,  nA_in, con_A_M3, con_b_M3);	% parameter estimation
[Yt_model_M3] = tf_ARX_pred_inc (Ut_model_i, A_M3, B_M3, zeros(1,n_var), Yt_model);
[Yv_model_M3] = tf_ARX_pred_inc (Uv_model_i, A_M3, B_M3, zeros(1,n_var), Yv_model);

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