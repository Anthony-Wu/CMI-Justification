% CMI_FIR_pt4.m
% Justification for constrained model identification - FIR models
% part 4: Model identification

nB_in = [100 100 100];	% total number of FIR coefficients (incl. dead time) for each input
dt_in = [0 0 0];	% dead time for each input
nA_in = [5];	% number of ARX coefficients

% align B coefficient vectors to the same size (for later comparison)
for i_var = 1:n_var;
	B_true{i_var} = [B{i_var}; zeros(nB_in(i_var) - length(B{i_var}),1)];	
	if  toggle_ascale == 1; % autoscaled
		B_true_mv_as{i_var} = filter(ones(1,mov_avg)/mov_avg, 1, B_true{i_var} * (Ut_std(i_var)/Yt_std)); 
	elseif toggle_ascale == 0; % not autoscaled
		B_true_mv{i_var} = filter(ones(1,mov_avg)/mov_avg, 1, B_true{i_var}); 
	end	
end

A_true{1} = [A{1}; zeros(nA_in(1) - length(A{1}),1)];

%%%
	% M1: UNC1: Unconstrained using WaterMV
%%%
dt_M1 = dt_in;			% specify a priori dead-time

% Parameter estimation carried out in WaterMV, coefficients c/p here
B_M1nz{1} = [-0.011987588	-0.0057527414	0.007253445	0.011857908	0.0023792558	0.007016613	0.003936901	0.0034470006	-0.0016745538	0.00289686	6.901158E-4	0.0028059082	-0.0057364143	0.013294969	-0.006560619	0.0068609	-0.0010447135	-0.0013204764	-0.011152418	-3.559586E-5	0.002687548	0.005546793	-0.005689342	-1.7746993E-4	9.920894E-4	0.011503805	0.010282977	0.017424198	0.02686762	0.017682565	0.018068839	0.029436272	0.03258679	0.030305712	0.03852943	0.054252107	0.037873425	0.06303217	0.068449676	0.06721381	0.055926893	0.05877128	0.052766625	0.055816717	0.04993967	0.06298661	0.05827245	0.051200964	0.05962377	0.05477502	0.05784205	0.053453524	0.056549508	0.03638006	0.04268573	0.030603511	0.035144072	0.036302913	0.045949753	0.02184338	0.03342671	0.02059183	0.0364352	0.029513989	0.031333674	0.024428977	0.02876124	0.01941057	0.026958166	0.024837151	0.018126165	0.010359048	0.024690477	0.01802576	0.019912051	0.020023234	0.021872934	0.005603928	0.018073022	0.013392293	0.025120705	0.026186347	0.029918116	0.020677412	0.025449554	0.022899633	0.013838954	0.017127136	0.017586462	0.011259005	0.02050263	0.021764575	0.007209273	0.026318552	0.014904868	0.016579378	0.014635498	0.004848116	0.0082448125	0.014560787]';	% c/p coefficients from WaterMV as appropriate
B_M1nz{2} = [-0.0017205264	-0.004123831	-8.5481803E-4	0.0070157605	0.008030039	0.0013454156	6.7636545E-4	0.015164872	0.01314317	0.002864274	0.02308978	0.015815217	0.014329067	0.014632141	0.02371124	0.01841551	0.031480134	0.019535197	0.016028501	0.024796568	0.031542845	0.023779638	0.02418676	0.013444123	0.009227057	0.027098669	0.014924913	0.022896292	0.011579713	0.021485105	0.026049253	0.01865381	2.7321477E-4	0.005862249	0.008518895	-1.20851866E-4	0.0018446365	0.007337647	0.007869093	0.0054891147	0.012093862	0.007727409	0.0022909204	0.02031871	0.0050762817	0.0069372538	0.0028981455	0.002793185	0.0034153825	0.0038934734	-0.0018185938	-0.0014117463	0.008432572	0.0077188914	0.005740926	-0.009829618	0.008017758	-0.0041511613	0.003243161	0.0031283307	-0.008337792	-5.9995404E-4	-0.012281572	0.0026312356	-0.0056024473	-0.0032375285	-0.01052555	0.0062576826	-0.009172198	-0.0024251817	-3.2086924E-4	-0.0017519693	0.009421743	0.0033800995	0.0037626163	0.0069932393	0.007585263	-0.0017850754	0.0043319417	-0.010062177	0.0015615138	-0.008932238	-0.0025059292	-3.1803868E-4	-0.0011864385	0.011637554	-0.008856302	-0.0016866146	-0.0065541184	-0.0021252197	-0.0036870546	0.0027007414	-0.016701236	-0.0036788261	-0.004190152	-0.01076423	0.0014028964	-0.004970243	-0.003971945	0.007063645]';
B_M1nz{3} = [-0.0033694347	-0.005107778	-0.0018675411	-0.008170664	-0.0019260517	-0.013768894	0.010463336	-0.008283604	-0.0027691973	-0.0020794538	-0.009448676	-0.004575997	-0.002172275	-0.009161916	-0.009467038	-0.007644529	-0.02219704	-0.010629399	-0.013429573	-0.0036405646	-0.0061624073	-0.0096675735	-0.016159989	-0.014042245	-0.022514142	-0.01562619	-0.012420389	-0.036834747	-0.026422659	-0.02321607	-0.022381706	-0.031031212	-0.035229728	-0.0235681	-0.016340178	-0.024721917	-0.023409229	-0.031020079	-0.034269728	-0.030713363	-0.03346466	-0.03348044	-0.026210634	-0.024094017	-0.021263357	-0.027441662	-0.030208394	-0.019396896	-0.024189416	-0.028918993	-0.026369749	-0.010154233	-0.019164495	-0.014041859	-0.017412733	-0.015535995	-0.023765983	-0.02057881	-0.027597418	-0.020526987	-0.009123457	-0.016843155	-0.015246755	-0.019801334	-0.0072647617	-0.020142369	-0.0068951966	-0.0139184585	-0.009042887	-0.012938097	-0.023289563	-0.011692046	-0.024703044	-0.0062562847	-0.019427326	-0.0032177547	-0.020074936	-0.0064843185	-0.0030704858	-0.0050291363	-0.012710491	-0.0017562486	-0.013808743	-0.015901262	-0.014233466	-0.0055513154	-0.0074711023	-0.01202635	-0.0016139116	-0.0054719397	-0.004121956	-0.0045022685	0.008235204	-0.012411097	-0.020993259	-0.005142883	-0.0028117846	-0.0014762131	-0.011289429	-0.0064363074]';

A_M1nz{1} = -[-0.07345935	-0.08155435	-0.008728936	-0.052472826	-0.027449515]';

% align B coefficient vectors to the same size (for later comparison)
for i_var = 1:n_var
	B_M1{i_var} = [zeros(1,dt_M1(i_var)) [B_M1nz{i_var}]' zeros(1,nB_in(i_var) - length([zeros(1,dt_M1(i_var)) [B_M1nz{i_var}]']))]';
end
A_M1{1} = [[A_M1nz{1}]' zeros(1,nA_in(1) - length(A_M1nz{1}))]';

% Output Prediction
if toggle_inc == 1;		% incremental model
	[Yt_model_M1] = tf_ARX_pred_inc (Ut_model_i, A_M1, B_M1, dt_M1, Yt_model);
	[Yv_model_M1] = tf_ARX_pred_inc (Uv_model_i, A_M1, B_M1, dt_M1, Yv_model);	
elseif toggle_inc == 0;	% absolute model
	[Yt_model_M1] = tf_ARX_pred (Ut_model, A_M1, B_M1, dt_M1);
	[Yv_model_M1] = tf_ARX_pred (Uv_model, A_M1, B_M1, dt_M1);
end

if toggle_ascale == 1;	% if auto-scaling is used, convert to Eng units
	Yt_M1 = Yt_model_M1 .* Yt_std + Yt_m;
	Yv_M1 = Yv_model_M1 .* Yt_std + Yt_m;
else
	Yt_M1 = Yt_model_M1;
	Yv_M1 = Yv_model_M1;
end
[RMSE_Yt_M1] = stat_RMSE(Yt_val(max(nB_in)+1:end), Yt_M1(max(nB_in)+1:end));
[RMSE_Yv_M1] = stat_RMSE(Yv_val(max(nB_in)+1:end), Yv_M1(max(nB_in)+1:end));

%%%
	% M2: CON1: Constrained (quadprog, specify: K_sign)
%%%
% set up constraints [B A]
dt_M2 = dt_in;			% specify a priori dead-time
nB_M2 = nB_in - dt_M2;
nZ_M2 = [nB_M2 nA_in];
K_sign_M2 = [1 1 -1 -1];	% specify sign direction vector
K_bound_M2 = [[zeros(1,n_var)] -0.5];
min_phase_M2 = zeros(1,n_var+1);	% 0 = not minimum phase

[con_A_M2 con_b_M2] = CMI_con2(nZ_M2, K_sign_M2, K_bound_M2, min_phase_M2);

% Model Identification and Output Prediction
if toggle_inc == 1;		% incremental model
	[B_M2 A_M2] = QP_ARX_MISO (Ut_model_i, Yt_model_i, nB_M2, dt_M2, nA_in, con_A_M2, con_b_M2);	% parameter estimation
	[Yt_model_M2] = tf_ARX_pred_inc (Ut_model_i, A_M2, B_M2, dt_M2, Yt_model);
	[Yv_model_M2] = tf_ARX_pred_inc (Uv_model_i, A_M2, B_M2, dt_M2, Yv_model);	
elseif toggle_inc == 0;	% absolute model
	[B_M2 A_M2] = QP_ARX_MISO (Ut_model, Yt_model, nB_M2, dt_M2, nA_in, con_A_M2, con_b_M2);	% parameter estimation
	[Yt_model_M2] = tf_ARX_pred (Ut_model, A_M2, B_M2, dt_M2);
	[Yv_model_M2] = tf_ARX_pred (Uv_model, A_M2, B_M2, dt_M2);
end

if toggle_ascale == 1;	% if auto-scaling is used, convert to Eng units
	Yt_M2 = Yt_model_M2 .* Yt_std + Yt_m;
	Yv_M2 = Yv_model_M2 .* Yt_std + Yt_m;
else
	Yt_M2 = Yt_model_M2;
	Yv_M2 = Yv_model_M2;
end
[RMSE_Yt_M2] = stat_RMSE(Yt_val(max(nB_in)+1:end), Yt_M2(max(nB_in)+1:end));
[RMSE_Yv_M2] = stat_RMSE(Yv_val(max(nB_in)+1:end), Yv_M2(max(nB_in)+1:end));
	
% M3: CON2: Constrained (quadprog, specify: K_sign min_phase)
% set up constraints
dt_M3 = dt_in;			% specify a priori dead-time
nB_M3 = nB_in - dt_M3;
nZ_M3 = [nB_M3 nA_in];
K_sign_M3 = [1 1 -1 -1];	% specify sign direction vector
K_bound_M3 = [[zeros(1,n_var)] -0.5];
min_phase_M3 = ones(1,n_var+1);	% 0 = not minimum phase

[con_A_M3 con_b_M3] = CMI_con2(nZ_M3, K_sign_M3, K_bound_M3, min_phase_M3);

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

if toggle_ascale == 1;	% if auto-scaling is used, convert to Eng units
	Yt_M3 = Yt_model_M3 .* Yt_std + Yt_m;
	Yv_M3 = Yv_model_M3 .* Yt_std + Yt_m;
else
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