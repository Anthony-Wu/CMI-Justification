% CMI_FIR_pt4.m
% Justification for constrained model identification - FIR models
% part 4: Model identification

nB_in = [100];	% total number of FIR coefficients (incl. dead time) for each input
dt_in = [0];	% dead time for each input
nA_in = [2];	% number of ARX coefficients

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
B_M1nz{1} = [-0.0025808548	0.013990402	0.016046807	0.0076725315	0.013818264	-0.02224317	-0.018345894	-0.0077741425	0.0016621134	0.02122031	0.029917236	-0.0171584	-0.006186414	-2.8886358E-4	-0.022286566	0.018192047	0.024405353	-0.001371922	0.022093866	-0.0067928135	-0.007630798	0.03013663	0.034939628	0.022182556	0.05026287	0.06018417	0.04318147	0.05500429	0.020111185	0.0875324	0.068798624	0.074690014	0.054420665	0.05929142	0.042838812	0.023785856	0.054993764	0.046548452	0.06514031	0.0036244579	0.026752854	0.014733824	-0.009007311	0.006005513	0.030603796	0.03335254	0.025345108	0.038781565	0.01110093	0.010752616	0.025224846	-0.0023830275	0.03866991	-0.0056319553	-0.010114491	-0.013590777	0.0028348197	-0.02671513	-6.0289825E-4	0.006202912	-0.038216017	-0.0013936176	-0.01672872	0.016334994	0.0027323093	0.019627392	-0.012720628	-0.0015307181	-0.003857929	0.019006817	0.013013185	0.016967693	0.003917471	0.01694666	0.0147841675	-0.005410764	-0.0070814495	0.0062767323	-0.028751835	-0.023305725	0.040213205	0.009831418	-0.0036822944	-0.008085642	0.0043403753	-0.024692442	-0.0017340357	0.028330687	0.036696456	8.439366E-4	-0.03544168	-0.032425504	-8.059525E-4	0.0068018665	0.006438884	0.02808579	0.045883257	-0.027090732	0.0051605576	-0.0033507505]';	% c/p coefficients from WaterMV as appropriate
%B_M1nz{2} = [2.5613592E-6	3.3532745E-6	3.284219E-6	6.6354494E-7	-6.021898E-7	-6.696386E-8	1.0366144E-7	5.83986E-7	6.789514E-7	8.421473E-7	2.842958E-7	5.572515E-7	-5.846483E-7	-3.060748E-7	5.7489666E-7	1.5921942E-6	1.4127825E-6	3.0548636E-6	1.5522578E-6	2.935479E-6	0.0077340635	0.008634036	0.008997529	0.009410912	0.009185033	0.008880123	0.00875175	0.008588336	0.008427757	0.008285991	0.00814496	0.008006543	0.007872757	0.0077425027	0.007613881	0.007487144	0.0073618987	0.0072400942	0.00711972	0.00700207	0.006885839	0.006772967	0.006662464	0.0065517584	0.0064420262	0.0063362233	0.0062307534	0.0061278734	0.006025811	0.005925554	0.0058287354	0.0057325065	0.0056378744	0.0055453493	0.0054544685	0.0053645247	0.0052750204	0.0051863785	0.005100284	0.005016728	0.004933666	0.004852871	0.004772598	0.0046952446	0.0046186466	0.004541454	0.004466844	0.004393119	0.0043187826	0.004246779	0.004175366	0.0041038315	0.0040362533	0.0039695227	0.003905126	0.0038424344	0.0037805082	0.003717745	0.0036567345	0.0035962476	0.0035361764	0.00347684	0.0034173436	0.0033592456	0.0033047984	0.0032517056	0.0031995208	0.0031479779	0.0030973295	0.0030465815	0.002995385	0.002943395	0.0028950996	0.0028470745	0.002799713	5.393168E-4	2.3570792E-4	8.820382E-5	-7.367972E-5	-5.39022E-5]';
%B_M1nz{3} = [2.507147E-7	-4.0059132E-7	-3.535021E-7	2.535453E-7	6.870822E-7	-1.2493934E-8	-5.936394E-7	-5.269842E-7	-9.6218116E-8	2.8944555E-7	-4.361731E-7	1.3434286E-7	-7.238902E-7	-2.0799084E-6	-1.72928E-6	-0.045139987	-0.04822752	-0.04806858	-0.048262104	-0.04474238	-0.04095992	-0.038428795	-0.03584219	-0.033430386	-0.031254802	-0.029214013	-0.027312702	-0.025543965	-0.023889594	-0.0223453	-0.020901956	-0.019552363	-0.01829121	-0.01711076	-0.016007766	-0.014975158	-0.014008489	-0.013104921	-0.012259579	-0.011469212	-0.002204335	-9.2775753E-4	-3.1037876E-4	3.3138672E-4	2.3260205E-4	4.8008074E-5	6.5979584E-5	4.085873E-5	1.9708918E-5	1.45616095E-5	9.951198E-6	6.8298777E-6	3.276344E-6	7.336694E-7	-2.3511237E-8	-1.4639858E-6	-1.8790179E-6	-1.9097838E-6	-1.4921928E-6	-9.640256E-7	-9.251409E-7	-1.2267969E-6	3.7330526E-7	1.8701534E-6	5.2908905E-7	-1.4980868E-7	-3.286104E-7	-7.2963843E-7	-1.03995966E-7	-5.720807E-7	-1.5128115E-6	-1.6867011E-6	-2.0221798E-6	-6.8605533E-7	-9.883147E-7	-1.143614E-7	1.3132421E-6	1.1160269E-6	1.5235787E-6	2.1523505E-8	-1.8618907E-6	-3.117529E-6	-3.0495178E-6	-2.4042872E-6	-1.0631715E-6	1.3430044E-7	7.171819E-7	1.7448627E-6	3.3933597E-6	2.0124992E-6	8.0332023E-7	-7.0053767E-7	-1.3147782E-6	-4.0148495E-7	-5.529108E-7	-1.5932477E-8	2.3128958E-7	1.2132157E-6	2.0484083E-6	1.800551E-6]';

A_M1nz{1} = -[-0.041426413	-0.012485306]';

% align B coefficient vectors to the same size (for later comparison)
for i_var = 1:n_var
	B_M1{i_var} = [zeros(1,dt_M1(i_var)) [B_M1nz{i_var}]' zeros(1,nB_in(i_var) - length([zeros(1,dt_M1(i_var)) [B_M1nz{i_var}]']))]';
end
A_M1{1} = [[A_M1nz{1}]' zeros(1,nA_in(1) - length(A_M1nz{1}))]';

% Output Prediction
if toggle_inc == 1;		% incremental model
	[Yt_model_M1] = tf_ARX_pred_inc2 (Ut_model_i, A_M1, B_M1, dt_M1, Yt_model);
	[Yv_model_M1] = tf_ARX_pred_inc2 (Uv_model_i, A_M1, B_M1, dt_M1, Yv_model);	
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
K_sign_M2 = [1 -1];	% specify sign direction vector
K_bound_M2 = [[zeros(1,n_var)] -0.5];
min_phase_M2 = zeros(1,n_var+1);	% 0 = not minimum phase

%[con_A_M2 con_b_M2] = CMI_con2(nZ_M2, K_sign_M2, K_bound_M2, min_phase_M2);
[con_A_M2 con_b_M2] = CMI_con(nZ_M2, K_sign_M2, min_phase_M2);	% set up inequality constraint matrices

% Model Identification and Output Prediction
if toggle_inc == 1;		% incremental model
	[B_M2 A_M2] = QP_ARX_MISO (Ut_model_i, Yt_model_i, nB_M2, dt_M2, nA_in, con_A_M2, con_b_M2);	% parameter estimation
	[Yt_model_M2] = tf_ARX_pred_inc2 (Ut_model_i, A_M2, B_M2, dt_M2, Yt_model);
	[Yv_model_M2] = tf_ARX_pred_inc2 (Uv_model_i, A_M2, B_M2, dt_M2, Yv_model);	
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
K_sign_M3 = [1 -1];	% specify sign direction vector
K_bound_M3 = [[zeros(1,n_var)] -0.5];
min_phase_M3 = ones(1,n_var+1);	% 0 = not minimum phase

% [con_A_M3 con_b_M3] = CMI_con2(nZ_M3, K_sign_M3, K_bound_M3, min_phase_M3);
[con_A_M3 con_b_M3] = CMI_con(nZ_M3, K_sign_M3, min_phase_M3);	% set up inequality constraint matrices

% Model Identification and Output Prediction
if toggle_inc == 1;		% incremental model
	[B_M3 A_M3] = QP_ARX_MISO (Ut_model_i, Yt_model_i, nB_M3, dt_M3,  nA_in, con_A_M3, con_b_M3);	% parameter estimation
	[Yt_model_M3] = tf_ARX_pred_inc2 (Ut_model_i, A_M3, B_M3, dt_M3, Yt_model);
	[Yv_model_M3] = tf_ARX_pred_inc2 (Uv_model_i, A_M3, B_M3, dt_M3, Yv_model);
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