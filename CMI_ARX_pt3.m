% CMI_ARX_pt3.m
% Justification for constrained model identification - ARX models
% part 3: Data Pre-treatment

% variables from earlier parts: Ut, Yt_measured

% set pre-treatment settings
mov_avg = 20; 	% 1-d digital filter (set 1 to switch off)
toggle_ascale = 1; 	% 1 = apply auto-scaling, 0 = do not apply
toggle_inc = 1; 	% 0 = use absolute model, 1 = use incremental model

% moving average
Yt_model = filter(ones(1,mov_avg)/mov_avg, 1, Yt_measured);
Yv_model = filter(ones(1,mov_avg)/mov_avg, 1, Yv);

Yt_val = Yt_model;	% crisp output in Eng units for model validation (pt5)
Yv_val = Yv_model;

% auto-scaling
if 	toggle_ascale == 1;
	[Yt_model, Yt_m, Yt_std] = ascale(Yt_model);
	Yv_model = (Yv_model - Yt_m) ./ Yt_std;
	[Ut_model, Ut_m, Ut_std] = ascale(Ut);
	Uv_model = (Uv - Ut_m) ./ Ut_std;
else
	Ut_model = Ut;
	Uv_model = Uv;
end

% if incremental model is used, sample values become change
if toggle_inc == 1;
	Ut_model_i = abs2inc(Ut_model);
	Yt_model_i = abs2inc(Yt_model);
	Uv_model_i = abs2inc(Uv_model);
	Yv_model_i = abs2inc(Yv_model);	
	
	Yt_val_i = abs2inc(Yt_val);
	Yv_val_i = abs2inc(Yv_val);
end

Yt_model_s1 = shiftdwn(Yt_model,1); % ARX component for pharmaMV

% export to WaterMV
MV_export = [Ut_model Yt_model_s1 Yt_model];