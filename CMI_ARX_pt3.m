% CMI_ARX_pt3.m
% Justification for constrained model identification - ARX models
% part 3: Data Pre-treatment

% variables from earlier parts: Ut, Yt_measured

% set toggles
toggle_filter = 1;	% 1 = apply filter, 0 = do not apply
	mov_avg = 1; 	% 1-d digital filter
	
toggle_ascale = 1; 	% 1 = apply auto-scaling, 0 = do not apply

% moving average
if 	toggle_filter == 1;
	Yt_model = filter(ones(1,mov_avg)/mov_avg, 1, Yt_measured);
	Yv_model = filter(ones(1,mov_avg)/mov_avg, 1, Yv);
	
	% using a moving average affects the FIR/ARX coefficients, this is to account for that
	for i_var = 1:n_var;
		B_mv{i_var} = filter(ones(1,mov_avg)/mov_avg, 1, B{i_var});
		A_mv{1} = filter(ones(1,mov_avg)/mov_avg, 1, A{1});
	end
else 
	Yt_model = Yt_measured;
	Yv_model = Yv;
end

% auto-scaling
if 	toggle_ascale == 1;
	[Yt_model, Yt_m, Yt_std] = ascale(Yt_model);
	Yv_model = (Yv_model - Yt_m) ./ Yt_std;
	[Ut_model, Ut_m, Ut_std] = ascale(Ut);
	Uv_model = (Uv - Ut_m) ./ Ut_std;
	
	% approximate ARX coefficients
	for i_var = 1:n_var
		B_as{i_var} = B{i_var} * (Ut_std(i_var)/Yt_std);
		if 	toggle_filter == 1;
			B_mv_as{i_var} = B_mv{i_var} * (Ut_std(i_var)/Yt_std);
		end
	end
else
	Ut_model = Ut;
	Uv_model = Uv;
end

% export to WaterMV
MV_export = [Ut_model Yt_model];