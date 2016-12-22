% CMI_FIR_pt3.m
% Justification for constrained model identification - FIR models
% part 3: Model Identification

% variables from earlier parts: Ut, Yt_measured

% Data pre-treatment
toggle_filter = 1;	% 1 = apply filter, 0 = do not apply
	mov_avg = 20; 	% 1-d digital filter
	
toggle_ascale = 1; 	% 1 = apply auto-scaling, 0 = do not apply

if 	toggle_filter == 1;
	Yt_model = filter(ones(1,mov_avg)/mov_avg, 1, Yt_measured);
	Yv_model = filter(ones(1,mov_avg)/mov_avg, 1, Yv);
	
	for i_var = 1:n_var;
		B_mv{i_var} = filter(ones(1,mov_avg)/mov_avg, 1, B{i_var});
	end
else 
	Yt_model = Yt_measured;
	Yv_model = Yv;
end

if 	toggle_ascale == 1;
	[Yt_model, Yt_m, Yt_std] = ascale(Yt_model);
	Yv_model = (Yv_model - Yt_m) ./ Yt_std;
	[Ut_model, Ut_m, Ut_std] = ascale(Ut);
	Uv_model = (Uv - Ut_m) ./ Ut_std;
	
	% approximate FIR coefficients
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

MV_export = [Ut_model Yt_model];