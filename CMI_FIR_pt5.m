% CMI_FIR_pt5.m
% Justification for constrained model identification - FIR models
% part 5: Model Validation

close all;	% close all figures
i_plot = 0;	% reset plot index

colour_true = [0 0 0];
colour_M1 = [0.4 0.8 0.4];
colour_M2 = [1 0.4 0.4];
colour_M3 = [0.6 0.6 1];

% plot Input - Output contribution training
for i_var = 1:n_var;
	i_plot = i_plot +1;
	figure (i_plot);
	
	subplot(2,1,1);
	stairs(Ut(:,i_var),'LineWidth',1.5);
	title(sprintf('Training Data Input %d (Eng)',i_var),'FontWeight','bold');
	xlabel('Sample','FontWeight','bold') % x-axis label
	ylabel('Value','FontWeight','bold') % y-axis label
	
	subplot(2,1,2);
	if n_var == 1;
		plot([1:n_sam_t]',Yt,'LineWidth',1.5,'Color',colour_true);
		title(sprintf('Training Data Output (no noise, Eng)'),'FontWeight','bold');	
	else
		plot([1:n_sam_t]',Yt_cont(:,i_var),'LineWidth',1.5,'Color',[0.8 0.4 0],'Linestyle','--');
		title(sprintf('Output contribution from Input %d (Eng)',i_var),'FontWeight','bold');
	end
	xlabel('Sample','FontWeight','bold') % x-axis label
	ylabel('Value','FontWeight','bold') % y-axis label
end

% plot Output cont. - Output training
if 	n_var > 1;	% no need to generate figures if SISO problem
	ppf = 2;	% plots-per-figure
	n_cyc = ceil((n_var)/ppf);	% number of iteration cycles
	
	% set common axis for better visual comparison
	ymax = ceil(max(max([Yt_cont Yt]))/10)*10;
	ymin = floor(min((min([Yt_cont Yt]))/10)*10);
	
	i_var = 0;
	for i_cyc = 1:n_cyc;
		i_plot = i_plot+1;
		figure (i_plot);
		
		for i_ppf = 1:ppf;
			if	(i_cyc-1) * ppf + i_ppf  > n_var % if the number of inputs are covered, stop
			else
				subplot(ppf,1,i_ppf);
				plot([1:n_sam_t]',Yt_cont(:,i_var+i_ppf),'LineWidth',1.5,'Color',[0.8 0.4 0],'Linestyle','--');
				title(sprintf('Output contribution from Input %d (Eng)',i_var+i_ppf),'FontWeight','bold');
				xlabel('Sample','FontWeight','bold') % x-axis label
				ylabel('Value','FontWeight','bold') % y-axis label
				ylim([ymin ymax]);
				
				i_var = i_var + 1;
			end
		end
	end
	
	if mod(n_var,ppf) == 0;	% if ppf is a factor of n_var, figures would be completely filled by subplots, need overall output in need subplot
		i_plot = i_plot+1;
		figure (i_plot);
		
		subplot(ppf,1,1);
		plot([1:n_sam_t]',Yt,'LineWidth',1.5,'Color',colour_true);
		title(sprintf('Training Data Output (no noise, Eng)'),'FontWeight','bold');	
		xlabel('Sample','FontWeight','bold') % x-axis label
		ylabel('Value','FontWeight','bold') % y-axis label
		ylim([ymin ymax]);
	else % ppf is not a factor of n_var, but overall output in a subplot of the current figure
		subplot(ppf,1,i_ppf);
		plot([1:n_sam_t]',Yt,'LineWidth',1.5,'Color',[colour_true]);
		title(sprintf('Training Data Output (no noise, Eng)'),'FontWeight','bold');	
		xlabel('Sample','FontWeight','bold') % x-axis label
		ylabel('Value','FontWeight','bold') % y-axis label
		ylim([ymin ymax]);
	end
end

% model vs model with noise
i_plot = i_plot+1;
figure (i_plot);

if 	toggle_Y_noise == 1;
	pl1 = plot([1:n_sam_t]',Yt_n,'Color',[0.9 0 0],'LineWidth', 1.5,'Linestyle',':');
	hold on;
	pl2 = plot([1:n_sam_t]',Yt,'Color',[colour_true],'LineWidth', 1.5);
	hold off;
	legend('Measured (noise added)', 'Original (no noise)')
	title(sprintf('Measured Output (Eng)'),'FontWeight','bold');	
	xlabel('Sample','FontWeight','bold') % x-axis label
	ylabel('Value','FontWeight','bold') % y-axis label	
else
	pl2 = plot([1:n_sam_t]',Yt,'Color',[colour_true],'LineWidth', 1.5);
	title(sprintf('Measured Output (Eng)'),'FontWeight','bold');	
	xlabel('Sample','FontWeight','bold') % x-axis label
	ylabel('Value','FontWeight','bold') % y-axis label		
end

% FIR coefficients
% M1
for i_var = 1:n_var
	i_plot = i_plot + 1;
	figure (i_plot);
	if toggle_ascale == 1;
		bar_ref = bar([B_true_as{i_var}, B_M1{i_var}, B_true_mv_as{i_var}]);
		title(sprintf('FIR coefficients of input %d (AS)',i_var),'FontWeight','bold'); % add title
	elseif toggle_ascale == 0;
		bar_ref = bar([B_true{i_var}, B_M1{i_var}, B_true_mv{i_var}]);
		title(sprintf('FIR coefficients of input %d (Eng)',i_var),'FontWeight','bold'); % add title
	end
	bar_ref(1).FaceColor = colour_true;
	bar_ref(1).BarWidth = 1.0;
	bar_ref(2).FaceColor = colour_M1;
	bar_ref(2).BarWidth = 1.0;
	
	legend('Actual','UNC');
	xlabel('Coefficient number','FontWeight','bold') % x-axis label
	ylabel('Value','FontWeight','bold') % y-axis label
	xlim([0.5 nB_in(i_var)+0.5]);
end

% M2
for i_var = 1:n_var
	i_plot = i_plot + 1;
	figure (i_plot);
	if toggle_ascale == 1;
		bar_ref = bar([B_true_as{i_var}, B_M2{i_var}, B_true_mv_as{i_var}]);
		title(sprintf('FIR coefficients of input %d (AS)',i_var),'FontWeight','bold'); % add title
	elseif toggle_ascale == 0;
		bar_ref = bar([B_true{i_var}, B_M2{i_var}, B_true_mv{i_var}]);
		title(sprintf('FIR coefficients of input %d (Eng)',i_var),'FontWeight','bold'); % add title
	end
	bar_ref(1).FaceColor = colour_true;
	bar_ref(1).BarWidth = 1.0;
	bar_ref(2).FaceColor = colour_M2;
	bar_ref(2).BarWidth = 1.0;
	
	legend('Actual','CON1');
	xlabel('Coefficient number','FontWeight','bold') % x-axis label
	ylabel('Value','FontWeight','bold') % y-axis label
	xlim([0.5 nB_in(i_var)+0.5]);
end

% M3
for i_var = 1:n_var
	i_plot = i_plot + 1;
	figure (i_plot);
	if toggle_ascale == 1;
		bar_ref = bar([B_true_as{i_var}, B_M3{i_var}, B_true_mv_as{i_var}]);
		title(sprintf('FIR coefficients of input %d (AS)',i_var),'FontWeight','bold'); % add title
	elseif toggle_ascale == 0;
		bar_ref = bar([B_true{i_var}, B_M3{i_var}, B_true_mv{i_var}]);
		title(sprintf('FIR coefficients of input %d (Eng)',i_var),'FontWeight','bold'); % add title
	end
	bar_ref(1).FaceColor = colour_true;
	bar_ref(1).BarWidth = 1.0;
	bar_ref(2).FaceColor = colour_M3;
	bar_ref(2).BarWidth = 1.0;
	
	legend('Actual','CON2');
	xlabel('Coefficient number','FontWeight','bold') % x-axis label
	ylabel('Value','FontWeight','bold') % y-axis label
	xlim([0.5 nB_in(i_var)+0.5]);
end

% plot Yt_model
i_plot = i_plot +1;
figure (i_plot);
for i_var = 1:n_var;
	pl1 = plot([1:n_sam_t]',Yt_M1,'Color',colour_M1,'LineWidth', 1.5,'Linestyle','-.');
	hold on;
	pl2 = plot([1:n_sam_t]',Yt_M2,'Color',colour_M2,'LineWidth', 1.5,'Linestyle','--');
	pl3 = plot([1:n_sam_t]',Yt_M3,'Color',colour_M3,'LineWidth', 1.5,'Linestyle',':');
	pl10 = plot([1:n_sam_t]',Yt_val,'Color',colour_true,'LineWidth', 1.5);
	hold off;
	title(sprintf('Output Prediction (Training, Eng) RMSE: UNC=%.2f; CON1=%.2f; CON2=%.2f', RMSE_Yt_M1, RMSE_Yt_M2, RMSE_Yt_M3),'FontWeight','bold'); % add title
	legend('UNC','CON1', 'CON2','Actual');	% add legend
	xlabel('Sample','FontWeight','bold') % x-axis label
	ylabel('Value','FontWeight','bold') % y-axis label
end

% plot Yv_model
i_plot = i_plot +1;
figure (i_plot);
for i_var = 1:n_var;
	pl1 = plot([1:n_sam_v]',Yv_M1,'Color',colour_M1,'LineWidth', 1.5,'Linestyle','-.');
	hold on;
	pl2 = plot([1:n_sam_v]',Yv_M2,'Color',colour_M2,'LineWidth', 1.5,'Linestyle','--');
	pl3 = plot([1:n_sam_v]',Yv_M3,'Color',colour_M3,'LineWidth', 1.5,'Linestyle',':');
	pl10 = plot([1:n_sam_v]',Yv_val,'Color',colour_true,'LineWidth', 1.5);
	hold off;
	title(sprintf('Output Prediction (Validation, Eng) RMSE: UNC=%.2f; CON1=%.2f; CON2=%.2f', RMSE_Yv_M1, RMSE_Yv_M2, RMSE_Yv_M3),'FontWeight','bold'); % add title
	legend('UNC','CON1', 'CON2','Actual');	% add legend
	xlabel('Sample','FontWeight','bold') % x-axis label
	ylabel('Value','FontWeight','bold') % y-axis label
end