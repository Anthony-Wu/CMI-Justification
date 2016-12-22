% AbsvsInc_pt2
colour_true = [0.0 0.0 0.0];
colour_i = [0.9 0.2 0.0];
colour_a = [0.2 0.6 0.6];
i_plot = 0;

% Engineering Units (Eng)
U_E_i = abs2inc(U_E);
Y_E_i = abs2inc(Y_E);

% Input-Output Plot (Abs)
i_plot = i_plot + 1;
figure(i_plot)
set(gca,'fontsize',24);
subplot(2,1,1);
	stairs(U_E,'LineWidth',2.0);
	title(sprintf('Input (Eng)'),'FontWeight','bold');
	xlabel('Sample','FontWeight','bold') % x-axis label
	ylabel('Value (Eng)','FontWeight','bold') % y-axis label
subplot(2,1,2);
	plot([1:n_sam]',Y_E,'LineWidth',2.0,'Color',colour_true);
	title(sprintf('Output (Eng)'),'FontWeight','bold');
	xlabel('Sample','FontWeight','bold') % x-axis label
	ylabel('Value (Eng)','FontWeight','bold') % y-axis label

% Input-Output Plot (Inc)
i_plot = i_plot + 1;
figure(i_plot)
set(gca,'fontsize',24);
subplot(2,1,1);
	stairs(U_E_i,'LineWidth',2.0);
	title(sprintf('Input (Eng)'),'FontWeight','bold');
	xlabel('Sample','FontWeight','bold') % x-axis label
	ylabel('Value (Eng)','FontWeight','bold') % y-axis label
subplot(2,1,2);
	plot([1:n_sam]',Y_E_i,'LineWidth',2.0,'Color',colour_true);
	title(sprintf('Output (Eng)'),'FontWeight','bold');
	xlabel('Sample','FontWeight','bold') % x-axis label
	ylabel('Value (Eng)','FontWeight','bold') % y-axis label	
	
% Model Identification (Eng)
[B_E_p_i] = QP_FIR_MISO (U_E_i, Y_E_i, nB, 0, [], []);
Y_E_p_i = tf_FIR_pred(U_E_i, B_E_p_i, 0);

Y_E_p_i2 = zeros(n_sam,1);
Y_E_p_i2(nB,1) = Y_E(nB,1);
for i = nB+1:n_sam;
Y_E_p_i2(i) = Y_E_p_i2(i-1) + Y_E_p_i(i);
end

[B_E_p_a] = QP_FIR_MISO (U_E, Y_E, nB, 0, [] ,[]);
Y_E_p_a = tf_FIR_pred(U_E, B_E_p_a, 0);

RMSE_E_i = stat_RMSE (Y_E(nB+1:end), Y_E_p_i2(nB+1:end));
RMSE_E_a = stat_RMSE (Y_E(nB+1:end), Y_E_p_a(nB+1:end));

% FIR coefficient plot (Eng)
i_plot = i_plot + 1;
figure(i_plot)
bar_ref = bar([B_E{1}, B_E_p_i{1} , B_E_p_a{1}]);
	bar_ref(1).FaceColor = colour_true;
	bar_ref(1).BarWidth = 1.0;	
	bar_ref(2).FaceColor = colour_i;
	bar_ref(2).BarWidth = 1.0;
	bar_ref(3).FaceColor = colour_a;
	bar_ref(3).BarWidth = 1.0;
set(gca,'fontsize',24);
title(sprintf('FIR Coefficients (Eng)'),'FontWeight','bold');
legend('Actual', 'Incremental', 'Absolute');
xlabel('Coefficient number','FontWeight','bold') % x-axis label
ylabel('Value','FontWeight','bold') % y-axis label

% Output Prediction (Eng)
i_plot = i_plot + 1;
figure(i_plot)
pl1 = plot([1:n_sam]',Y_E,'Color',colour_true,'LineWidth', 2.0);
hold on;
pl2 = plot([1:n_sam]', Y_E_p_i2,'Color',colour_i,'LineWidth', 2.0,'Linestyle','--');
pl3 = plot([1:n_sam]', Y_E_p_a,'Color',colour_a,'LineWidth', 2.0,'Linestyle',':');
hold off;
title(sprintf('Output Prediction (Eng) RMSE: Inc: %.2f, Abs: %.2f',RMSE_E_i, RMSE_E_a));
legend('Measured','Incremental','Absolute');
set(gca,'fontsize',24);
xlabel('Sample','FontWeight','bold') % x-axis label
ylabel('Value (Eng)','FontWeight','bold') % y-axis label


% Auto-Scaled Units (AS)
U_m = mean(U_E);
U_s = std(U_E);
Y_m = mean(Y_E);
Y_s = std(Y_E);

U_AS = (U_E - U_m) ./ U_s;
Y_AS = (Y_E - Y_m) ./ Y_s;

U_AS_i = abs2inc(U_AS);
Y_AS_i = abs2inc(Y_AS);

% Input-Output Plot (Abs)
i_plot = i_plot + 1;
figure(i_plot)
set(gca,'fontsize',24);
subplot(2,1,1);
	stairs(U_AS,'LineWidth',2.0);
	title(sprintf('Input(AS)'),'FontWeight','bold');
	xlabel('Sample','FontWeight','bold') % x-axis label
	ylabel('Value (AS)','FontWeight','bold') % y-axis label
subplot(2,1,2);
	plot([1:n_sam]',Y_AS,'LineWidth',2.0,'Color',colour_true);
	title(sprintf('Output (AS)'),'FontWeight','bold');
	xlabel('Sample','FontWeight','bold') % x-axis label
	ylabel('Value (AS)','FontWeight','bold') % y-axis label

% Input-Output Plot (Inc)
i_plot = i_plot + 1;
figure(i_plot)
set(gca,'fontsize',24);
subplot(2,1,1);
	stairs(U_AS_i,'LineWidth',2.0);
	title(sprintf('Input (AS)'),'FontWeight','bold');
	xlabel('Sample','FontWeight','bold') % x-axis label
	ylabel('Value (AS)','FontWeight','bold') % y-axis label
subplot(2,1,2);
	plot([1:n_sam]',Y_AS_i,'LineWidth',2.0,'Color',colour_true);
	title(sprintf('Output (AS)'),'FontWeight','bold');
	xlabel('Sample','FontWeight','bold') % x-axis label
	ylabel('Value (AS)','FontWeight','bold') % y-axis label
	
[B_AS_p_a] = QP_FIR_MISO (U_AS, Y_AS, nB, 0,[],[]);
Y_AS_p_a = tf_FIR_pred(U_AS, B_AS_p_a, 0);

[B_AS_p_i] = QP_FIR_MISO (U_AS_i, Y_AS_i, nB, 0,[],[]);
Y_AS_p_i = tf_FIR_pred(U_AS_i, B_AS_p_i, 0);

Y_AS_p_i2 = zeros(n_sam,1);
Y_AS_p_i2(nB,1) = Y_AS(nB,1);
for i = nB+1:n_sam;
Y_AS_p_i2(i) = Y_AS_p_i2(i-1) + Y_AS_p_i(i);
end

RMSE_AS_i = stat_RMSE (Y_AS(nB+1:end), Y_AS_p_i2(nB+1:end));
RMSE_AS_a = stat_RMSE (Y_AS(nB+1:end), Y_AS_p_a(nB+1:end));

% FIR coefficient plot (AS)
i_plot = i_plot + 1;
figure(i_plot)
bar_ref = bar([B_AS_p_i{1} , B_AS_p_a{1}]);
	bar_ref(1).FaceColor = colour_i;
	bar_ref(1).BarWidth = 1.0;
	bar_ref(2).FaceColor = colour_a;
	bar_ref(2).BarWidth = 1.0;
set(gca,'fontsize',24);
title(sprintf('FIR Coefficients (AS)'),'FontWeight','bold');
legend('Incremental','Absolute');
xlabel('Coefficient number','FontWeight','bold') % x-axis label
ylabel('Value','FontWeight','bold') % y-axis label

% Output Prediction (AS)
i_plot = i_plot + 1;
figure(i_plot)
pl1 = plot([1:n_sam]',Y_AS,'Color',colour_true,'LineWidth', 2.0);
hold on;
pl2 = plot([1:n_sam]', Y_AS_p_i2,'Color',colour_i,'LineWidth', 2.0,'Linestyle','--');
pl3 = plot([1:n_sam]', Y_AS_p_a,'Color',colour_a,'LineWidth', 2.0,'Linestyle',':');
hold off;
title(sprintf('Output Prediction (AS) RMSE: Inc: %.2f, Abs: %.2f',RMSE_AS_i, RMSE_AS_a));
legend('Measured','Incremental','Absolute');
set(gca,'fontsize',24);
xlabel('Sample','FontWeight','bold') % x-axis label
ylabel('Value (AS)','FontWeight','bold') % y-axis label