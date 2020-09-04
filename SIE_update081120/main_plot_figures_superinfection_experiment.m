% function void = main_plot_figures_superinfection_experiment(void)

clear all; close all; clc;

%load('coinfection_parameter_estimates');
% load('superinfection_parameter_estimates');
% load('superinfection_parameter_estimates_withCIs');
load('superinfection_parameter_estimates_update081120');

estimated_actual_coinfection_moi_rH1N2 = paramsfit_co(1);
prob_gene_segment_expressed_PB2_PB1_PA = paramsfit_co(2);
prob_gene_segment_expressed_HA_rH3N1 = paramsfit_co(3);
prob_gene_segment_expressed_HA_rH1N2 = paramsfit_co(4);
fraction_x = paramsfit_co(5);
fraction_y = paramsfit_co(6);

%subplot(3,2,1);
figure(1); subplot(1,3,1);
semilogx(data_co.coinfection_actual_moi_rH3N1, data_co.coinfection_H3_positive_percent, 'ro'); hold on;
semilogx(data_co.coinfection_actual_moi_rH3N1, data_co.coinfection_H1_positive_percent, 'bo');
semilogx(data_co.coinfection_actual_moi_rH3N1, data_co.coinfection_H3_positive_H1_positive_percent, 'go'); 
semilogx(results_co.est_actual_moi_rH3N1, results_co.est_coinfection_H3_positive_percent, 'r');
semilogx(results_co.est_actual_moi_rH3N1, results_co.est_coinfection_H1_positive_percent, 'b');
semilogx(results_co.est_actual_moi_rH3N1, results_co.est_coinfection_H3_positive_H1_positive_percent, 'g');
title('coinfection experiment'); xlabel('actual moi rH3N1'); ylabel('proportion of cells expressing'); axis([10^-2 10^2 0 100]);
%legend('H3+ (data)', 'H1+ (data)', 'H3+ H1+ (data)', 'H3+ (fit)', 'H1+ (fit)', 'H3+ H1+ (fit)');

%subplot(3,2,2);
figure(2);
semilogx(data_co.coinfection_actual_moi_rH3N1, data_co.coinfection_percent, 'kd'); hold on;
results_co.coinfection_percent = 100*(results_co.est_coinfection_H3_positive_H1_positive_percent./results_co.est_coinfection_H3_positive_percent);
semilogx(results_co.est_actual_moi_rH3N1, results_co.coinfection_percent, 'k'); hold on;
xlabel('actual moi rH3N1'); ylabel('coinfection %'); axis([10^-2 10^2 0 100]);
%legend('data', 'fit');


%subplot(3,2,3);
figure(1); subplot(1,3,2);
semilogx(data_super.superinfection_actual_moi_rH3N1, data_super.superinfection_H3_positive_percent, 'ro'); hold on;
semilogx(data_super.superinfection_actual_moi_rH3N1, data_super.superinfection_H1_positive_percent, 'bo');
semilogx(data_super.superinfection_actual_moi_rH3N1, data_super.superinfection_H3_positive_H1_positive_percent, 'go'); 
semilogx(results_super.est_actual_moi_rH3N1, results_super.est_superinfection_H3_positive_percent_const, 'r--');
semilogx(results_super.est_actual_moi_rH3N1, results_super.est_superinfection_H1_positive_percent_const, 'b--');
semilogx(results_super.est_actual_moi_rH3N1, results_super.est_superinfection_H3_positive_H1_positive_percent_const, 'g--');
title('superinfection experiment - susceptibility lower for infected compared to uninfected cells'); 
xlabel('actual moi rH3N1'); ylabel('proportion of cells expressing'); axis([10^-2 10^2 0 100]);
%legend('H3+ (data)', 'H1+ (data)', 'H3+ H1+ (data)', 'H3+ (fit)', 'H1+ (fit)', 'H3+ H1+ (fit)');

%subplot(3,2,4);
figure(2);
semilogx(data_super.superinfection_actual_moi_rH3N1, data_super.superinfection_percent, 'kd'); hold on;
results_super.superinfection_percent_const = 100*(results_super.est_superinfection_H3_positive_H1_positive_percent_const./results_super.est_superinfection_H3_positive_percent_const);
semilogx(results_super.est_actual_moi_rH3N1, results_super.superinfection_percent_const, 'k'); hold on;
xlabel('actual moi rH3N1'); ylabel('superinfection %'); axis([10^-2 10^2 0 100]);
%legend('data', 'fit');

%subplot(3,2,5);
figure(1); subplot(1,3,2);
semilogx(data_super.superinfection_actual_moi_rH3N1, data_super.superinfection_H3_positive_percent, 'ro'); hold on;
semilogx(data_super.superinfection_actual_moi_rH3N1, data_super.superinfection_H1_positive_percent, 'bo');
semilogx(data_super.superinfection_actual_moi_rH3N1, data_super.superinfection_H3_positive_H1_positive_percent, 'go'); 
semilogx(results_super.est_actual_moi_rH3N1, results_super.est_superinfection_H3_positive_percent_prop, 'r');
semilogx(results_super.est_actual_moi_rH3N1, results_super.est_superinfection_H1_positive_percent_prop, 'b');
semilogx(results_super.est_actual_moi_rH3N1, results_super.est_superinfection_H3_positive_H1_positive_percent_prop, 'g');
title('superinfection experiment - susceptibility decreases with cellular rH3N1 MOI'); 
xlabel('actual moi rH3N1'); ylabel('proportion of cells expressing'); axis([10^-2 10^2 0 100]);
legend('H3+ (data)', 'H1+ (data)', 'H3+ H1+ (data)', 'H3+ (fit)', 'H1+ (fit)', 'H3+ H1+ (fit)');

%subplot(3,2,6);
figure(2);
semilogx(data_super.superinfection_actual_moi_rH3N1, data_super.superinfection_percent, 'kd'); hold on;
results_super.superinfection_percent_prop = 100*(results_super.est_superinfection_H3_positive_H1_positive_percent_prop./results_super.est_superinfection_H3_positive_percent_prop);
semilogx(results_super.est_actual_moi_rH3N1, results_super.superinfection_percent_prop, 'k'); hold on;
legend('% superinfection - data', '% superinfection - estimated');
xlabel('actual moi rH3N1'); ylabel('superinfection %'); axis([10^-2 10^2 0 100]);
legend('data', 'fit');

%figure; 

i = 0:1:10;

s = paramsfit_super(1)
r = paramsfit_super(2)

%subplot(3,2,3);
figure(1); subplot(1,3,3);
susc_level_input_indep = ones(size(i))*s; susc_level_input_indep(1) = 1;
susc_level_input_dep = r.^i;
plot(i, susc_level_input_indep, 'k--'); hold on;
plot(i, susc_level_input_dep, 'k');
legend('input independent model', 'input dependent model');
xlabel('cellular input i'); ylabel('susceptibility')

