% function void = main_estimate_coinfection(void)

% updated 05/12/20 to include 95% confidence intervals for actual rH1N2 MOI, x, and y 

% iterations display - max 0f 100 iterations; can change options in 'GetEstimated_percents_coinfection_fit_fraction'

% Two types of cells - low S and high S (susceptibility).
% The low S class is y fraction of the cell population, the high S class is (1-y) fraction.
% The low S class carries fraction x of the viral population, the high S class carries the remaining fraction (1-x).
% Then, MOI of entire population = V/C, where V is number of virions and C is total number of cells
% MOI of low S population = (V*x)/(C*y) = MOI*(x/y)
% MOI of high S population = (V*(1-x)/(C*(1-y)) = MOI*(1-x)/(1-y)

% rH3N1 is the first virus, doing the excluding
% rH1N2 is the second virus, subjected to exclusion

clear all; close all; clc;

data_co.coinfection_intended_moi_rH3N1 = [ones(1,3)*0.05 ones(1,3)*0.25 ones(1,3)*1 ones(1,3)*2.5 ones(1,3)*10];
data_co.coinfection_intended_moi_rH1N2 = 0.5;
data_co.coinfection_actual_moi_rH3N1 = [0.256683124 0.238288221 0.260457452 1.113392459 1.146690531 1.13953782 2.374098681 3.575872797 3.025801283 13.47316799 17.98389388 14.81802829 49.67490602 63.91834109 47.56707034];

data_co.coinfection_H3_positive_percent = [26.48 27.77 27.47 58.64 59.18 56.78 75.6 77.5 76.3 88.7 87.6 87.8 97.7 97.6 97.4];
data_co.coinfection_H1_positive_percent = [61.8 64.8 64.5 69.9 69.6 67.7 67.65 68 67.25 68.58 64.89 66.53 63.99 62.73 62.95];
data_co.coinfection_H3_positive_H1_positive_percent = [24 25.4 25.1 51.4 51.9 49.5 60.3 61.5 60.5 65.4 61.9 63.4 63.7 62.5 62.7];
data_co.coinfection_percent = 100*data_co.coinfection_H3_positive_H1_positive_percent./data_co.coinfection_H3_positive_percent;

% estimates: exogenous same pp constraints!
estimated_actual_coinfection_moi_rH1N2 = 1.6;
prob_gene_segment_expressed_PB2_PB1_PA = 0.95;
prob_gene_segment_expressed_HA_rH3N1 = prob_gene_segment_expressed_PB2_PB1_PA;
prob_gene_segment_expressed_HA_rH1N2 = prob_gene_segment_expressed_PB2_PB1_PA;
% this should be replaced by values that correspond to the r value in Table 1 (if possible):
fraction_x = 0.001; % this fraction of the viral population...
fraction_y = 0.06; % ...goes into this proportion of the cells

% k = 0.3076, CI[0.1576, 0.6001] 
% b = 0.4472, CI[0.2674, 0.7478] 
% mu = 0.0146, CI[0.0099, 0.0215] 

params_fixed.death_k = 0.307563294635010; % shape parameter of Weibull hazard
params_fixed.death_b = 0.447205463276431; % scale parameter of Weibull hazard - alpha(t,i) = bk(t^(k-1))
params_fixed.death_mu = 0.0145642125211200;


params_fixed.lambda = exp(log(params_fixed.death_b)/(-params_fixed.death_k));
params_fixed.prob_infected_cell_dead_19hpi = 1-exp(-(19/params_fixed.lambda)^params_fixed.death_k - 19*params_fixed.death_mu);
params_fixed.prob_infected_cell_dead_6hpi = 1-exp(-(6/params_fixed.lambda)^params_fixed.death_k - 6*params_fixed.death_mu);


% for fitting, only estimate at observed actual moi of rH3N1:
params_init = [estimated_actual_coinfection_moi_rH1N2 prob_gene_segment_expressed_PB2_PB1_PA prob_gene_segment_expressed_HA_rH3N1 prob_gene_segment_expressed_HA_rH1N2 fraction_x fraction_y];

[paramsfit_co, est_coinfection_H3_positive_percent_fit, est_coinfection_H1_positive_percent_fit, est_coinfection_H3_positive_H1_positive_percent_fit] = GetEstimated_percents_coinfection_fit_fraction(data_co, params_init, params_fixed);

estimated_actual_coinfection_moi_rH1N2 = paramsfit_co(1);
prob_gene_segment_expressed_PB2_PB1_PA = paramsfit_co(2);
prob_gene_segment_expressed_HA_rH3N1 = paramsfit_co(3);
prob_gene_segment_expressed_HA_rH1N2 = paramsfit_co(4);
fraction_x = paramsfit_co(5);
fraction_y = paramsfit_co(6);

results_co.actual_coinfection_moi_rH1N2 = estimated_actual_coinfection_moi_rH1N2;
results_co.fraction_x = fraction_x;
results_co.fraction_y = fraction_y;

% get 95% confidence intervals for the six parameters being fit
results_co.CI_estimates = get_confidinterval_co(data_co, paramsfit_co, params_fixed);

results_co.est_actual_moi_rH3N1 = transpose(10.^[-2:0.1:2]);

[est_coinfection_H3_positive_percent, est_coinfection_H1_positive_percent, est_coinfection_H3_positive_H1_positive_percent] = GetEstimated_percents_coinfection_sim_fraction(transpose(results_co.est_actual_moi_rH3N1), paramsfit_co, params_fixed);

results_co.est_coinfection_H3_positive_percent = est_coinfection_H3_positive_percent';
results_co.est_coinfection_H1_positive_percent = est_coinfection_H1_positive_percent';
results_co.est_coinfection_H3_positive_H1_positive_percent = est_coinfection_H3_positive_H1_positive_percent';

if 1
    
    save('coinfection_parameter_estimates_update081120', 'data_co', 'paramsfit_co', 'results_co', 'params_fixed')
    
end