%   function void = main_estimate_superinfection_update081120(void)

% updated 05/12/20 to include 95% confidence intervals for s and r

% Two types of cells - low S and high S (susceptibility).
% The low S class is y fraction of the cell population, the high S class is (1-y) fraction.
% The low S class carries fraction x of the viral population, the high S class carries the remaining fraction (1-x).
% Then, MOI of entire population = V/C, where V is number of virions and C is total number of cells
% MOI of low S population = (V*x)/(C*y) = MOI*(x/y)
% MOI of high S population = (V*(1-x)/(C*(1-y)) = MOI*(1-x)/(1-y)

% rH3N1 is the first virus, doing the excluding
% rH1N2 is the second virus, subjected to exclusion

clear all; close all; clc;


load('coinfection_parameter_estimates_update081120'); %coinfection_parameter_estimates

data_super.superinfection_intended_moi_rH3N1 = [ones(1,3)*0.05 ones(1,3)*0.25 ones(1,3)*1 ones(1,3)*2.5 ones(1,3)*10];
data_super.superinfection_intended_moi_rH1N2 = 0.5;
data_super.superinfection_actual_moi_rH3N1 = [0.258087545 0.258600291 0.259107147 1.254055654 1.133384838 1.025532747 3.678579616 3.16619692 3.093577427 14.59138541 17.13254373 15.39865285 48.78388203 39.38645332 60.49241813];

data_super.superinfection_H3_positive_percent = [27.9 29 31.5 59.7 60.3 62.4 78.96 80.72 82.88 90.71 92.41 91.69 98.43 98.18 98.52];
data_super.superinfection_H1_positive_percent = [56 56.4 58.9 33.3 32.6 31.6 7.86 8.87 6.36 3.27 4.2 3.04 0.24 0.194 0.225];
data_super.superinfection_H3_positive_H1_positive_percent = [13.4 14.7 15.6 14.3 13.5 14.3 4.26 4.92 3.68 2.21 3.01 2.09 0.23 0.18 0.22];
data_super.superinfection_percent = 100*data_super.superinfection_H3_positive_H1_positive_percent./data_super.superinfection_H3_positive_percent;

n_pts = length(data_super.superinfection_H3_positive_percent)+length(data_super.superinfection_H1_positive_percent)+length(data_super.superinfection_H3_positive_H1_positive_percent);

% params_fixed.death_k = 0.307563294635010; % shape parameter of Weibull hazard
% params_fixed.death_b = 0.447205463276431; % scale parameter of Weibull hazard - alpha(t,i) = bk(t^(k-1))
% params_fixed.death_mu = 0.0145642125211200;

params_fixed.lambda = exp(log(params_fixed.death_b)/(-params_fixed.death_k));
params_fixed.prob_infected_cell_dead_19hpi = 1-exp(-(19/params_fixed.lambda)^params_fixed.death_k - 19*params_fixed.death_mu);
params_fixed.prob_infected_cell_dead_6hpi = 1-exp(-(6/params_fixed.lambda)^params_fixed.death_k - 6*params_fixed.death_mu);

disp('fitting input-independent susceptibility model to superinfection data...'); fprintf('\n');

% fitting constant model:
model_num = 0;                    % 0 = constant model; 1 = proportional reduction
susc_reduction_const = 0.04;      % in constant model, infected cells have a constant level of reduced susceptibility relative to uninfected cells
n_parsfit_m0 = length(susc_reduction_const);

[susc_reduction_const_fit, results_super.LSE_m0, est_superinfection_H3_positive_percent_fit, est_superinfection_H1_positive_percent_fit, est_superinfection_H3_positive_H1_positive_percent_fit] = GetEstimated_percents_superinfection_fit_fraction(data_super, paramsfit_co, model_num, susc_reduction_const, params_fixed);

% results_super.est_superinfection_H3_positive_percent_fit = est_superinfection_H3_positive_percent_fit';
% results_super.est_superinfection_H1_positive_percent_fit = est_superinfection_H1_positive_percent_fit';
% results_super.est_superinfection_H3_positive_H1_positive_percent_fit = est_superinfection_H3_positive_H1_positive_percent_fit';

%susc_reduction_const_fit = susc_reduction_const;
paramsfit_super(1) = susc_reduction_const_fit;
results_super.susc_reduction_const_fit = susc_reduction_const_fit;
disp(['input-independent susceptibility model estimate: s = ', num2str(susc_reduction_const_fit,'%2.4f')]); fprintf('\n');

AIC_value_m0 = n_pts*log(results_super.LSE_m0)+2*n_parsfit_m0;

disp('getting 95% confidence intervals...');

% get 95% confidence intervals
results_super.CI_estimate_m0 = get_confidinterval_super(data_super, paramsfit_co, model_num, log(susc_reduction_const_fit), params_fixed)

% get percent positive curves with best-fit parameters
results_super.est_actual_moi_rH3N1 = transpose(10.^[-2:0.1:2]);
[est_superinfection_H3_positive_percent_const, est_superinfection_H1_positive_percent_const, est_superinfection_H3_positive_H1_positive_percent_const] = GetEstimated_percents_superinfection_sim_fraction(transpose(results_super.est_actual_moi_rH3N1), paramsfit_co, model_num, paramsfit_super(1), params_fixed);

results_super.est_superinfection_H3_positive_percent_const = est_superinfection_H3_positive_percent_const';
results_super.est_superinfection_H1_positive_percent_const = est_superinfection_H1_positive_percent_const';
results_super.est_superinfection_H3_positive_H1_positive_percent_const = est_superinfection_H3_positive_H1_positive_percent_const';

pause;

disp('fitting input-dependent susceptibility model to superinfection data..'); fprintf('\n');

model_num = 1;                   % 0 = constant model; 1 = proportional reduction
susc_reduction_prop = 0.3;       % in proportional model, each virion reduced susceptibility by a certain %
n_parsfit_m1 = length(susc_reduction_prop);

[susc_reduction_prop_fit, results_super.LSE_m1, est_superinfection_H3_positive_percent_fit, est_superinfection_H1_positive_percent_fit, est_superinfection_H3_positive_H1_positive_percent_fit] = GetEstimated_percents_superinfection_fit_fraction(data_super, paramsfit_co, model_num, susc_reduction_prop, params_fixed);

%susc_reduction_prop_fit = susc_reduction_prop;
paramsfit_super(2) = susc_reduction_prop_fit;
results_super.susc_reduction_prop_fit = susc_reduction_prop_fit;
disp(['input-independent susceptibility model estimate: r = ', num2str(susc_reduction_prop_fit,'%2.3f')]); fprintf('\n');

AIC_value_m1 = n_pts*log(results_super.LSE_m1)+2*n_parsfit_m1;

disp('getting 95% confidence intervals...');

% get 95% confidence intervals
results_super.CI_estimate_m1 = get_confidinterval_super(data_super, paramsfit_co, model_num, log(susc_reduction_prop_fit), params_fixed)

% get percent positive curves with best-fit parameters
[est_superinfection_H3_positive_percent_prop, est_superinfection_H1_positive_percent_prop, est_superinfection_H3_positive_H1_positive_percent_prop] = GetEstimated_percents_superinfection_sim_fraction(transpose(results_super.est_actual_moi_rH3N1), paramsfit_co, model_num, paramsfit_super(2), params_fixed);

results_super.est_superinfection_H3_positive_percent_prop = est_superinfection_H3_positive_percent_prop';
results_super.est_superinfection_H1_positive_percent_prop = est_superinfection_H1_positive_percent_prop';
results_super.est_superinfection_H3_positive_H1_positive_percent_prop = est_superinfection_H3_positive_H1_positive_percent_prop';

% calculate deltaAIC
results_super.AIC_values = [AIC_value_m0, AIC_value_m1];
results_super.delta_AIC_values = results_super.AIC_values - min(results_super.AIC_values)


if 1
    
%     save('superinfection_parameter_estimates');
%     save('superinfection_parameter_estimates_update051320');
    save('superinfection_parameter_estimates_update081120');
    
end



