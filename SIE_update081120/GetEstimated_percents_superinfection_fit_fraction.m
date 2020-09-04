function [susc_reduction_fit, fval, est_superinfection_H3_positive_percent_fit, est_superinfection_H1_positive_percent_fit, est_superinfection_H3_positive_H1_positive_percent_fit] = GetEstimated_percents_superinfection_fit_fraction(data_super, paramsfit_co, model_num, susc_reduction, params_fixed)

x0 = log(susc_reduction); % transform variables to log scale
options = optimset('Display','iter','MaxIter',200,'TolFun',1e-1, 'TolX', 1e-3);
% options = optimset('Display','iter','MaxIter',100,'TolFun',1e-1);
[xfit, fval, exitflag, output] = fminsearch(@(x)GetEstimated_percents_superinfection_min_fraction(x, data_super.superinfection_actual_moi_rH3N1, data_super.superinfection_H3_positive_percent, data_super.superinfection_H1_positive_percent, data_super.superinfection_H3_positive_H1_positive_percent, paramsfit_co, model_num, params_fixed), x0, options); % this function searches for the value that minimimizes the sum of squared error between data and model prediction
susc_reduction_fit = exp(xfit);

[est_superinfection_H3_positive_percent_fit, est_superinfection_H1_positive_percent_fit, est_superinfection_H3_positive_H1_positive_percent_fit] = GetEstimated_percents_superinfection_sim_fraction(data_super.superinfection_actual_moi_rH3N1, paramsfit_co, model_num, susc_reduction_fit, params_fixed);
