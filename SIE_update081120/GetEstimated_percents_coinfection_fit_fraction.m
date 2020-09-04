function [params_fit, est_coinfection_H3_positive_percent_fit, est_coinfection_H1_positive_percent_fit, est_coinfection_H3_positive_H1_positive_percent_fit] = GetEstimated_percents_coinfection_fit_fraction(data, params_init, params_fixed)

x0 = [params_init(1) log(1-params_init(2)) log(params_init(5)) log(params_init(6))]; % force pp to be the same for PB2/PB1/PAB and both H3 and H1
options = optimset('Display','iter', 'MaxIter',200,'TolFun',1e-1);
[xfit, fval, exitflag, output] = fminsearch(@(x)GetEstimated_percents_coinfection_min_fraction_constraint(x, data.coinfection_actual_moi_rH3N1, data.coinfection_H3_positive_percent, data.coinfection_H1_positive_percent, data.coinfection_H3_positive_H1_positive_percent, params_fixed), x0, options); % this function searches for the value that minimimizes the sum of squared error between data and model prediction

% only four unique parameters being fit, since xfit(2) = prob(PB1) = prob(PB2) = prob(PA)  
params_fit = [xfit(1) (1-exp(xfit(2))) (1-exp(xfit(2))) (1-exp(xfit(2))) exp(xfit(3)) exp(xfit(4))]; 

[est_coinfection_H3_positive_percent_fit, est_coinfection_H1_positive_percent_fit, est_coinfection_H3_positive_H1_positive_percent_fit] = GetEstimated_percents_coinfection_sim_fraction(data.coinfection_actual_moi_rH3N1, params_fit, params_fixed);
