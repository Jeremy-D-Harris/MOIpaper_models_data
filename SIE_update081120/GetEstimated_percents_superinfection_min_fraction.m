function f = GetEstimated_percents_superinfection_min_fraction(x, superinfection_actual_moi_rH3N1, superinfection_H3_positive_percent, superinfection_H1_positive_percent, superinfection_H3_positive_H1_positive_percent, paramsfit_co, model_num, params_fixed)

susc_reduction = exp(x(1));
 
[est_superinfection_H3_positive_percent_fit, est_superinfection_H1_positive_percent_fit, est_superinfection_H3_positive_H1_positive_percent_fit] = GetEstimated_percents_superinfection_sim_fraction(superinfection_actual_moi_rH3N1, paramsfit_co, model_num, susc_reduction, params_fixed);

% minimize sum squared error simultaneously
f = sum((est_superinfection_H3_positive_percent_fit - superinfection_H3_positive_percent).^2 + (est_superinfection_H1_positive_percent_fit - superinfection_H1_positive_percent).^2 + (est_superinfection_H3_positive_H1_positive_percent_fit - superinfection_H3_positive_H1_positive_percent).^2);

% disp(['susceptible reduction estimate = ', num2str(susc_reduction,'%2.2f')]);
% disp(['sum squared error = ', num2str(f,'%2.0f')]);
