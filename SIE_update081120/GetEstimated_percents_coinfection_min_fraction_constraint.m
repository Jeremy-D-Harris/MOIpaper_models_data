function f = GetEstimated_percents_coinfection_min_fraction_constraint(x, coinfection_actual_moi_rH3N1, coinfection_H3_positive_percent, coinfection_H1_positive_percent, coinfection_H3_positive_H1_positive_percent, params_fixed)

params = [x(1) (1-exp(x(2))) (1-exp(x(2))) (1-exp(x(2))) exp(x(3)) exp(x(4))];

[est_coinfection_H3_positive_percent_fit, est_coinfection_H1_positive_percent_fit, est_coinfection_H3_positive_H1_positive_percent_fit] = GetEstimated_percents_coinfection_sim_fraction(coinfection_actual_moi_rH3N1, params, params_fixed);

f = sum((est_coinfection_H3_positive_percent_fit - coinfection_H3_positive_percent).^2 + (est_coinfection_H1_positive_percent_fit - coinfection_H1_positive_percent).^2 + (est_coinfection_H3_positive_H1_positive_percent_fit - coinfection_H3_positive_H1_positive_percent).^2);

