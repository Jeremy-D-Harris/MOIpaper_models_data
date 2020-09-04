function CI_estimate = get_confidinterval_super(data_super, paramsfit_co, model_num, log_susc_reduction, params_fixed)

% for superinfection parameter estimates

% updated 05/13/20 using the sandwich method and asymptotic variance
% estimate: see https://stats.stackexchange.com/questions/238378/how-to-estimate-confidence-interval-of-a-least-squares-fit-parameters-by-means-o

moi_rH3N1_vector = data_super.superinfection_actual_moi_rH3N1;

superinfection_H3_positive_percent = data_super.superinfection_H3_positive_percent;
superinfection_H1_positive_percent = data_super.superinfection_H1_positive_percent;
superinfection_H3_positive_H1_positive_percent = data_super.superinfection_H3_positive_H1_positive_percent;

% adapted from variance estimate function for virus production -- 05/12/20
z_score = 1.96;

log_theta = log_susc_reduction;
theta = exp(log_theta);

% small value to estimate gradient wrt parameters
h_small = 10^-9;

% regularization parameter for H_inv = inverse(X^T*X+ b_small*Identity)
b_small = 0;

% var_estimate = []; % consistent estimate of the variance in residual error
Jf_logtheta = []; % Jacobian (rows: # data points) x (columns: # parameters)

% get percent positive curves: best-fit theta
[est_superinfection_H3_positive_percent_bestfit, est_superinfection_H1_positive_percent_bestfit, est_superinfection_H3_positive_H1_positive_percent_bestfit] = GetEstimated_percents_superinfection_sim_fraction(moi_rH3N1_vector, paramsfit_co, model_num, theta, params_fixed);

% F = (est_superinfection_H3_positive_percent_bestfit - superinfection_H3_positive_percent).^2 + (est_superinfection_H1_positive_percent_bestfit - superinfection_H1_positive_percent).^2 + (est_superinfection_H3_positive_H1_positive_percent_bestfit - superinfection_H3_positive_H1_positive_percent).^2;
F = transpose([(est_superinfection_H3_positive_percent_bestfit - superinfection_H3_positive_percent).^2, (est_superinfection_H1_positive_percent_bestfit - superinfection_H1_positive_percent).^2 , (est_superinfection_H3_positive_H1_positive_percent_bestfit - superinfection_H3_positive_H1_positive_percent).^2]);

var_estimate = sum(F);

for m = 1:length(theta)

% get percent positive curves: theta+h
this_theta_plus_h =  theta(m) + h_small;
[est_superinfection_H3_positive_percent_plus_h, est_superinfection_H1_positive_percent_plus_h, est_superinfection_H3_positive_H1_positive_percent_plus_h] = GetEstimated_percents_superinfection_sim_fraction(moi_rH3N1_vector, paramsfit_co, model_num, this_theta_plus_h, params_fixed);

% F_plus_h = (est_superinfection_H3_positive_percent_plus_h - superinfection_H3_positive_percent).^2 + (est_superinfection_H1_positive_percent_plus_h - superinfection_H1_positive_percent).^2 + (est_superinfection_H3_positive_H1_positive_percent_plus_h - superinfection_H3_positive_H1_positive_percent).^2;
F_plus_h = transpose([(est_superinfection_H3_positive_percent_plus_h - superinfection_H3_positive_percent).^2, (est_superinfection_H1_positive_percent_plus_h - superinfection_H1_positive_percent).^2, (est_superinfection_H3_positive_H1_positive_percent_plus_h - superinfection_H3_positive_H1_positive_percent).^2]);

% multiply by theta to account for log transform
Jf_logtheta(:,m) = theta(m)*(F_plus_h-F)/h_small;

end


sigma_hat = diag(var_estimate);

variance_sandwich = Jf_logtheta'*sigma_hat*Jf_logtheta;

% size(variance_sandwich)

Hessian_inv_estimate = inv(Jf_logtheta'*Jf_logtheta + b_small*eye(size(Jf_logtheta'*Jf_logtheta)));

C_matrix = Hessian_inv_estimate*variance_sandwich*Hessian_inv_estimate;

standard_error_estimate = transpose(sqrt(diag(C_matrix)));

CI_estimate(1) = exp((log_theta-z_score*standard_error_estimate));
CI_estimate(2) = exp((log_theta+z_score*standard_error_estimate));

