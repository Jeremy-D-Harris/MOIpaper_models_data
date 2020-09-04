function CI_estimate = get_confidinterval_co(data_co, all_params_estimates, params_fixed)

% for coinfection parameter estimates

% only four of the six parameters were fit
paramsfit_estimates = [all_params_estimates(1) all_params_estimates(2) all_params_estimates(5) all_params_estimates(6)];
% params_for_sim = [params_fit_estimates(1) params_fit_estimates(2) params_fit_estimates(2) params_fit_estimates(2) params_fit_estimates(2) exp(xfit(4))];


% updated 05/13/20 using the sandwich method and asymptotic variance
% estimate: see https://stats.stackexchange.com/questions/238378/how-to-estimate-confidence-interval-of-a-least-squares-fit-parameters-by-means-o

coinfection_H3_positive_percent = data_co.coinfection_H3_positive_percent;
coinfection_H1_positive_percent = data_co.coinfection_H1_positive_percent;
coinfection_H3_positive_H1_positive_percent = data_co.coinfection_H3_positive_H1_positive_percent;

% adapted from variance estimate function for virus production -- 05/12/20
z_score = 1.96;

theta = [paramsfit_estimates(1), paramsfit_estimates(3), paramsfit_estimates(4)]; % 1x3
% scale on which parameters were fit
log_theta = log(theta);
% size(theta)

% small value to estimate gradient wrt parameters
h_small = 10^-9;

% regularization parameter for H_inv = inverse(X^T*X+ b_small*Identity)
b_small = 0;

% var_estimate = []; % consistent estimate of the variance in residual error
Jf_theta_transform = []; % Jacobian (rows: # data_co points) x (columns: # parameters)

% [est_H3_positive_percent_vector, est_H1_positive_percent_vector, est_H3_positive_H1_positive_percent_vector] = GetEstimated_percents_coinfection_sim_fraction(moi_rH3N1_vector, params, params_fixed)

% get percent positive curves: best-fit theta
[est_coinfection_H3_positive_percent_bestfit, est_coinfection_H1_positive_percent_bestfit, est_coinfection_H3_positive_H1_positive_percent_bestfit] = GetEstimated_percents_coinfection_sim_fraction(data_co.coinfection_actual_moi_rH3N1, all_params_estimates, params_fixed);

% F = (est_superinfection_H3_positive_percent_bestfit - superinfection_H3_positive_percent).^2 + (est_superinfection_H1_positive_percent_bestfit - superinfection_H1_positive_percent).^2 + (est_superinfection_H3_positive_H1_positive_percent_bestfit - superinfection_H3_positive_H1_positive_percent).^2;
F = transpose([(est_coinfection_H3_positive_percent_bestfit - coinfection_H3_positive_percent).^2, (est_coinfection_H1_positive_percent_bestfit - coinfection_H1_positive_percent).^2 , (est_coinfection_H3_positive_H1_positive_percent_bestfit - coinfection_H3_positive_H1_positive_percent).^2]);

var_estimate = sum(F);

for m = 1:length(theta)
    
    % get percent positive curves: theta+h
    this_theta_plus_h = theta;
    this_theta_plus_h(m) =  theta(m) + h_small;
    
    all_params_estimates_plus_h = [this_theta_plus_h(1) paramsfit_estimates(2) paramsfit_estimates(2) paramsfit_estimates(2) this_theta_plus_h(2) this_theta_plus_h(3)];
    
    [est_coinfection_H3_positive_percent_plus_h, est_coinfection_H1_positive_percent_plus_h, est_coinfection_H3_positive_H1_positive_percent_plus_h] = GetEstimated_percents_coinfection_sim_fraction(data_co.coinfection_actual_moi_rH3N1, all_params_estimates_plus_h, params_fixed);
    
    % F_plus_h = (est_superinfection_H3_positive_percent_plus_h - superinfection_H3_positive_percent).^2 + (est_superinfection_H1_positive_percent_plus_h - superinfection_H1_positive_percent).^2 + (est_superinfection_H3_positive_H1_positive_percent_plus_h - superinfection_H3_positive_H1_positive_percent).^2;
    F_plus_h = transpose([(est_coinfection_H3_positive_percent_plus_h - coinfection_H3_positive_percent).^2, (est_coinfection_H1_positive_percent_plus_h - coinfection_H1_positive_percent).^2, (est_coinfection_H3_positive_H1_positive_percent_plus_h - coinfection_H3_positive_H1_positive_percent).^2]);
    
    % multiply by theta to account for log transform
    Jf_theta_transform(:,m) = theta(m)*(F_plus_h-F)/h_small;
    
end
% n x n matrix, where n is the number of data_co points
sigma_hat = diag(var_estimate);

% (mxn)*(nxn)*(nxm) = mxm, where m is the number of parameters
variance_sandwich = Jf_theta_transform'*sigma_hat*Jf_theta_transform;

% size(variance_sandwich)

% Hessian matrix is mxm
Hessian_inv_estimate = inv(Jf_theta_transform'*Jf_theta_transform + b_small*eye(size(Jf_theta_transform'*Jf_theta_transform)));

C_matrix = Hessian_inv_estimate*variance_sandwich*Hessian_inv_estimate;

standard_error_estimate = sqrt(diag(C_matrix));

% size(standard_error_estimate)

CI_estimate_transform(:,1) = log_theta'-z_score*standard_error_estimate;
CI_estimate_transform(:,2) = log_theta'+z_score*standard_error_estimate;

% bring back to original scale
CI_estimate =  exp(CI_estimate_transform);

