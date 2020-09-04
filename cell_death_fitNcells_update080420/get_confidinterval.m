function confidence_interval = get_confidinterval(mean_log_theta,var_log_theta,z_score)

% standard error based on numerical estimates of Fisher information --
% https://stats.stackexchange.com/questions/238378/how-to-estimate-confidence-interval-of-a-least-squares-fit-parameters-by-means-o

% the scripts 'get_var_nb.m' and 'get_var_poisson.m' numerically estimate the
% Hessian and its inverse to get the asymptotic variance estimate Var(log(theta))

% standard error based on variance estimates
SE_log_pars = sqrt(var_log_theta);

% the Wald confidence interval for theta on log scale -- 
% log(theta) +- 1.96*SE(log(theta))
confidence_interval_ub_log = mean_log_theta  + z_score*SE_log_pars;
confidence_interval_lb_log = mean_log_theta  - z_score*SE_log_pars;

% back transform interval (note: variance estimate already accounts for log transform)
confidence_interval = [exp(confidence_interval_lb_log), exp(confidence_interval_ub_log)];