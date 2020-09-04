function confidence_interval = get_confidinterval(mean_log_theta,var_log_theta,z_score)

SE_log_pars = sqrt(var_log_theta);

confidence_interval_ub_log = mean_log_theta  + z_score*SE_log_pars;
confidence_interval_lb_log = mean_log_theta  - z_score*SE_log_pars;

confidence_interval = [exp(confidence_interval_lb_log), exp(confidence_interval_ub_log)];