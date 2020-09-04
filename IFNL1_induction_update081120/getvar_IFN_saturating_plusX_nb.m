function sample_par_var_estimate = getvar_IFN_saturating_plusX_nb(log_lambda,moi_vector,log2_IFNdata,params_inputindependent,params_linear,params_saturating,alpha,hpi,model_ind,k_dispersion)

lambda = exp(log_lambda);

num_pars = length(lambda);

% small value to estimate gradient wrt parameters
h_small = 10^-9;

% regularization parameter for H_inv = inverse(X^T*X+ b_small*Identity)
b_small = 10^-10;

var_estimate = [];
Jf_logtheta = []; % Jacobian (rows: # data points) x (columns: # parameters)

counter = 1;
counter_var = 1;
for k = 1:length(moi_vector)
    
    this_moi = moi_vector(k);
    
    this_8hpi_data = log2_IFNdata(1,k);
    this_18hpi_data = log2_IFNdata(2,k);
    
    
    for n = 1:length(hpi)
        
        if n==1
            % saturating (0-8 hpi)
            
            if isnan(this_8hpi_data)~=1
                
                
                params_saturating.r = lambda(1);
                params_saturating.K = lambda(2);
                IFNmodel_lambda_foldchange_saturating_8hpi = get_IFN_saturating_nb(this_moi,params_saturating,alpha,hpi(1),k_dispersion);
                log2_lambda_IFNmodel_8hpi = log2(IFNmodel_lambda_foldchange_saturating_8hpi);
                
                for m = 1:2

                lambda_plus_h = lambda;
                lambda_plus_h(m) = lambda(m) + h_small;

                params_saturating.m = lambda_plus_h(1);
                params_saturating.K = lambda_plus_h(2);

                IFNmodel_lambda_plus_h_foldchange_saturating_8hpi = get_IFN_saturating_nb(this_moi,params_saturating,alpha,hpi(1),k_dispersion);
                log2_lambda_plus_h_IFNmodel_8hpi = log2(IFNmodel_lambda_plus_h_foldchange_saturating_8hpi);
                
                Jf_logtheta(counter,m) = (lambda(m)/h_small)*(log2_lambda_plus_h_IFNmodel_8hpi - log2_lambda_IFNmodel_8hpi);
                
                
                end
                Jf_logtheta(counter,3:num_pars) = 0;
                var_estimate(counter_var,1) = (this_8hpi_data - log2_lambda_IFNmodel_8hpi)^2;

                counter = counter + 1;
                counter_var = counter_var + 1;
            end
            
        else
            
            if isnan(this_18hpi_data)~=1
                
                if model_ind == 1
                    
                    % input-independent (8-18 hpi)
                    params_inputindependent.r = lambda(3);
                    IFNmodel_lambda_foldchange_inputindependent_0to8hpi = get_IFN_inputindependent_nb(this_moi,params_inputindependent,alpha,hpi(1),k_dispersion);
                    IFNmodel_lambda_foldchange_inputindependent_0to18hpi = get_IFN_inputindependent_nb(this_moi,params_inputindependent,alpha,hpi(2),k_dispersion);
                    IFNmodel_lambda_foldchange_inputindependent_8to18hpi = IFNmodel_lambda_foldchange_inputindependent_0to18hpi-IFNmodel_lambda_foldchange_inputindependent_0to8hpi;
                    
                    log2_lambda_IFNmodel_saturatingplusX_18hpi = log2(IFNmodel_lambda_foldchange_saturating_8hpi+IFNmodel_lambda_foldchange_inputindependent_8to18hpi);
                    
                    lambda_plus_h = lambda;
                    lambda_plus_h(3) = lambda(3)+h_small;
                    params_inputindependent.r = lambda_plus_h(3);
                    IFNmodel_lambda_plus_h_foldchange_inputindependent_0to8hpi = get_IFN_inputindependent_nb(this_moi,params_inputindependent,alpha,hpi(1),k_dispersion);
                    IFNmodel_lambda_plus_h_foldchange_inputindependent_0to18hpi = get_IFN_inputindependent_nb(this_moi,params_inputindependent,alpha,hpi(2),k_dispersion);
                    IFNmodel_lambda_plus_h_foldchange_inputindependent_8to18hpi = IFNmodel_lambda_plus_h_foldchange_inputindependent_0to18hpi-IFNmodel_lambda_plus_h_foldchange_inputindependent_0to8hpi;
                    
                    log2_lambda_plus_h_IFNmodel_saturatingplusX_18hpi = log2(IFNmodel_lambda_foldchange_saturating_8hpi+IFNmodel_lambda_plus_h_foldchange_inputindependent_8to18hpi);
                    
                    Jf_logtheta(counter,1:2) = 0;
                    Jf_logtheta(counter,3) = (lambda(3)/h_small)*(log2_lambda_plus_h_IFNmodel_saturatingplusX_18hpi - log2_lambda_IFNmodel_saturatingplusX_18hpi);
                    
                    var_estimate(counter_var,1) = (this_18hpi_data - log2_lambda_IFNmodel_saturatingplusX_18hpi)^2;
                    
                    counter = counter + 1;
                    counter_var = counter_var + 1;
                    
                elseif model_ind == 2
                    % linear (8-18 hpi)
                    
                    params_linear.r = lambda(3);
                    IFNmodel_lambda_foldchange_linear_0to8hpi = get_IFN_linear_nb(this_moi,params_linear,alpha,hpi(1),k_dispersion);
                    IFNmodel_lambda_foldchange_linear_0to18hpi = get_IFN_linear_nb(this_moi,params_linear,alpha,hpi(2),k_dispersion);
                    IFNmodel_lambda_foldchange_linear_8to18hpi = IFNmodel_lambda_foldchange_linear_0to18hpi-IFNmodel_lambda_foldchange_linear_0to8hpi;
                    
                    log2_lambda_IFNmodel_saturatingplusX_18hpi = log2(IFNmodel_lambda_foldchange_saturating_8hpi+IFNmodel_lambda_foldchange_linear_8to18hpi);
                    
                    lambda_plus_h = lambda;
                    lambda_plus_h(3) = lambda(3)+h_small;
                    params_linear.r = lambda_plus_h(3);
                    IFNmodel_lambda_plus_h_foldchange_linear_0to8hpi = get_IFN_linear_nb(this_moi,params_linear,alpha,hpi(1),k_dispersion);
                    IFNmodel_lambda_plus_h_foldchange_linear_0to18hpi = get_IFN_linear_nb(this_moi,params_linear,alpha,hpi(2),k_dispersion);
                    IFNmodel_lambda_plus_h_foldchange_linear_8to18hpi = IFNmodel_lambda_plus_h_foldchange_linear_0to18hpi-IFNmodel_lambda_plus_h_foldchange_linear_0to8hpi;
                    
                    log2_lambda_plus_h_IFNmodel_saturatingplusX_18hpi = log2(IFNmodel_lambda_foldchange_saturating_8hpi+IFNmodel_lambda_plus_h_foldchange_linear_8to18hpi);
                    
                    Jf_logtheta(counter,1:2) = 0;
                    Jf_logtheta(counter,3) = (lambda(3)/h_small)*(log2_lambda_plus_h_IFNmodel_saturatingplusX_18hpi - log2_lambda_IFNmodel_saturatingplusX_18hpi);
                    
                    var_estimate(counter_var,1) = (this_18hpi_data - log2_lambda_IFNmodel_saturatingplusX_18hpi)^2;
                    
                    counter = counter + 1;
                    counter_var = counter_var + 1;
                    
                else
                    % saturating (8-18 hpi)
                    params_saturating.r = lambda(3);
                    params_saturating.K = lambda(4);
                    
                    IFNmodel_lambda_foldchange_saturating_0to8hpi = get_IFN_saturating_nb(this_moi, params_saturating,alpha,hpi(1),k_dispersion);
                    IFNmodel_lambda_foldchange_saturating_0to18hpi = get_IFN_saturating_nb(this_moi, params_saturating,alpha,hpi(2),k_dispersion);
                    IFNmodel_lambda_foldchange_saturating_8to18hpi = IFNmodel_lambda_foldchange_saturating_0to18hpi-IFNmodel_lambda_foldchange_saturating_0to8hpi;
                    
                    log2_lambda_IFNmodel_saturatingplusX_18hpi = log2(IFNmodel_lambda_foldchange_saturating_8hpi+IFNmodel_lambda_foldchange_saturating_8to18hpi);
                    
                    for m = 3:num_pars
                        
                        lambda_plus_h = lambda;
                        lambda_plus_h(m) =  lambda(m) + h_small;
                        
                        params_saturating.r = lambda_plus_h(3);
                        params_saturating.K = lambda_plus_h(4);
                        
                        IFNmodel_lambda_plus_h_foldchange_saturating_0to8hpi = get_IFN_saturating_nb(this_moi, params_saturating,alpha,hpi(1),k_dispersion);
                        IFNmodel_lambda_plus_h_foldchange_saturating_0to18hpi = get_IFN_saturating_nb(this_moi, params_saturating,alpha,hpi(2),k_dispersion);
                        IFNmodel_lambda_plus_h_foldchange_saturating_8to18hpi = IFNmodel_lambda_plus_h_foldchange_saturating_0to18hpi-IFNmodel_lambda_plus_h_foldchange_saturating_0to8hpi;
                        
                        log2_lambda_plus_h_IFNmodel_saturatingplusX_18hpi = log2(IFNmodel_lambda_foldchange_saturating_8hpi+IFNmodel_lambda_plus_h_foldchange_saturating_8to18hpi);
                        
                        Jf_logtheta(counter,m) = (lambda(m)/h_small)*(log2_lambda_plus_h_IFNmodel_saturatingplusX_18hpi - log2_lambda_IFNmodel_saturatingplusX_18hpi);
                        
                    end
                    
                    Jf_logtheta(counter,1:2) = 0;
                    
                    var_estimate(counter_var,1) = (this_18hpi_data - log2_lambda_plus_h_IFNmodel_saturatingplusX_18hpi)^2;
                    
                    counter = counter + 1;
                    counter_var = counter_var + 1;
                    
                end
                
            end
            
        end
    end
    
end



sigma_hat = diag(var_estimate);

variance_sandwich = Jf_logtheta'*sigma_hat*Jf_logtheta;

Hessian_inv_estimate = inv(Jf_logtheta'*Jf_logtheta + diag(b_small)*eye(size(Jf_logtheta'*Jf_logtheta)));

C_matrix = Hessian_inv_estimate*variance_sandwich*Hessian_inv_estimate;

sample_par_var_estimate = diag(C_matrix);