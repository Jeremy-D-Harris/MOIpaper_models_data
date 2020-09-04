function sample_par_var_estimate = get_var_timeindependent(log_nu, moi_vector, data_allhpi, all_hpi, alpha, k_dispersion, Ncells, cell_line, model_num)

which_time_model = 0;

nu = exp(log_nu);

% updated on 02/20/20
num_pars = length(nu);

% small value to estimate gradient wrt parameters
h_small = 10^-9;

% regularization parameter for H_inv = inverse(X^T*X+ b_small*Identity)
% b_small = 0;
if model_num == 3
    
    if cell_line == 1
        
        b_small = [1; 1];
        
    else
        
        b_small = [0; 0];
        
    end
    
else
    
    b_small = 0;
    
end


var_estimate = [];
Jf_logtheta = []; % Jacobian (rows: # data points) x (columns: # parameters)

counter = 1;
counter_var = 1;
for n = 1:length(all_hpi)
    
    % this hour post infection
    this_hpi = all_hpi(n);
    
    % data at this hpi
    data_this_hpi = data_allhpi(n,:);
    
    % MOI vals for 6, 12, 18 hpi
    these_moi_vals = moi_vector;
    
    % virus production
    for k = 1:length(these_moi_vals)
        
        this_moi = these_moi_vals(k);
        this_datum = data_this_hpi(1,k);
        
        if isnan(this_datum) ~= 1
            
            log_this_datum = log(this_datum);
            log_GE_output_nu = log(get_GEmodel(nu,this_moi,Ncells,alpha,k_dispersion,this_hpi,cell_line,model_num,which_time_model));
            
            for m = 1:num_pars
                
                nu_plus_h = nu;
                nu_plus_h(m) =  nu(m) + h_small;
                
                log_GE_output_nu_plus_h = log(get_GEmodel(nu_plus_h,this_moi,Ncells,alpha,k_dispersion,this_hpi,cell_line,model_num,which_time_model));
                
                
                Jf_logtheta(counter,m) = (nu(m)/h_small)*(log_GE_output_nu_plus_h - log_GE_output_nu);
                
            end
            
            var_estimate(counter_var,1) = (log_this_datum - log_GE_output_nu)^2;
            counter_var = counter_var+1;
            
            counter=counter+1;
        end
    end
    
end


sigma_hat = diag(var_estimate);

variance_sandwich = Jf_logtheta'*sigma_hat*Jf_logtheta;

Hessian_inv_estimate = inv(Jf_logtheta'*Jf_logtheta + diag(b_small)*eye(size(Jf_logtheta'*Jf_logtheta)));

C_matrix = Hessian_inv_estimate*variance_sandwich*Hessian_inv_estimate;

sample_par_var_estimate = diag(C_matrix);

