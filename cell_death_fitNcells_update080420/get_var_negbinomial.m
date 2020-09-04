function sample_param_var_estimate = get_var_negbinomial(log_alpha, moi_vector_3hpi, moi_vector, moi_vector_FACs_data, all_hpi, data_allhpi, FACs_data_18hpi, model_num, cell_line)

% variance estimates based on numerical estimates of Fisher information --
% https://stats.stackexchange.com/questions/238378/how-to-estimate-confidence-interval-of-a-least-squares-fit-parameters-by-means-o

which_distribution = 2;

alpha = exp(log_alpha);
Ncells_init = alpha(end-1);

these_moi_vals_FACs = moi_vector_FACs_data;

% updated on 02/20/20
num_pars = length(alpha);
% n_pts_data = n_pts_allphi_FACs; % number of rows of Jacobian
% n_degrees =  n_pts_data - num_pars; % if using student's t-test, or unbiased estimator of variance


% for 95 percent confidence intervals
% significance_level = 0.05;
% p = 1-significance_level/2;
% prefactor = 1.96; %tinv(p,n_degrees); % based on z-score 95% CI or t-test(p,n_degrees);

% small value to estimate gradient wrt parameters
h_small = 10^-6;

% regularization parameter for H_inv = inverse(X^T*X+ b_small*Identity)
if model_num == 3
    
    b_small = 1e-10;
    
elseif model_num == 4
    
    if cell_line == 1
        b_small = 1e-3;
    else
        b_small = 1e-10;
    end
    
else
    
    b_small = 0;
    
end

var_estimate = []; % consistent estimate of the variance in residual error
Jf_logtheta = []; % Jacobian (rows: # data points) x (columns: # parameters)

counter = 1;
counter_var = 1;
for n = 1:length(all_hpi)
    
    % this hour post infection
    this_hpi = all_hpi(n);
    
    % data at this hpi
    data_this_hpi = data_allhpi(n,:);
    
    % different MOI vals for 3 vs. 6, 12, 18 hpi
    if n == 1
        
        these_moi_vals = moi_vector_3hpi;
        
    else
        
        these_moi_vals = moi_vector;
        
    end
    
    % at 18 hpi, include flow data also
    if n ==4
        
        
        % percent cells remaining - trypan blue staining
        for k = 1:length(these_moi_vals)
            
            this_moi = these_moi_vals(k);
            this_datum = data_this_hpi(1,k);
            
            
            for m = 1:num_pars
                
                alpha_plus_h = alpha;
                alpha_plus_h(m) =  alpha(m) + h_small;
                
                if m == (num_pars-1) % different for N cells init
                % partial derivative with chain rule to account for the transformation log(alpha)
                Jf_logtheta(counter,m) = -100*this_datum/Ncells_init; % have to multiply by alpha(m)
                
                else 
                % partial derivative with chain rule to account for the transformation log(alpha)
                Jf_logtheta(counter,m) = (alpha(m)/h_small)*(get_percentalive(alpha_plus_h, this_hpi, this_moi, alpha_plus_h(end), model_num, which_distribution) - get_percentalive(alpha, this_hpi, this_moi, alpha(end), model_num, which_distribution));
                end
                
            end
            
            transform_this_datum = 100*this_datum/Ncells_init;
            var_estimate(1,counter_var) = (transform_this_datum - get_percentalive(alpha, this_hpi, this_moi, alpha(end), model_num, which_distribution))^2;
            counter_var = counter_var+1;
            
            counter=counter+1;
            
        end
        
        
        % percent cells infected - flow data
        for k = 1:length(these_moi_vals_FACs)
            
            this_moi = these_moi_vals_FACs(k);
            this_datum = FACs_data_18hpi(1,k);
            
            
            for m = 1:num_pars
                
                alpha_plus_h = alpha;
                alpha_plus_h(m) =  alpha(m) + h_small;
                
                % partial derivative with chain rule to account for the transformation log(alpha)
                Jf_logtheta(counter,m) = (alpha(m)/h_small)*(get_percentaliveinfected(alpha_plus_h, this_hpi, this_moi, alpha_plus_h(end), model_num, which_distribution) - get_percentaliveinfected(alpha, this_hpi, this_moi, alpha(end), model_num, which_distribution));
                
            end
            
            
            var_estimate(1,counter_var) = (this_datum - get_percentaliveinfected(alpha, this_hpi, this_moi, alpha(end), model_num, which_distribution))^2;
            counter_var = counter_var+1;
            
            counter=counter+1;
            
        end
        
        
    % 3, 6, 12 hpi data
    else
        
        
        % percent cells remaining - trypan blue staining
        for k = 1:length(these_moi_vals)
            
            this_moi = these_moi_vals(k);
            this_datum = data_this_hpi(1,k);
            
            
            for m = 1:num_pars
                
                alpha_plus_h = alpha;
                alpha_plus_h(m) =  alpha(m) + h_small;
                
                if m == (num_pars-1) % different for N cells init
                % partial derivative with chain rule to account for the transformation log(alpha)
                Jf_logtheta(counter,m) = -100*this_datum/Ncells_init; % have to multiply by alpha(m)
                
                else 
                % partial derivative with chain rule to account for the transformation log(alpha)
                Jf_logtheta(counter,m) = (alpha(m)/h_small)*(get_percentalive(alpha_plus_h, this_hpi, this_moi, alpha_plus_h(end), model_num, which_distribution) - get_percentalive(alpha, this_hpi, this_moi, alpha(end), model_num, which_distribution));
                end
                
            end
            
            transform_this_datum = 100*this_datum/Ncells_init;
            var_estimate(1,counter_var) = (transform_this_datum - get_percentalive(alpha, this_hpi, this_moi, alpha(end), model_num, which_distribution))^2;
            counter_var = counter_var+1;
            
            
            counter=counter+1;
            
        end
        
    end
    
end


sigma_hat = diag(var_estimate);

variance_sandwich = Jf_logtheta'*sigma_hat*Jf_logtheta;

Hessian_inv_estimate = inv(Jf_logtheta'*Jf_logtheta + b_small*eye(size(Jf_logtheta'*Jf_logtheta)));

C_matrix = Hessian_inv_estimate*variance_sandwich*Hessian_inv_estimate;


sample_param_var_estimate = diag(C_matrix);

