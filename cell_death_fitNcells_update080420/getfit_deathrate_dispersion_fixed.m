function sumsq_error = getfit_deathrate_dispersion_fixed(log_alpha, r_dispersion, moi_vector_3hpi, moi_vector, moi_vector_FACs_data, data_allhpi, FACs_data_18hpi, hpi, model_num, which_distribution)

alpha_vals = exp(log_alpha);

% only fit to the first 4 treatment groups of FACs data
these_moi_FACs = moi_vector_FACs_data;

% switch-case over models: 1-4
switch model_num
    
    case 1
        
        % 'model 1': constant rate of cell death, input-independent
        celldeath_rates = alpha_vals(1:2); % includes mu
        Ncells_init = alpha_vals(3);
%         mu = alpha_dispersion(2); % constant death rate of uninfected cells
%         k_dispersion = alpha_dispersion(end);
        
    case 2
        % 'model 2': Weibull hazard function cell death rate, input-independent
        celldeath_rates = alpha_vals(1:3); % includes mu
        Ncells_init = alpha_vals(4);
%         mu = alpha_dispersion(3); % constant death rate of uninfected cells
%         k_dispersion = alpha_dispersion(end);
        
    case 3
        % 'model 3': constant rate of cell death (time-independent), input-dependent
        celldeath_rates = alpha_vals(1:3); % includes mu
        Ncells_init = alpha_vals(4);
%         mu = alpha_dispersion(3); % constant death rate of uninfected cells
%         k_dispersion = alpha_dispersion(end);
        
    case 4
        % 'model 4': time-dependent + input-dependent
        celldeath_rates = alpha_vals(1:4); % includes mu
        Ncells_init = alpha_vals(5);
%         mu = alpha_dispersion(4); % constant death rate of uninfected cells
%         k_dispersion = alpha_dispersion(end);
        
end


sumsq_error = 0;
for n = 1:length(hpi)
    
    % this hour post infection
    this_hpi = hpi(n);
    
    % data at this hpi
    data_this_hpi = data_allhpi(n,:);
    
    if n == 1    
        these_moi = moi_vector_3hpi;
    else
        these_moi = moi_vector;
    end
    
    
    if n == 4 % 18 hpi
        
        transform_data_this_hpi = 100*data_this_hpi/Ncells_init;
        
        percentcells_alive = get_percentalive(celldeath_rates,this_hpi,these_moi,r_dispersion,model_num,which_distribution);
        sumsq_error = sumsq_error + sum((transform_data_this_hpi - percentcells_alive).^2); % linear
%         sumsq_error = sumsq_error + sum((log10(data_this_hpi) - log10(percentcells_alive)).^2); % log
        
        percent_aliveinfected = get_percentaliveinfected(celldeath_rates,this_hpi,these_moi_FACs,r_dispersion,model_num,which_distribution);
        sumsq_error = sumsq_error + sum((FACs_data_18hpi - percent_aliveinfected).^2); % linear
%         sumsq_error = sumsq_error + sum((log10(FACs_data_18hpi) - log10(percent_aliveinfected)).^2); % log
        
    else
        
        transform_data_this_hpi = 100*data_this_hpi/Ncells_init;
        
        percentcells_alive = get_percentalive(celldeath_rates,this_hpi,these_moi,r_dispersion,model_num,which_distribution);
        sumsq_error = sumsq_error + sum((transform_data_this_hpi - percentcells_alive).^2); % linear
%         sumsq_error = sumsq_error + sum((log10(data_this_hpi) - log10(percentcells_alive)).^2); % log

        
    end
    
    
    
end

