function cell_deathrate = get_celldeathrate_model(alpha,time_hpi,model_num,i_val)

% to plot cell death rate models vs. hpi under best fit parameter estimates

switch model_num
    
    case 1
        
        % cell death rate: time-independent, input-independent - a
        a = alpha(1);
        
        cell_deathrate = a*ones(size(time_hpi));
        cell_deathrate(1) = 0;
        
    case 2
        
        % cell death rate: time-dependent, input-independent - b*k*hpi^(k-1)
        k = alpha(1); % shape parameter
        b = alpha(2); % scale parameters
        
        cell_deathrate = b*k*time_hpi.^(k-1);
        
    case 3
        
        % cell death rate: time-independent, input-dependent - c*i^eps
        c = alpha(1);
        eps = alpha(2);

        cell_deathrate = c*ones(size(time_hpi))*i_val.^eps;
        cell_deathrate(1) = 0;
        
    case 4
        
        % cell death rate: time-dependent, input-dependent - b*k*hpi^(k-1)*i^eps,
        k = alpha(1);
        d = alpha(2);
        eps = alpha(3);
        
        cell_deathrate = d*k*time_hpi.^(k-1)*i_val.^eps;
        
end