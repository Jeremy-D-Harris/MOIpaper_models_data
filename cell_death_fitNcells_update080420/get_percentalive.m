function percent_alive = get_percentalive(alpha, hpi, moi, k_dispersion, model_num, which_distribution)


i_vals = 0:1000;

switch model_num
    
    case 1
        
        % cell death rate: time-independent, input-independent - a
        a = alpha(1);
        mu = alpha(2);
        int_psi = a*hpi; % integrating
        
    case 2
        
        % cell death rate: time-dependent, input-independent - b*k*hpi^(k-1)
        k = alpha(1); % shape parameter
        b = alpha(2); % scale parameters
        mu = alpha(3);
        int_psi = b*hpi^k; % integrating Weibull hazard function
        
    case 3
        
        % cell death rate: time-independent, input-dependent - c*i^eps
        c = alpha(1);
        eps = alpha(2);
        mu = alpha(3);
        int_psi = c*hpi*i_vals(2:end).^eps; % time integration
        
    case 4
        
        % cell death rate: time-dependent, input-dependent - d*k*hpi^(k-1)*i^eps,
        k = alpha(1);
        d = alpha(2);
        eps = alpha(3);
        mu = alpha(4);
        int_psi = d*hpi^k*i_vals(2:end).^eps; % time integration
        
end

for n=1:length(moi)
    
    this_moi = moi(n);
    
    switch which_distribution
        
        case 1 % Poisson
            
            % probability cell not infected
            init_prob_moi0 = poisspdf(i_vals(1), this_moi);
            
            % probability of i virions entering - according to Poisson
            init_prob_moi1plus = poisspdf(i_vals(2:end), this_moi);
            
        case 2 % negative binomial
            
            p_extrazeros = k_dispersion;
            
            matlab_p = 1-this_moi/(p_extrazeros + this_moi); % formulation of negative binomial in matlab speak
            
            % probability cell not infected
            init_prob_moi0 = nbinpdf(i_vals(1),p_extrazeros,matlab_p);
            
            % probability of i virions entering - according to negative binomial
            init_prob_moi1plus = nbinpdf(i_vals(2:end),p_extrazeros,matlab_p);
            
        case 3 % zero-inflated Poisson
            
            % based on wikipedia parametrization:
            %             https://en.wikipedia.org/wiki/Zero-inflated_model
            
            p_extrazeros = k_dispersion;
            
            lambda_val = this_moi/(1-p_extrazeros);
            
            % probability cell not infected
            init_prob_moi0 = p_extrazeros + (1-p_extrazeros)*poisspdf(i_vals(1), lambda_val);
            
            % probability of i virions entering - according to Poisson
            init_prob_moi1plus = (1-p_extrazeros)*poisspdf(i_vals(2:end), lambda_val);
            
            
            %             p_extrazeros = k_dispersion;
            %
            %             lambda_val = this_moi;%this_moi/(1-p_extrazeros);
            % %             matlab_p = 1-this_moi/(p_extrazeros + this_moi); % formulation of negative binomial in matlab speak
            %
            %             % probability cell not infected
            %             init_prob_moi0 = p_extrazeros + (1-p_extrazeros)*exp(-lambda_val); %nbinpdf(i_vals(1),p_extrazeros,matlab_p);
            %
            %             % probability of i virions entering - according to zero-inflated Poisson
            %             init_prob_moi1plus = (1-p_extrazeros)*exp(-lambda_val)*((lambda_val.^i_vals(2:end))./factorial(i_vals(2:end))); %nbinpdf(i_vals(2:end),p_extrazeros,matlab_p);
            
    end
    
    init_prob_moi1plus(end) = init_prob_moi1plus(end) + (1 - (sum(init_prob_moi1plus)+init_prob_moi0));
    
    % int_psi: from integrating cell death rate model (above)
    % also add the background death rate mu: integration gives mu*hpi
    prob_infectedcells_alive = exp(-mu*hpi)*sum(init_prob_moi1plus.*exp(-int_psi));
    
    % some background death of uninfected cells, mu: integration gives mu*hpi
    prob_uninfectedcells_alive = exp(-mu*hpi)*init_prob_moi0;
    
    % total prob alive
    prob_alive = prob_uninfectedcells_alive + prob_infectedcells_alive;
    
    % assume that cells only die if they are infected
    percent_alive(n) = 100*prob_alive;
    
end
