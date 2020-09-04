function GE_output = get_GEmodel(nu,moi_vector,Ncells,alpha,k_dispersion,this_hpi,cell_line,model_num, which_time_model)

matlab_r = k_dispersion;

dt = 0.01;
time = 0:dt:this_hpi;
i_vals = 0:1000;


% cell death rate: now have the same model with most support!! -- jdh 08/11/20
% 1:MDCK (time-independent, input-independent)
% 2:A549 (time-independent, input-independent)
switch cell_line
    
    case 1 % MDCK
        
        % cell death rate: time-dependent, input-independent: k,b,mu
        k = alpha(1); % shape parameter
        b = alpha(2); % scale parameters
        mu = alpha(3);
        int_psi = b*time.^k + mu*time; % integrating Weibull hazard function
        
    case 2 % A549
        
        % cell death rate: time-dependent, input-independent: k,b,mu
        k = alpha(1); % shape parameter
        b = alpha(2); % scale parameters
        mu = alpha(3);
        int_psi = b*time.^k + mu*time; % integrating Weibull hazard function
        
end

% probability that an infected cell is still alive
prob_cellsremain_t = exp(-int_psi);


% need to include time-independent models
switch which_time_model
    
    case 0 % time-independent
        
        % which virus production model
        switch model_num
            
            case 1
                
                % virus production rate: input-independent - l*max(t-s,0)
                l = nu(1);
                
                % total cellular output as a function of i WT at time t
                nu_i_vector_t = l;
                viraloutput_giveninfected(1:length(i_vals(2:end))) = sum(nu_i_vector_t.*prob_cellsremain_t)*dt;
                
            case 2
                
                % virus production rate: linear input-dependent - r*max(t-s,0)*i
                r = nu(1);
                counter=1;
                for n = 2:length(i_vals)
                    
                    this_i = i_vals(n);
                    
                    % total cellular output as a function of i WT at time t
                    nu_i_vector_t = r*this_i;
                    viraloutput_giveninfected(1,counter) = sum(nu_i_vector_t.*prob_cellsremain_t)*dt;
                    counter=counter+1;
                end
                
            case 3
                
                % virus production rate: saturating input-dependent - m*max(t-s,0)*i/(K+i)
                m = nu(1);
                K_MM = nu(2);
                
                counter=1;
                for n = 2:length(i_vals)
                    
                    this_i = i_vals(n);
                    
                    % total cellular output as a function of i WT at time t
                    nu_i_vector_t = m*this_i/(this_i+K_MM);
                    viraloutput_giveninfected(1,counter) = sum(nu_i_vector_t.*prob_cellsremain_t)*dt;
                    counter=counter+1;
                end
                
        end
        
        
    case 1 % time-dependent models
        
        % which virus production model
        switch model_num
            
            case 1
                
                % virus production rate: input-independent - l*max(t-s,0)
                l = nu(1);
                s = nu(2);
                
                % total cellular output as a function of i WT at time t
                nu_i_vector_t = l*max((time-s),0);
                viraloutput_giveninfected(1:length(i_vals(2:end))) = sum(nu_i_vector_t.*prob_cellsremain_t)*dt;   
                
            case 2
                
                % virus production rate: linear input-dependent - r*max(t-s,0)*i
                r = nu(1);
                s = nu(2);
                counter=1;
                for n = 2:length(i_vals)
                    
                    this_i = i_vals(n);
                    
                    % total cellular output as a function of i WT at time t
                    nu_i_vector_t = r*(max((time-s),0))*this_i;
                    viraloutput_giveninfected(1,counter) = sum(nu_i_vector_t.*prob_cellsremain_t)*dt;
                    counter=counter+1;
                end
                
            case 3
                
                % virus production rate: saturating input-dependent - m*max(t-s,0)*i/(K+i)
                m = nu(1);
                s = nu(2);
                K_MM = nu(3);
                
                counter=1;
                for n = 2:length(i_vals)
                    
                    this_i = i_vals(n);
                    
                    % total cellular output as a function of i WT at time t
                    nu_i_vector_t = m*(max((time-s),0))*(this_i/(this_i+K_MM));
                    viraloutput_giveninfected(1,counter) = sum(nu_i_vector_t.*prob_cellsremain_t)*dt;
                    counter=counter+1;
                end
                
        end
        
end

for n=1:length(moi_vector)
    
    this_moi = moi_vector(n);
    
    
    matlab_p = 1-this_moi/(matlab_r + this_moi); % formulation of negative binomial in matlab speak
    
    % probability cell not infected
    init_prob_moi0 = nbinpdf(i_vals(1),matlab_r,matlab_p);
    
    % probability of i virions entering - according to negative binomial
    init_prob_moi1plus = nbinpdf(i_vals(2:end),matlab_r,matlab_p);
    
    
    init_prob_moi1plus(end) = init_prob_moi1plus(end) + (1 - (sum(init_prob_moi1plus)+init_prob_moi0));
    
    % total virus production
    GE_output(n) = Ncells*sum(init_prob_moi1plus.*viraloutput_giveninfected);
    
    
end

