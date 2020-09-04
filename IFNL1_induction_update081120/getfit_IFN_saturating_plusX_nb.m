function f = getfit_IFN_saturating_plusX_nb(log_lambda,moi,log2_IFNdata,params_inputindependent,params_linear,params_saturating,alpha,hpi,model_ind,k_dispersion)

lambda = exp(log_lambda);

n_pts = length(log2_IFNdata(1,:));

sum_sqerror = 0;
for n = 1:n_pts
    
    this_moi = moi(1,n);
    this_8hpi_data = log2_IFNdata(1,n);
    this_18hpi_data = log2_IFNdata(2,n);
    
    if isnan(this_8hpi_data)~=1

            % saturating (0-8 hpi)               
            params_saturating.r = lambda(1);
            params_saturating.K = lambda(2);
            IFNmodel_foldchange_saturating_8hpi = get_IFN_saturating_nb(this_moi,params_saturating,alpha,hpi(1),k_dispersion);
            log2_IFNmodel_foldchange_saturating_8hpi = log2(IFNmodel_foldchange_saturating_8hpi);
            
            sum_sqerror = sum_sqerror + (this_8hpi_data - log2_IFNmodel_foldchange_saturating_8hpi)^2;


        if model_ind == 1
            
            % input-independent (8-18 hpi)
            params_inputindependent.r = lambda(3);
            IFNmodel_foldchange_constant_0to8hpi = get_IFN_inputindependent_nb(this_moi, params_inputindependent,alpha,hpi(1),k_dispersion);
            IFNmodel_foldchange_constant_0to18hpi = get_IFN_inputindependent_nb(this_moi, params_inputindependent,alpha,hpi(2),k_dispersion);
            IFNmodel_foldchange_constant_8to18hpi = IFNmodel_foldchange_constant_0to18hpi-IFNmodel_foldchange_constant_0to8hpi;
            
            log2_IFNmodel_foldchange_saturatingplusX_18hpi = log2(IFNmodel_foldchange_saturating_8hpi+IFNmodel_foldchange_constant_8to18hpi);
            
        elseif model_ind == 2
            
            % linear (8-18 hpi)
            params_linear.r = lambda(3);
            IFNmodel_foldchange_linear_0to8hpi = get_IFN_linear_nb(this_moi, params_linear,alpha,hpi(1),k_dispersion);
            IFNmodel_foldchange_linear_0to18hpi = get_IFN_linear_nb(this_moi, params_linear,alpha,hpi(2),k_dispersion);
            IFNmodel_foldchange_linear_8to18hpi = IFNmodel_foldchange_linear_0to18hpi-IFNmodel_foldchange_linear_0to8hpi;
            
            log2_IFNmodel_foldchange_saturatingplusX_18hpi = log2(IFNmodel_foldchange_saturating_8hpi+IFNmodel_foldchange_linear_8to18hpi);
            
        else
            
            % saturating (8-18 hpi)
            params_saturating.r = lambda(3);
            params_saturating.K = lambda(4);
            IFNmodel_foldchange_saturating_0to8hpi = get_IFN_saturating_nb(this_moi,params_saturating,alpha,hpi(1),k_dispersion);
            IFNmodel_foldchange_saturating_0to18hpi = get_IFN_saturating_nb(this_moi,params_saturating,alpha,hpi(2),k_dispersion);
            IFNmodel_foldchange_saturating_8to18hpi = IFNmodel_foldchange_saturating_0to18hpi-IFNmodel_foldchange_saturating_0to8hpi;
            
            log2_IFNmodel_foldchange_saturatingplusX_18hpi = log2(IFNmodel_foldchange_saturating_8hpi+IFNmodel_foldchange_saturating_8to18hpi);
            
            
        end
        
        sum_sqerror = sum_sqerror + (this_18hpi_data - log2_IFNmodel_foldchange_saturatingplusX_18hpi)^2;

    end
    
end

f = sum_sqerror;