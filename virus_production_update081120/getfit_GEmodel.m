function sumsq_error = getfit_GEmodel(log_nu, moi_vector, data_allhpi, all_hpi, alpha, k_dispersion, Ncells, cell_line, model_num, which_time_model)

nu = exp(log_nu);

sumsq_error = 0;
for n = 1:length(all_hpi)
    
    this_hpi = all_hpi(n);
    these_moi = [];
    these_data = [];
    
    % data at this hpi
    these_data_thishpi = data_allhpi(n,:);
    
    [vals, nonzero_indices] = find(1-isnan(these_data_thishpi));
    these_moi = moi_vector(nonzero_indices);
    these_data = these_data_thishpi(nonzero_indices);
    
    these_GE_outputs = get_GEmodel(nu,these_moi,Ncells,alpha,k_dispersion,this_hpi,cell_line,model_num,which_time_model);
    
    sumsq_error = sumsq_error + sum((log(these_data) - log(these_GE_outputs)).^2); % log-scale
    
end


