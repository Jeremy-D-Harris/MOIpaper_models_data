function [est_H3_positive_percent_vector, est_H1_positive_percent_vector, est_H3_positive_H1_positive_percent_vector] = GetEstimated_percents_superinfection_sim_fraction(moi_rH3N1_vector, params_coinfection, model_num, susc_reduction, params_fixed)

rH1N2_moi_scalar = params_coinfection(1);
prob_gene_segment_expressed_PB2_PB1_PA = params_coinfection(2);
prob_gene_segment_expressed_HA_rH3N1 = params_coinfection(3);
prob_gene_segment_expressed_HA_rH1N2 = params_coinfection(4);
fraction_x = params_coinfection(5);
fraction_y = params_coinfection(6);

%display('assuming that if x percent of cells alive by 6 hpi, only x percent of H1N2 virus gets in')

prob_infected_cell_alive_19hpi = 1 - params_fixed.prob_infected_cell_dead_19hpi;
prob_infected_cell_alive_6hpi = 1 - params_fixed.prob_infected_cell_dead_6hpi;

max_cellular_moi = 10; %ceil(2*max(moi_rH3N1_vector));
X = 0:max_cellular_moi; % different degrees of cellular input considered

cntr = 1; 
for rH3N1_moi_scalar = moi_rH3N1_vector
    
    rH3N1_low_S_MOI = rH3N1_moi_scalar*(fraction_x/fraction_y);
    poiss_rH3N1_low_S = poisspdf(X, rH3N1_low_S_MOI); poiss_rH3N1_low_S(end) = poiss_rH3N1_low_S(end) + (1-sum(poiss_rH3N1_low_S));
    rH3N1_high_S_MOI = rH3N1_moi_scalar*((1-fraction_x)/(1-fraction_y));
    poiss_rH3N1_high_S = poisspdf(X, rH3N1_high_S_MOI); poiss_rH3N1_high_S(end) = poiss_rH3N1_high_S(end) + (1-sum(poiss_rH3N1_high_S));
   
    % possibility 1: only a fraction p of H1N2 enters the remaining percentage p of alive cells.
    % possibility 2: only a fraction of H1N2 enters the remaining alive cells.
    %rH1N2_moi_scalar
    rH1N2_low_S_MOI = rH1N2_moi_scalar*(fraction_x/fraction_y);
    poiss_rH1N2_low_S = poisspdf(X, rH1N2_low_S_MOI); poiss_rH1N2_low_S(end) = poiss_rH1N2_low_S(end) + (1-sum(poiss_rH1N2_low_S));
    rH1N2_high_S_MOI = rH1N2_moi_scalar*((1-fraction_x)/(1-fraction_y));
    poiss_rH1N2_high_S = poisspdf(X, rH1N2_high_S_MOI); poiss_rH1N2_high_S(end) = poiss_rH1N2_high_S(end) + (1-sum(poiss_rH1N2_high_S));
   
    joint_poisson_matrix_lowS = poiss_rH3N1_low_S'*poiss_rH1N2_low_S*fraction_y;
    joint_poisson_matrix_highS = poiss_rH3N1_high_S'*poiss_rH1N2_high_S*(1-fraction_y);
    
    est_H3_positive_percent = 0; 
    est_H1_positive_percent = 0; 
    est_H3_positive_H1_positive_percent = 0;

    for X_rH3N1 = 0:max_cellular_moi
        for X_rH1N2 = 0:max_cellular_moi
            
            % loop around from 0 to X_rH1N2 to determine probabilities of 0
            % getting in getting susc_reduction, 1 getting in, etc. up to
            % X_rH1N2
            for n_H1N2_entered = 0:X_rH1N2
                if model_num == 0  % constant model:
                    if X_rH3N1 == 0
                        pmf_entered = binopdf(n_H1N2_entered,  X_rH1N2, 1);
                    else
                        pmf_entered = binopdf(n_H1N2_entered,  X_rH1N2, susc_reduction);
                    end
                elseif model_num == 1 % proportional model:
                    pmf_entered = binopdf(n_H1N2_entered,  X_rH1N2, susc_reduction^X_rH3N1);
                else
                    error('proportional_reduction needs to be 0 or 1');
                end
             
                % add pmf entered below
                prob_PB2 = 1 - binopdf(0,  X_rH3N1, prob_gene_segment_expressed_PB2_PB1_PA)*binopdf(0,  n_H1N2_entered, prob_gene_segment_expressed_PB2_PB1_PA); % probability that at least one PB2!
                prob_PB1 = 1 - binopdf(0,  X_rH3N1, prob_gene_segment_expressed_PB2_PB1_PA)*binopdf(0,  n_H1N2_entered, prob_gene_segment_expressed_PB2_PB1_PA); % probability that at least one PB1!
                prob_PA = 1 - binopdf(0,  X_rH3N1, prob_gene_segment_expressed_PB2_PB1_PA)*binopdf(0,  n_H1N2_entered, prob_gene_segment_expressed_PB2_PB1_PA); % probability that at least one PA!
             
                prob_good_polymerase = prob_PB2*prob_PB1*prob_PA;
                %prob_good_polymerase = prob_PB2*prob_PB1*prob_PA*prob_NP;
            
                prob_H3 = 1 - binopdf(0,  X_rH3N1, prob_gene_segment_expressed_HA_rH3N1); 
                prob_H1 = 1 - binopdf(0,  n_H1N2_entered, prob_gene_segment_expressed_HA_rH1N2); 
                prob_H3_H1 = prob_H3*prob_H1; 
                prob_H3_positive = prob_good_polymerase*prob_H3;
                prob_H1_positive = prob_good_polymerase*prob_H1;
                prob_H3_positive_H1_positive = prob_good_polymerase*prob_H3_H1;
            
                est_H3_positive_percent = est_H3_positive_percent + pmf_entered*prob_H3_positive*joint_poisson_matrix_lowS(X_rH3N1 + 1, X_rH1N2 + 1) + pmf_entered*prob_H3_positive*joint_poisson_matrix_highS(X_rH3N1 + 1, X_rH1N2 + 1);
                est_H1_positive_percent = est_H1_positive_percent + pmf_entered*prob_H1_positive*joint_poisson_matrix_lowS(X_rH3N1 + 1, X_rH1N2 + 1) + pmf_entered*prob_H1_positive*joint_poisson_matrix_highS(X_rH3N1 + 1, X_rH1N2 + 1);
                est_H3_positive_H1_positive_percent = est_H3_positive_H1_positive_percent + pmf_entered*prob_H3_positive_H1_positive*joint_poisson_matrix_lowS(X_rH3N1 + 1, X_rH1N2 + 1) + pmf_entered*prob_H3_positive_H1_positive*joint_poisson_matrix_highS(X_rH3N1 + 1, X_rH1N2 + 1);
            end
            
        end
    end
    
    percent_cells_uninfected = joint_poisson_matrix_lowS(1,1) + joint_poisson_matrix_highS(1,1);
    
    est_H3_positive_percent_vector(cntr) = 100*(prob_infected_cell_alive_19hpi*est_H3_positive_percent/(prob_infected_cell_alive_19hpi*(1-percent_cells_uninfected) + percent_cells_uninfected));
    est_H1_positive_percent_vector(cntr) = 100*(prob_infected_cell_alive_19hpi*est_H1_positive_percent/(prob_infected_cell_alive_19hpi*(1-percent_cells_uninfected) + percent_cells_uninfected));
    est_H3_positive_H1_positive_percent_vector(cntr) = 100*(prob_infected_cell_alive_19hpi*est_H3_positive_H1_positive_percent/(prob_infected_cell_alive_19hpi*(1-percent_cells_uninfected) + percent_cells_uninfected));
    
    cntr = cntr + 1;
end
