function [est_H3_positive_percent_vector, est_H1_positive_percent_vector, est_H3_positive_H1_positive_percent_vector] = GetEstimated_percents_coinfection_sim_fraction(moi_rH3N1_vector, params, params_fixed)

% 0 hpi
% superinfection at 6 hpi
% measured at 19 hpi

rH1N2_moi_scalar = params(1);
prob_gene_segment_expressed_PB2_PB1_PA = params(2);
prob_gene_segment_expressed_HA_rH3N1 = params(3);
prob_gene_segment_expressed_HA_rH1N2 = params(4);
fraction_x = params(5);
fraction_y = params(6);

max_cellular_moi = 10;  %ceil(2*max(moi_rH3N1_vector));
X = 0:max_cellular_moi; % different degrees of cellular input considered

prob_infected_cell_alive_by19hpi = 1 - params_fixed.prob_infected_cell_dead_19hpi;
prob_infected_cell_alive_by6hpi = 1 - params_fixed.prob_infected_cell_dead_6hpi;

cntr = 1; 
for rH3N1_moi_scalar = moi_rH3N1_vector
    
    rH3N1_low_S_MOI = rH3N1_moi_scalar*(fraction_x/fraction_y);
    poiss_rH3N1_low_S = poisspdf(X, rH3N1_low_S_MOI); poiss_rH3N1_low_S(end) = poiss_rH3N1_low_S(end) + (1-sum(poiss_rH3N1_low_S));
    rH3N1_high_S_MOI = rH3N1_moi_scalar*((1-fraction_x)/(1-fraction_y));
    poiss_rH3N1_high_S = poisspdf(X, rH3N1_high_S_MOI); poiss_rH3N1_high_S(end) = poiss_rH3N1_high_S(end) + (1-sum(poiss_rH3N1_high_S));
    
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
            
            prob_PB2 = 1 - binopdf(0,  X_rH3N1, prob_gene_segment_expressed_PB2_PB1_PA)*binopdf(0,  X_rH1N2, prob_gene_segment_expressed_PB2_PB1_PA); % probability that at least one PB2!
            prob_PB1 = 1 - binopdf(0,  X_rH3N1, prob_gene_segment_expressed_PB2_PB1_PA)*binopdf(0,  X_rH1N2, prob_gene_segment_expressed_PB2_PB1_PA); % probability that at least one PB1!
            prob_PA = 1 - binopdf(0,  X_rH3N1, prob_gene_segment_expressed_PB2_PB1_PA)*binopdf(0,  X_rH1N2, prob_gene_segment_expressed_PB2_PB1_PA); % probability that at least one PA!
            
            prob_good_polymerase = prob_PB2*prob_PB1*prob_PA;
            
            prob_H3 = 1 - binopdf(0,  X_rH3N1, prob_gene_segment_expressed_HA_rH3N1); 
            prob_H1 = 1 - binopdf(0,  X_rH1N2, prob_gene_segment_expressed_HA_rH1N2); 
            prob_H3_H1 = prob_H3*prob_H1; 
            prob_H3_positive = prob_good_polymerase*prob_H3;
            prob_H1_positive = prob_good_polymerase*prob_H1;
            prob_H3_positive_H1_positive = prob_good_polymerase*prob_H3_H1;
            
            est_H3_positive_percent = est_H3_positive_percent + prob_H3_positive*joint_poisson_matrix_lowS(X_rH3N1 + 1, X_rH1N2 + 1) + prob_H3_positive*joint_poisson_matrix_highS(X_rH3N1 + 1, X_rH1N2 + 1);
            est_H1_positive_percent = est_H1_positive_percent + prob_H1_positive*joint_poisson_matrix_lowS(X_rH3N1 + 1, X_rH1N2 + 1) + prob_H1_positive*joint_poisson_matrix_highS(X_rH3N1 + 1, X_rH1N2 + 1);
            est_H3_positive_H1_positive_percent = est_H3_positive_H1_positive_percent + prob_H3_positive_H1_positive*joint_poisson_matrix_lowS(X_rH3N1 + 1, X_rH1N2 + 1) + prob_H3_positive_H1_positive*joint_poisson_matrix_highS(X_rH3N1 + 1, X_rH1N2 + 1);

        end
    end
    percent_cells_uninfected = joint_poisson_matrix_lowS(1,1) + joint_poisson_matrix_highS(1,1);
    
    % now account for infected cell death!
    %est_H3_positive_percent_vector(cntr) = 100*est_H3_positive_percent;
    %est_H1_positive_percent_vector(cntr) = 100*est_H1_positive_percent;
    %est_H3_positive_H1_positive_percent_vector(cntr) = 100*est_H3_positive_H1_positive_percent;
    
    est_H3_positive_percent_vector(cntr) = 100*(prob_infected_cell_alive_by19hpi*est_H3_positive_percent/(prob_infected_cell_alive_by19hpi*(1-percent_cells_uninfected) + percent_cells_uninfected));
    est_H1_positive_percent_vector(cntr) = 100*(prob_infected_cell_alive_by19hpi*est_H1_positive_percent/(prob_infected_cell_alive_by19hpi*(1-percent_cells_uninfected) + percent_cells_uninfected));
    est_H3_positive_H1_positive_percent_vector(cntr) = 100*(prob_infected_cell_alive_by19hpi*est_H3_positive_H1_positive_percent/(prob_infected_cell_alive_by19hpi*(1-percent_cells_uninfected) + percent_cells_uninfected));
    cntr = cntr + 1;
end
