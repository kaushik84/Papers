%% This function evaluates the capacities at the SR with and without intereference from the PR

function[C_1] = calc_capacities(snr_p2, alpha_s_true, P_reg, noise_power, N_p2, N_s)
    
    %% density of the P_reg * |h_s|^2 -- non-central chi2 distribution
    lambda_s =  alpha_s_true / noise_power;
    mean_ana_s = noise_power * P_reg * (1 + lambda_s);
    var_ana_s = (noise_power)^2 * P_reg^2/N_s * (2 + 4 * lambda_s);
    
    %% Gammma Approximation 
    b_gamma_s = var_ana_s/mean_ana_s;
    a_gamma_s = mean_ana_s/b_gamma_s; 
    

    %% density for (|\hp2| * P_tran + noise) -- non-central chi2 distribution
    mean_ana_p_ = noise_power * ( 1  + snr_p2);
    var_ana_p_ = noise_power^2/N_p2 * ( 2  + 4 * snr_p2);


    %% Gammma Approximation     
    b_gamma_p_ = var_ana_p_/mean_ana_p_;
    a_gamma_p_ = mean_ana_p_/b_gamma_p_;
    
    %% Evaluate C_1 = log2(1 + P_reg * |h_s|^2/ ((|\hp2| * P_tran + noise)) )
    func_exp_C1 = @(t) t .* log(2) .* exp(log(2.^t)  + (a_gamma_s - 1) .* log(2.^t - 1)...
            - (a_gamma_s + a_gamma_p_) * log(1/b_gamma_p_ + (2.^t - 1)/b_gamma_s)...
            + gammaln(a_gamma_p_ + a_gamma_s) - gammaln(a_gamma_p_) - gammaln(a_gamma_s)...
            - a_gamma_p_ * log(b_gamma_p_) - a_gamma_s * log(b_gamma_s));
    C_1 = integral(func_exp_C1, 0, 100);
    if isnan(C_1)
        error('C_1 evaluation resulted in an error: Check for numerical intergral');
    end

end