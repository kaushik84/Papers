%% This function evaluates the capacities at the SR with and without intereference from the PR

function[C_0, C_1] = calc_capacities(snr_p2, snr_s, noise_power, N_est)

    % Parameters concerning the ED and Matched filertering Models
    
    N_s = 10;                               % Number of Pilot Symobols 
    sigma_est_error = noise_power/N_s;      % Variance of the Estimation error for h_s 
    dof_s = 1;                              % Degree of freedom = 1
    K = N_est;                              % Number of samples used for energy detection for the signal at 

    g_s = snr_s * noise_power;              % Power gain for the channel g_s    
    tran_snr_ST = snr_s / g_s;              % Transmitted SNR = Transmitted power at PT/Noise at SR  = received SNR / g_s            


    %% Evaluate C_0
    
    %% density of the SNR_s
    lambda_s = dof_s *  g_s / sigma_est_error;
    mean_ana_s = sigma_est_error * tran_snr_ST * (dof_s + lambda_s);
    var_ana_s = (sigma_est_error)^2 * tran_snr_ST^2 * (2 * dof_s + 4 * lambda_s);
    
    %% Gammma Approximation 
    b_gamma_s = var_ana_s/mean_ana_s;
    a_gamma_s = mean_ana_s/b_gamma_s;
    
    func_exp_C0 = @(t) t .* 2.^(t) * log(2) .* exp(a_gamma_s * log(1/b_gamma_s)...
                    + (a_gamma_s - 1) * log(2.^t - 1) - gammaln(a_gamma_s)...
                    - (2.^t - 1)/b_gamma_s);
    C_0 = integral(func_exp_C0, 0, 25);
    if isnan(C_0)
        error('C_0 evaluation resulted in an error: Check for numerical intergral');
    end
    
    
    %% Evaluate C_1

    %% density for SNR_p_
    mean_ana_p_ = ( 1  + snr_p2);
    var_ana_p_ = 1/K * ( 2  + 4 * snr_p2);

    %% Gammma Approximation     
    b_gamma_p_ = var_ana_p_/mean_ana_p_;
    a_gamma_p_ = mean_ana_p_/b_gamma_p_;
    
    func_exp_C1 = @(t) t .* log(2) .* exp(log(2.^t)  + (a_gamma_s - 1) .* log(2.^t - 1)...
            - (a_gamma_s + a_gamma_p_) * log(1/b_gamma_p_ + (2.^t - 1)/b_gamma_s)...
            + gammaln(a_gamma_p_ + a_gamma_s) - gammaln(a_gamma_p_) - gammaln(a_gamma_s)...
            - a_gamma_p_ * log(b_gamma_p_) - a_gamma_s * log(b_gamma_s));
    C_1 = integral(func_exp_C1, 0, 25);
    if isnan(C_1)
        error('C_1 evaluation resulted in an error: Check for numerical intergral');
    end

end