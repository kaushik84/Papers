%% This function evaluates the capacities at the SR with and without intereference from the PR

function[C_1] = calc_capacities_fad(snr_p2, g_s, Preg, np, N_p2, N_s, m_p, m_s)
   

    if 0 %T This is done to bring the right expressions in the de

        %% density of the P_reg * |h_s|^2
        %lambda_s = dof_s *  g_s / sigma_est_error;
        mean_ana_s = np * Preg * (1 + g_s/np);
        var_ana_s = np^2 * Preg^2/N_s * (2 + 4*g_s/np);

        %% Gammma Approximation 
        b_gamma_s = np*Preg/N_s*(2+u*4*g_s/np)./(1+u* g_s/np);
        %var_ana_s/mean_ana_s;
        a_gamma_s = N_s*(1+u*g_s/np).^2./(2+u*4*g_s/np); 
        %mean_ana_s/b_gamma_s;    


        %% density for SNR_p_
        mean_ana_p_ = np*(1+snr_p2);
        var_ana_p_ = np^2/N_p2*(2+4*snr_p2);

        %% Gammma Approximation     
        b_gamma_p_ = np/N_p2*(2+v*4*snr_p2)./(1+v*snr_p2);
        %var_ana_p_/mean_ana_p_;
        a_gamma_p_ = N_p2*(1+v*snr_p2).^2./(2+v*4*snr_p2);
        %mean_ana_p_/b_gamma_p_;

    end
    %u = 1,v = 1;
    
    func_exp_C1 = @(t, v, u) t .* log(2) .* exp(log(2.^t)  + (N_s*(1+u*g_s/np).^2./(2+u*4*g_s/np) - 1)...
        .* log(2.^t - 1) - (N_s*(1+u*g_s/np).^2./(2+u*4*g_s/np) + N_p2*(1+v*snr_p2).^2./(2+v*4*snr_p2)) .*...
        log(1./(np/N_p2*(2+v*4*snr_p2)./(1+v*snr_p2)) + (2.^t - 1)./(np*Preg/N_s * (2+u*4*g_s/np)./(1+u*g_s/np)))...
        + gammaln(N_p2*(1+v*snr_p2).^2./(2+v*4*snr_p2) + N_s*(1+u*g_s/np).^2./(2+u*4*g_s/np))...
        - gammaln(N_p2*(1+v*snr_p2).^2./(2+v*4*snr_p2)) - gammaln(N_s*(1+u*g_s/np).^2./(2+u*4*g_s/np))...
        - N_p2*(1+v*snr_p2).^2./(2+v*4*snr_p2) .* log(np/N_p2*(2+v*4*snr_p2)./(1+v*snr_p2))...
        - N_s*(1+u*g_s/np).^2./(2+u*4*g_s/np) .* log(np*Preg/N_s*(2+u*4*g_s/np)./(1+u*g_s/np))) *...
        1/gamma(m_p) * m_p^(m_p) .* v.^(m_p - 1) .* exp(-m_p * v) *... 
        1/gamma(m_s) * m_s^(m_s) .* u.^(m_s - 1) .* exp(-m_s * u);
    %C_1 = integral(func_exp_C1, 0, 100);
    %C_1 = integral2(func_exp_C1, 0, 25, 0 , 100);

    C_1 = integral3(func_exp_C1, 0, 100, 0 , 100, 0, 100);
    if isnan(C_1)
        error('C_1 evaluation resulted in an error: Check for numerical intergral');
    end

end