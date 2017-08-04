%% This function evaluates the capacities at the SR with and without intereference from the PR

function[C_0, C_1] = calc_capacities_fad(P_p, P_s, g_p2_true, g_s_true, noise_power, est_er, test, N_s, m_p2, m_s)

    %% Evaluate C_0    
    %% density of the SNR_s
    l_s = P_s * g_s_true / est_er;                                                     % Lambda parameter of non-central chi square distribution
   
    %% Gammma Approximation
    if 0
        mean_ana_s = (1/N_s) * (1*N_s + l_s);
        var_ana_s = (1/N_s)^2 * (2*N_s + 4*l_s);  
        u = 1;
        b_gamma_s = 1/N_s * (2*N_s + u*4*l_s)./(1*N_s + u*l_s);
        a_gamma_s = (1*N_s + u*l_s).^2./(2*N_s + u*4*l_s);
    end
    func_C0 = @(t, u) t .* 2.^(t) .* log(2) .* exp(((1*N_s + u*l_s).^2./(2*N_s + 4*u*l_s)) .* log(1./(1/N_s * (2*N_s + u*4*l_s)./(1*N_s + u*l_s)))...
        + ((1*N_s + u*l_s).^2./(2*N_s + 4*u*l_s) - 1) .* log(2.^t - 1) - gammaln((1*N_s + u*l_s).^2./(2*N_s + 4*u*l_s))...
        -(2.^t - 1)./(1/N_s * (2*N_s+ u*4*l_s)./(1*N_s + u*l_s))) .* ...
        1/gamma(m_s) * m_s^(m_s) .* u.^(m_s - 1) .* exp(-m_s * u);
    C_0 = integral2(func_C0, 0, 25, 0, 100);
    if isnan(C_0)
        error('C_0 evaluation resulted in an error: Check for numerical intergral');
    end
    
    
    %% Evaluate C_1    
    % Use a smaller notation
    gp2t = g_p2_true;                   
    np = noise_power;
        
    %% density for rcvd_power
    if 0
        v=1;
        mean_ana_p2 = test/2 * 2*(v*P_p*gp2t + np)/(test*np);
        var_ana_p2 = test/2 * (2*(v*P_p*gp2t + np)/(test*np))^2;
        mean_sim_p2 = mean(rcvd_energy/np);
        var_sim_p2 = var(rcvd_energy/np);  
    end
    func_C1 = @(t, u, v) t * log(2) .* exp(log(2.^t)  + ((1*N_s + u*l_s).^2./(2*N_s + u*4*l_s) - 1) .* log(2.^t - 1)...
        - ((1*N_s + u*l_s).^2./(2*N_s + u*4*l_s) + test/2) .* log(1./(2*(v*P_p*gp2t + np)/(test*np)) +...
         (2.^t - 1)./(1/N_s * (2*N_s + u*4*l_s)./(1*N_s + u*l_s)))...
        + gammaln(test/2 + (1*N_s + u*l_s).^2./(2*N_s + u*4*l_s)) - gammaln(test/2) -...
        gammaln((1*N_s + u*l_s).^2./(2*N_s + u*4*l_s))...
        - test/2 * log(2 * (v*P_p*gp2t + np)/(test*np)) - (1*N_s + u*l_s).^2./(2*N_s + u*4*l_s) .*...
        log(1/N_s * (2*N_s + u*4*l_s)./(1*N_s + u*l_s)))...
        .* 1/gamma(m_s) * m_s^(m_s) .* u.^(m_s - 1) .* exp(-m_s * u)...
        .* 1/gamma(m_p2) * m_p2^(m_p2) .* v.^(m_p2 - 1) .* exp(-m_p2 * v);
    C_1 = integral3(func_C1, 0, 25, 0, 100, 0, 100);
    if isnan(C_1)
        error('C_1 evaluation resulted in an error: Check for numerical intergral');
    end
end