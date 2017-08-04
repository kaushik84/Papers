%% This matlab script validates analytical expression for pdf C_0 charterized using
%% Gamma approximation for different value of SNR_s

clear all;
close all;

analysis = 0;                                           % Enable it for performing a complete new analysis
                                                        % otherwise a mat file is read 

if analysis
        M = 1e4;                                                % Number of realizations
        snr_s = 10.^([-10 0 10]/10);                            % received snr for the channel g_s  
        snr_p2 = 10.^([10]/10);                                 % snr for the channel g_p2  
        K = 1000;                                               % Number of samples used for energy detection 
        noise_power = 10^(-100/10);                             % Noise power at SR 
        E_s = 1;                                                % Pilot Energy 
        N_s = 10;                                               % Number of Pilot Symobols 
        sigma_est_error = noise_power/N_s;                      % Variance of the Estimation error for h_s 
        dof_s = 1;                                              % Degree of freedom = 1
        g_s_hat = zeros(1,M);                                   % Estimated power gain for the channel g_s        
        snr_p_ = zeros(1,M);                                    % snr_p_prime = 1 + snr_p  
        CDF_C1_sim = zeros(length(snr_s), length(snr_p2), M);   % Buffer to store the simulated CDF values 
        CDF_C1_ana = zeros(length(snr_s),length(snr_p2), M);    % Buffer to store the analytic CDF values 
       
        C1_pts = zeros(length(snr_s),length(snr_p2), M);         % Buffer to store the C1 points for the distribution function
        
        for k = 1:length(snr_s)
            for l = 1:length(snr_p2)
                %% Simulation
                %% g_hat_s is supposed to be non-central chi squared distributed with degree of freedom = 1 
                g_s = snr_s(k) * noise_power;        % Power gain for the channel g_s    
                tran_snr_ST = snr_s(k) / g_s;        % Transmitted SNR = Transmitted power at PT/Noise at SR  = received SNR / g_s            
                for i=1:M
                    g_s_hat(i) = sum(normrnd(sqrt(g_s), sqrt(sigma_est_error), 1, dof_s).^2);
                    snr_p_(i) = mean((sqrt(snr_p2(l)) * ones(1, K)...
                    + random('norm', 0, 1, 1, K)).^2);
                end

                %% density of the SNR_s
                lambda_s = dof_s *  g_s / sigma_est_error;
                snr_s_hat = tran_snr_ST * (g_s_hat);                
                mean_ana_s = sigma_est_error * tran_snr_ST * (dof_s + lambda_s);
                var_ana_s = (sigma_est_error)^2 * tran_snr_ST^2 * (2 * dof_s + 4 * lambda_s);
                mean_sim_s = mean(snr_s_hat);
                var_sim_s = var(snr_s_hat);        

                %% Gammma Approximation 
                b_gamma_s = var_ana_s/mean_ana_s;
                a_gamma_s = mean_ana_s/b_gamma_s;


                %% density for SNR_p_
                mean_ana_p_ = ( 1  + snr_p2(l));
                var_ana_p_ = 1/K * ( 2  + 4 * snr_p2(l));
                mean_sim_p_ = mean(snr_p_);
                var_sim_p_ = var(snr_p_);    

                b_gamma_p_ = var_ana_p_/mean_sim_p_;
                a_gamma_p_ = mean_sim_p_/b_gamma_p_;    

                %% CDF of Capacity C_1 = log2(1 + SNR_s)         

                C_1 = log2(1  + snr_s_hat./snr_p_);
                [CDF_C1_temp, C1_pts(k,l,:),~,~,eid] = cdfcalc(C_1);
                CDF_C1_sim(k,l,:) = CDF_C1_temp(1:M);

                for j=1:M
                    func_C1 = @(t) log(2) * exp(log(2.^t)  + (a_gamma_s - 1) .* log(2.^t - 1)...
                    - (a_gamma_s + a_gamma_p_) * log(1/b_gamma_p_ + (2.^t - 1)/b_gamma_s)...
                    + gammaln(a_gamma_p_ + a_gamma_s) - gammaln(a_gamma_p_) - gammaln(a_gamma_s)...
                    - a_gamma_p_ * log(b_gamma_p_) - a_gamma_s * log(b_gamma_s));
                    CDF_C1_ana(k,l,j) = integral(func_C1, 0, C1_pts(k,l,j));
                end
            end
        end        
        save('results_CDF_C1_diff_SNR_s_SNR_p2_10.mat');
end
load('results_CDF_C1_diff_SNR_s_SNR_p2_10.mat')
        
% Plotting Curves
figure(1);
diff = 10;
sim = 1:diff:M;

k = 1; 
plot(reshape(C1_pts(k, 1, :), 1, M), reshape(CDF_C1_ana(k, 1, :),1, M), 'c-', 'Linewidth',5);
hold on,
plot(reshape(C1_pts(k, 1, sim), 1, length(sim)), reshape(CDF_C1_sim(k, 1, sim), 1, length(sim)),...
    '-', 'Linewidth', 2);

k = 2;
hold on,
plot(reshape(C1_pts(k, 1, :), 1, M), reshape(CDF_C1_ana(k, 1, :),1, M), 'c-', 'Linewidth',5);
hold on,
plot(reshape(C1_pts(k, 1, sim), 1, length(sim)), reshape(CDF_C1_sim(k, 1, sim), 1, length(sim)),...
    '-', 'Linewidth', 2);

k = 3;
hold on,
plot(reshape(C1_pts(k, 1, :), 1, M), reshape(CDF_C1_ana(k, 1, :),1, M), 'c-', 'Linewidth',5);
hold on,
plot(reshape(C1_pts(k, 1, sim), 1, length(sim)), reshape(CDF_C1_sim(k, 1, sim), 1, length(sim)),...
    '-', 'Linewidth', 2);


Fontsize = 9;
grid on;
xlabel('$\text{C}_1$ [bits/sec/Hz]','FontSize',Fontsize);
ylabel('CDF','FontSize',Fontsize);
hl = legend('Theoretical', 'Simulated');
set(hl, 'position',[0.56 0.12 0.34 0.14]);
set(gca,'FontSize',Fontsize);
laprint(1, '../figures/fig_CDF_C1_diff_SNR_s_SNR_p2_10', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');        