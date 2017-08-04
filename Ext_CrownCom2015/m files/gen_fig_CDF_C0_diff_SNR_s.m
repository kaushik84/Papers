%% This matlab script validates analytical expression for pdf C_0 charterized using
%% Gamma approximation for different value of SNR_s

clear all;
close all;

analysis = 0;                                           % Enable it for performing a complete new analysis
                                                        % otherwise a mat file is read 

if analysis
        M = 1e4;                                % Number of realizations
        snr_s = 10.^([-10 0 10]/10);             % received snr for the channel g_s  
        noise_power = 10^(-100/10);             % Noise power at SR 
        E_s = 1;                                % Pilot Energy 
        N_s = 10;                               % Number of Pilot Symobols 
        sigma_est_error = noise_power/N_s;      % Variance of the Estimation error for h_s 
        dof_s = 1;                              % Degree of freedom = 1
        g_s_hat = zeros(1,M);                   % Estimated power gain for the channel g_s        
        CDF_C0_sim = zeros(length(snr_s), M);   % Buffer to store the simulated CDF values 
        CDF_C0_ana = zeros(length(snr_s), M);   % Buffer to store the analytic CDF values 
        
        C0_pts = zeros(length(snr_s),M);        % Buffer to store the C0 points for the distribution function
        
        for k = 1:length(snr_s)
           %% Simulation
           %% g_hat_s is supposed to be non-central chi squared distributed with degree of freedom = 1 
            g_s = snr_s(k) * noise_power;        % Power gain for the channel g_s (including path loss)   
            tran_snr_ST = snr_s(k) / g_s;        % Transmitted SNR = Transmitted power at PT/Noise at SR  = received SNR / g_s            
            for i=1:M
                g_s_hat(i) = sum(normrnd(sqrt(g_s), sqrt(sigma_est_error), 1, dof_s).^2);
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

           %% CDF of Capacity C_0 = log2(1 + SNR_s)         
            
           %% Emperical CDF
            C_0 = log2(ones(1,M) +  snr_s_hat);
            [CDF_C0_temp, C0_pts(k,:),~,~,eid] = cdfcalc(C_0);
            CDF_C0_sim(k,:) = CDF_C0_temp(1:M);

           %% Analytic CDF
            for j=1:M
                func_C0 = @(t) 2.^(t) * log(2) .* 1/b_gamma_s^(a_gamma_s)...
                    .* (2.^t - 1).^(a_gamma_s - 1)/gamma(a_gamma_s)...
                    .* exp(- (2.^t - 1)/b_gamma_s);
                CDF_C0_ana(k,j) = integral(func_C0, 0, C0_pts(k,j));
            end
        end
        save('results_CDF_C0_diff_SNR_s.mat');
end
load('results_CDF_C0_diff_SNR_s.mat');
        
% Plotting Curves
figure(1);
diff = 500;
sim = 1:diff:M;

k = 1;
plot(C0_pts(k, :), CDF_C0_ana(k, :), 'c-', 'Linewidth',1);
hold on,
plot(C0_pts(k, sim), CDF_C0_sim(k, sim), 'x', 'Linewidth',1);


k = 2;
hold on,
plot(C0_pts(k, :), CDF_C0_ana(k, :), 'c-', 'Linewidth',1);
hold on,
plot(C0_pts(k, sim), CDF_C0_sim(k, sim), 'x', 'Linewidth',1);

k = 3;
hold on,
plot(C0_pts(k, :), CDF_C0_ana(k, :), 'c-', 'Linewidth',1);
hold on,
plot(C0_pts(k, sim), CDF_C0_sim(k, sim), 'x', 'Linewidth',1);

Fontsize = 8;
grid on;
xlabel('$\text{C}_0$ [bits/sec/Hz]','FontSize',Fontsize);
ylabel('CDF','FontSize',Fontsize);
hl = legend('Theoretical', 'Simulated');
set(hl, 'position',[0.26 0.12 0.28 0.14]);
set(gca,'FontSize',Fontsize);
laprint(1, '../figures/fig_CDF_C0_diff_SNR_s', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');