%% This matlab script validates analytical expression for pdf C_0 charterized using
%% Gamma approximation for different value of SNR_sp2

clear all;
close all;

analysis = 0;                                                   % Enable it for performing a complete new analysis
                                                                % otherwise a mat file is read 

if analysis
        M = 1e4;                                                % Number of realizations
        P_reg = 1;                                              % Regulate power at ST
        P_tran = 0;                                             % Powe transmitted by the PT
        alpha_s_true = 10^(-80/10);                             % Access Channel
        alpha_p2_true = 10.^(-[90 100 110]/10);                 % Secondary Interference Channel 
        N_p2 = 1000;                                             % Number of samples used for estimation
        noise_power = 10^(-100/10);                             % Noise power at SR 
        N_s = 10;                                               % Number of Pilot Symobols 
        snr_p2 = alpha_p2_true/noise_power;                     % interference to noise power at the SR
        snr_s = alpha_s_true/noise_power;                       % signal to noise power the the SR
        alpha_s_est = zeros(1,M);                               % Estimated power gain for the channel g_s        
        P_rcvd_sr_sim = zeros(1,M);                             % snr_p_prime = 1 + snr_p  
        CDF_Cs_sim = zeros(length(N_p2), length(snr_p2), M);    % Buffer to store the simulated CDF values 
        CDF_Cs_ana = zeros(length(N_p2),length(snr_p2), M);     % Buffer to store the analytic CDF values 
       
        Cs_pts = zeros(length(N_p2),length(snr_p2), M);         % Buffer to store the C1 points for the distribution function
        
        for k = 1:length(N_p2)
            for l = 1:length(alpha_p2_true)
              %% Simulation
                for i=1:M
                    alpha_s_est(i) = mean(normrnd(sqrt(alpha_s_true), sqrt(noise_power),...
                        1, N_s).^2);
                    P_rcvd_sr_sim(i) = mean((sqrt(10^((P_tran)/10) * alpha_p2_true(l))...
                        * ones(1, N_p2(k)) + random('norm', 0, sqrt(noise_power), 1, N_p2(k))).^2);

                end
                

               %% density of the P_reg * |h_s|^2 -- non-central chi2 distribution
                 lambda_s =  alpha_s_true / noise_power;
                 mean_ana_s = noise_power * P_reg * (1 + lambda_s);
                 var_ana_s = (noise_power)^2 * P_reg^2/N_s * (2 + 4 * lambda_s);

               %% Gammma Approximation 
                 b_gamma_s = var_ana_s/mean_ana_s;
                 a_gamma_s = mean_ana_s/b_gamma_s; 


               %% density for (|\hp2| * P_tran + noise) -- non-central chi2 distribution
                 mean_ana_p_ = noise_power * ( 1  + snr_p2(l));
                 var_ana_p_ = noise_power^2/N_p2(k) * ( 2  + 4 * snr_p2(l));


               %% Gammma Approximation     
                 b_gamma_p_ = var_ana_p_/mean_ana_p_;
                 a_gamma_p_ = mean_ana_p_/b_gamma_p_;   

               %% CDF of Capacity C_1 = log2(1 + SNR_s)         

                 Cs = log2(1  + alpha_s_est./P_rcvd_sr_sim);
                 [CDF_Cs_temp, Cs_pts(k,l,:),~,~,eid] = cdfcalc(Cs);
                 CDF_Cs_sim(k,l,:) = CDF_Cs_temp(1:M);

                 for i = 1:length(Cs_pts(k,l,:))
                     func_den_Cs = @(t) log(2) .* exp(log(2.^t)  + (a_gamma_s - 1) .* log(2.^t - 1)...
                        - (a_gamma_s + a_gamma_p_) * log(1/b_gamma_p_ + (2.^t- 1)/b_gamma_s)...
                        + gammaln(a_gamma_p_ + a_gamma_s) - gammaln(a_gamma_p_) - gammaln(a_gamma_s)...
                        - a_gamma_p_ * log(b_gamma_p_) - a_gamma_s * log(b_gamma_s));                  
                     CDF_Cs_ana(k,l,i) = integral(func_den_Cs, 1, Cs_pts(k,l,i));
                 end
             end
        end        
        save('results_CDF_Cs_diff_SNR_p2_AWGN.mat');
end
load('results_CDF_Cs_diff_SNR_p2_AWGN.mat')

diss = 0;  %% 1 = generate plot for the dissertation, 0 = generate plot for the Journal
if diss 
    % Plotting Curves
    figure(1);
    diff = 500;
    sim = 1:diff:M;

    k = 1; 
    plot(reshape(Cs_pts(1, k, :), 1, M), reshape(CDF_Cs_ana(1, k, :),1, M), 'c-', 'Linewidth',1);
    hold on,
    plot(reshape(Cs_pts(1, k, sim), 1, length(sim)), reshape(CDF_Cs_sim(1, k, sim), 1, length(sim)),...
        'x', 'Linewidth', 1);

    k = 2;
    hold on,
    plot(reshape(Cs_pts(1, k, :), 1, M), reshape(CDF_Cs_ana(1, k, :),1, M), 'c-', 'Linewidth',1);
    hold on,
    plot(reshape(Cs_pts(1, k, sim), 1, length(sim)), reshape(CDF_Cs_sim(1, k, sim), 1, length(sim)),...
        'o', 'Linewidth', 1);

    k = 3;
    hold on,
    plot(reshape(Cs_pts(1, k, :), 1, M), reshape(CDF_Cs_ana(1, k, :),1, M), 'c-', 'Linewidth',1);
    hold on,
    plot(reshape(Cs_pts(1, k, sim), 1, length(sim)), reshape(CDF_Cs_sim(1, k, sim), 1, length(sim)),...
        '+', 'Linewidth', 1);


    Fontsize = 8;
    axis tight;
    grid on;
    xlabel('$\eca$ [bits/sec/Hz]','FontSize',Fontsize);
    ylabel('CDF','FontSize',Fontsize);
    hl = legend('Theoretical', 'Simulated');
    set(hl, 'position',[0.26 0.12 0.28 0.14]);
    set(gca,'FontSize',Fontsize);
    laprint(1, '../figures/fig_CDF_Cs_diff_SNR_p2_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on'); 
else
    % Plotting Curves
    figure(1);
    diff = 500;
    sim = 1:diff:M;

    k = 1; 
    plot(reshape(Cs_pts(1, k, :), 1, M), reshape(CDF_Cs_ana(1, k, :),1, M), 'b-', 'Linewidth',1);
    hold on,
    plot(reshape(Cs_pts(1, k, sim), 1, length(sim)), reshape(CDF_Cs_sim(1, k, sim), 1, length(sim)),...
        'x', 'Linewidth', 1);

    k = 2;
    hold on,
    plot(reshape(Cs_pts(1, k, :), 1, M), reshape(CDF_Cs_ana(1, k, :),1, M), 'b-', 'Linewidth',1);
    hold on,
    plot(reshape(Cs_pts(1, k, sim), 1, length(sim)), reshape(CDF_Cs_sim(1, k, sim), 1, length(sim)),...
        'o', 'Linewidth', 1);

    k = 3;
    hold on,
    plot(reshape(Cs_pts(1, k, :), 1, M), reshape(CDF_Cs_ana(1, k, :),1, M), 'b-', 'Linewidth',1);
    hold on,
    plot(reshape(Cs_pts(1, k, sim), 1, length(sim)), reshape(CDF_Cs_sim(1, k, sim), 1, length(sim)),...
        '+', 'Linewidth', 1);


    Fontsize = 8;
    axis tight;
    
        
    lp1 = plot(0,0,'-b','Marker','x','visible','off');
    lp2 = plot(0,0,'-b','Marker','o','visible','off');
    lp3 = plot(0,0,'-b','Marker','+','visible','off');
    
    grid on;
    xlabel('$\eca$ [bits/sec/Hz]','FontSize',Fontsize);
    ylabel('CDF','FontSize',Fontsize);
    hl = legend([lp1, lp2, lp3], '$\frac{\pgpt\ptran}{\nps}=10$dB', '$\frac{\pgpt\ptran}{\nps}=0$dB',...
        '$\frac{\pgpt\ptran}{\nps}=-10$dB');    
    set(hl, 'position',[0.235 0.12 0.34 0.28]);
    set(gca,'FontSize',Fontsize);
    laprint(1, '../figures/fig_CDF_Cs_diff_SNR_p2_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on'); 
end