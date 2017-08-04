%% This matlab script validates analytical expression for pdf C_0 charterized using
%% Gamma approximation for different value of SNR_sp2

clear all;
close all;

analysis = 0;                                                   % Enable it for performing a complete new analysis
                                                                % otherwise a mat file is read 

if analysis
        M = 2e3;                                                % Number of realizations
        P_reg = 1;                                              % Regulate power at ST
        P_tran = 0;                                             % Powe transmitted by the PT
        alpha_s_true = 10^(-80/10);                             % Access Channel
        alpha_p2_true = 10.^(-[100]/10);                        % Secondary Interference Channel 
        N_p2 = [1000];                                          % Number of samples used for estimation
        noise_power = 10^(-100/10);                             % Noise power at SR 
        N_s = 10;                                               % Number of Pilot Symobols 
        m_p = [1 2 5];                                          % Nakagam-m fading paramter    
        m_s = [1 2 5];        
        snr_p2 = alpha_p2_true/noise_power;                     % interference to noise power at the SR
        snr_s = alpha_s_true/noise_power;                       % signal to noise power the the SR
        alpha_s_est = zeros(1,M);                               % Estimated power gain for the channel g_s        
        P_rcvd_sr_sim = zeros(1,M);                             % snr_p_prime = 1 + snr_p  
        CDF_Cs_sim = zeros(length(N_p2), length(snr_p2), M);    % Buffer to store the simulated CDF values 
        CDF_Cs_ana = zeros(length(N_p2),length(snr_p2), M);     % Buffer to store the analytic CDF values 
       
        Cs_pts = zeros(length(N_p2),length(snr_p2), M);         % Buffer to store the C1 points for the distribution function

        % Reinitialize the variaibles to ease the length of the function
        % that computes the integral
        g_s = alpha_s_true;
        np = noise_power;
        Preg = P_reg;
        

        for k = 1:length(m_p)
            disp(strcat('m = ',num2str(m_p(k))));
            for l = 1:length(alpha_p2_true)
		g_p2_est = ones(1, M) .* random('gam',m_p(k), 1/m_p(k), 1, M);
                g_s_est = ones(1, M) .* random('gam', m_s(k), 1/m_s(k), 1, M); 
              %% Simulation
                for i=1:M
                    alpha_s_est(i) = mean(normrnd(sqrt(g_s_est(i) * alpha_s_true), sqrt(noise_power),...
                        1, N_s).^2);
                    P_rcvd_sr_sim(i) = mean((sqrt(g_p2_est(i) * 10^((P_tran)/10) * alpha_p2_true(l))...
                        * ones(1, N_p2) + random('norm', 0, sqrt(noise_power), 1, N_p2)).^2);

                end
                

               %% CDF of Capacity C_1 = log2(1 + SNR_s)         

                 Cs = log2(1  + alpha_s_est./P_rcvd_sr_sim);
                 [CDF_Cs_temp, Cs_pts(k,l,:),~,~,eid] = cdfcalc(Cs);
                 CDF_Cs_sim(k,l,:) = CDF_Cs_temp(1:M);
                
                 for i = 1:length(Cs_pts(k,l,:))
                     func_den_Cs = @(t, v, u) log(2) .* exp(log(2.^t)  + (N_s*(1+u*g_s/np).^2./(2+u*4*g_s/np) - 1)...
                        .* log(2.^t - 1) - (N_s*(1+u*g_s/np).^2./(2+u*4*g_s/np) + N_p2*(1+v*snr_p2(l)).^2./(2+v*4*snr_p2(l))) .*...
                        log(1./(np/N_p2*(2+v*4*snr_p2(l))./(1+v*snr_p2(l))) + (2.^t - 1)./(np*Preg/N_s * (2+u*4*g_s/np)./(1+u*g_s/np)))...
                        + gammaln(N_p2*(1+v*snr_p2(l)).^2./(2+v*4*snr_p2(l)) + N_s*(1+u*g_s/np).^2./(2+u*4*g_s/np))...
                        - gammaln(N_p2*(1+v*snr_p2(l)).^2./(2+v*4*snr_p2(l))) - gammaln(N_s*(1+u*g_s/np).^2./(2+u*4*g_s/np))...
                        - N_p2*(1+v*snr_p2(l)).^2./(2+v*4*snr_p2(l)) .* log(np/N_p2*(2+v*4*snr_p2(l))./(1+v*snr_p2(l)))...
                        - N_s*(1+u*g_s/np).^2./(2+u*4*g_s/np) .* log(np*Preg/N_s*(2+u*4*g_s/np)./(1+u*g_s/np))) *...
                        1/gamma(m_p(k)) * m_p(k)^(m_p(k)) .* v.^(m_p(k) - 1) .* exp(-m_p(k) * v) *... 
                        1/gamma(m_s(k)) * m_s(k)^(m_s(k)) .* u.^(m_s(k) - 1) .* exp(-m_s(k) * u);               
                     CDF_Cs_ana(k,l,i) = integral3(func_den_Cs, 1, Cs_pts(k,l,i), 0, 100, 0, 100);
                 end
             end
        end        
        save('results_CDF_Cs_diff_m_fading.mat');
        quit;
end
load('results_CDF_Cs_diff_m_fading.mat')
        
% Plotting Curves
figure(1);
diff = 50;
sim = 1:diff:M;

k = 1; 
plot(reshape(Cs_pts(k, 1, :), 1, M), reshape(CDF_Cs_ana(k, 1, :),1, M), 'c-', 'Linewidth',1);
hold on,
plot(reshape(Cs_pts(k, 1, sim), 1, length(sim)), reshape(CDF_Cs_sim(k, 1, sim), 1, length(sim)),...
    'x', 'Linewidth', 1);

k = 2;
hold on,
plot(reshape(Cs_pts(k, 1, :), 1, M), reshape(CDF_Cs_ana(k, 1, :),1, M), 'c-', 'Linewidth',1);
hold on,
plot(reshape(Cs_pts(k, 1, sim), 1, length(sim)), reshape(CDF_Cs_sim(k, 1, sim), 1, length(sim)),...
    'x', 'Linewidth', 1);

k = 3;
hold on,
plot(reshape(Cs_pts(k, 1, :), 1, M), reshape(CDF_Cs_ana(k, 1, :),1, M), 'c-', 'Linewidth',1);
hold on,
plot(reshape(Cs_pts(k, 1, sim), 1, length(sim)), reshape(CDF_Cs_sim(k, 1, sim), 1, length(sim)),...
    'x', 'Linewidth', 1);

axis tight;
Fontsize = 8;
grid on;
xlabel('$\text{C}_s$ [bits/sec/Hz]','FontSize',Fontsize);
ylabel('CDF','FontSize',Fontsize);
hl = legend('Theoretical', 'Simulated');
set(hl, 'position',[0.62 0.12 0.28 0.14]);
set(gca,'FontSize',Fontsize);
laprint(1, '../figures/fig_CDF_Cs_diff_m_fading', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');        
