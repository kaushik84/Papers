%% This matlab script validates analytical expression for pdf C_1 charterized using
%% Gamma approximation for different value of Nakgama fading parameter m

clear all;
close all;

analysis = 1;                                                                           % Enable it for performing a complete new analysis
                                                                                        % otherwise a mat file is read 

if analysis
        % System Parameters
        M = 1e4;                                                                        % Number of realizations
        P_s = 10.^([-10]/10);                                                         % received snr for the channel g_s  
        test = 100;                                                                     % Estimation time  
        np = 10^(-100/10);                                                              % Noise power at SR 
        N_s = 10;                                                                       % Number of Pilot Symobols 
        est_er = np/N_s;                                                                % Variance of the Estimation error for h_s 
        m_s = [0.5 1 2 5 100];                                                          % Nakagam-m parameter
        m_p2 = [0.5 1 2 5 100];                                                         % Nakagam-m parameter
        g_s_true = 10^(-100/10);                                                        % True value of the path loss
        g_p2_true = 10^(-100/10);                                                       % True Path Loss for the channel  PT - ST 
        P_p = 10^(+00/10);                                                              % Transmit power at the PT, PR 
        rp = (P_p * g_p2_true + np);                                                     % True value of received energy
              
        % Buffers
        g_s_hat = zeros(1,M);                                                           % Estimated power gain for the channel g_s        
        rcvd_energy = zeros(1, M);                                                      % Received Energy at the secondary transmitter     
        CDF_C1_sim = zeros(length(m_s), M);                                             % Buffer to store the simulated CDF values 
        CDF_C1_ana = zeros(length(m_s), M);                                             % Buffer to store the analytic CDF values 
        C1_pts = zeros(length(m_s), M);                                                 % Buffer to store the C0 points for the distribution function
      
        for k = 1:length(m_s)
          %% Simulation
          %% g_hat_s is supposed to be non-central chi squared distributed with degree of freedom = 1 
            disp(strcat('m_s = ',num2str(m_s(k))));
            g_p2_sim = random('gam',m_p2(k), g_p2_true/m_p2(k), 1, M);                  % Power gain for the channel g_p2 (including path loss)
            g_s_sim = random('gam', m_s(k), g_s_true/m_s(k), 1, M);                     % Power gain for the channel g_s (including path loss)   
            for i=1:M
                pilot_s_hat(i) = mean(normrnd(sqrt(P_s * g_s_sim(i)), sqrt(np), 1, N_s).^2);
                rcvd_energy(i) = mean(random('norm',...
                    random('norm',0, sqrt(g_p2_sim(i) * P_p),1, test), sqrt(np), 1, test).^2);   
            end
          %% density of the SNR_s
            l_s = N_s * P_s * g_s_true / np;                                             % Lambda parameter of non-central chi square distribution
            mean_ana_s = (1/N_s) * (1*N_s + l_s);
            var_ana_s = (1/N_s)^2 * (2*N_s + 4*l_s);
         
            snr_s_hat = pilot_s_hat/np;               
            mean_sim_s = mean(snr_s_hat);
            var_sim_s = var(snr_s_hat);      

          %% Gammma Approximation
           if 0
                b_gamma_s = 1/N_s * (2*N_s + u*4*l_s)./(1*N_s + u*l_s);
                a_gamma_s = (1*N_s + u*l_s).^2./(2*N_s + u*4*l_s);
           end

           %% density for SNR_p_
            l_p2 = g_p2_true/np;  
            gp2t = g_p2_true; % Use a smaller notation

            v=1;
            mean_ana_p2 = test/2 * 2*(v*P_p*gp2t + np)/(test*np);
            var_ana_p2 = test/2 * (2*(v*P_p*gp2t + np)/(test*np))^2;
            mean_sim_p2 = mean(rcvd_energy/np);
            var_sim_p2 = var(rcvd_energy/np);    

            if 0
                b_gamma_p2 = 2*(v*P_p*gp2t + np)/(np*test);
                a_gamma_p2 = test/2;    
            end
          %% CDF of Capacity C_1 = log2(1 + SNR_s)         

            C_1 = log2(1  + np*snr_s_hat./rcvd_energy);
            [CDF_C1_temp, C1_pts(k,:),~,~,eid] = cdfcalc(C_1);
            CDF_C1_sim(k,:) = CDF_C1_temp(1:M);

            parfor i=1:M
                warning off;
                func_C1 = @(t, u, v) log(2) * exp(log(2.^t)  + ((1*N_s + u*l_s).^2./(2*N_s + u*4*l_s) - 1) .* log(2.^t - 1)...
                - ((1*N_s + u*l_s).^2./(2*N_s + u*4*l_s) + test/2) .* log(1./(2*(v*P_p*gp2t + np)/(test*np)) +...
                 (2.^t - 1)./(1/N_s * (2*N_s + u*4*l_s)./(1*N_s + u*l_s)))...
                + gammaln(test/2 + (1*N_s + u*l_s).^2./(2*N_s + u*4*l_s)) - gammaln(test/2) -...
                gammaln((1*N_s + u*l_s).^2./(2*N_s + u*4*l_s))...
                - test/2 * log(2 * (v*P_p*gp2t + np)/(test*np)) - (1*N_s + u*l_s).^2./(2*N_s + u*4*l_s) .*...
                log(1/N_s * (2*N_s + u*4*l_s)./(1*N_s + u*l_s)))...
                .* 1/gamma(m_s(k)) * m_s(k)^(m_s(k)) .* u.^(m_s(k) - 1) .* exp(-m_s(k) * u)...
                .* 1/gamma(m_p2(k)) * m_p2(k)^(m_p2(k)) .* v.^(m_p2(k) - 1) .* exp(-m_p2(k) * v);
                CDF_C1_ana(k,i) = integral3(func_C1, 0, C1_pts(k,i), 0, 100, 0, 100);
            end
            save('results_CDF_C1_diff_m_fading.mat');
        end        
        quit;
end
load('results_CDF_C1_diff_m_fading.mat');
        
% Plotting Curves
figure(1);
diff = 50;
sim = 1:diff:M;

k = 1; 
plot(C1_pts(k, :), CDF_C1_ana(k, :), 'c-', 'Linewidth',1);
hold on,
plot(C1_pts(k, sim), CDF_C1_sim(k, sim), 'x', 'Linewidth',1);

k = 2;
hold on,
plot(C1_pts(k, :), CDF_C1_ana(k, :), 'c-', 'Linewidth',1);
hold on,
plot(C1_pts(k, sim), CDF_C1_sim(k, sim), 'x', 'Linewidth',1);

k = 3;
hold on,
plot(C1_pts(k, :), CDF_C1_ana(k, :), 'c-', 'Linewidth',1);
hold on,
plot(C1_pts(k, sim), CDF_C1_sim(k, sim), 'x', 'Linewidth',1);

k = 4;
hold on,
plot(C1_pts(k, :), CDF_C1_ana(k, :), 'c-', 'Linewidth',1);
hold on,
plot(C1_pts(k, sim), CDF_C1_sim(k, sim), 'x', 'Linewidth',1);

k = 5;
hold on,
plot(C1_pts(k, :), CDF_C1_ana(k, :), 'c-', 'Linewidth',1);
hold on,
plot(C1_pts(k, sim), CDF_C1_sim(k, sim), 'x', 'Linewidth',1);

axis([0 2 0 1]);
Fontsize = 8;
grid on;
xlabel('$\text{C}_1$ [bits/sec/Hz]','FontSize',Fontsize);
ylabel('CDF','FontSize',Fontsize);
hl = legend('Theoretical', 'Simulated');
set(hl, 'position',[0.62 0.12 0.28 0.14]);
set(gca,'FontSize',Fontsize);
laprint(1, '../figures/fig_CDF_C1_diff_m_fading', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');        
