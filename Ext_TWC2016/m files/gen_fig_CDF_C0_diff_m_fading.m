%% This matlab script validates analytical expression for pdf C_0 charterized using
%% Gamma approximation for different value of Nakgama fading parameter m

clear all;
close all;
analysis = 0;                                                                           % Enable it for performing a complete new analysis
                                                                                        % otherwise a mat file is read 

if analysis
        % System Parameters
        M = 1e3;                                                                        % Number of realizations
        P_s = 10.^([-10]/10);                                                           % snr of the received signal for the channel g_s  
        np = 10^(-100/10);                                                              % transmitted power
        N_s = 10;                                                                       % Number of Pilot Symobols 
        est_er = np/N_s;                                                                % Variance of the Estimation error for h_s 
        m_s = [0.5 1 2 5 100];                                                          % Nakagam-m parameter
        g_s_true = 10.^([-80]/10);                                                      % True value of the path loss

        % Buffers
        g_s_hat = zeros(1,M);                                                           % Estimated power gain for the channel g_s        
        CDF_C0_sim = zeros(length(m_s), M);                                             % Buffer to store the simulated CDF values 
        CDF_C0_ana = zeros(length(m_s), M);                                             % Buffer to store the analytic CDF values 
        C0_pts = zeros(length(m_s), M);                                                 % Buffer to store the C0 points for the distribution function
        

        for k = 1:length(m_s) 
           %% Simulation
           %% g_hat_s is supposed to be non-central chi squared distributed with degree of freedom = 1 
            disp(strcat('m_s = ',num2str(m_s(k))));
            g_s_sim = g_s_true * random('gam', m_s(k), 1/m_s(k), 1, M);                  % Power gain for the channel g_s (including path loss)   
            snr_t = P_s / g_s_true;                                                      % Transmitted SNR = Transmitted power at PT/Noise at SR  = received SNR / g_s            
            for i=1:M
                pilot_s_hat(i) = mean(normrnd(sqrt(P_s * g_s_sim(i)), sqrt(np), 1, N_s).^2);
            end

           %% density of the SNR_s
            l_s = N_s * P_s * g_s_true / np;                                                     % Lambda parameter of non-central chi square distribution
            mean_ana_s = (1/N_s) * (1 * N_s + l_s);
            var_ana_s = (1/N_s)^2 * (2 * N_s + 4 * l_s);
            
            snr_s_hat = pilot_s_hat/(np);               
            mean_sim_s = mean(snr_s_hat);
            var_sim_s = var(snr_s_hat);        

          %% Gammma Approximation
           if 0
                b_gamma_s = 1/N_s * (2*N_s + u*4*l_s)./(1*N_s + u*l_s);
                a_gamma_s = (1*N_s + u*l_s).^2./(2*N_s + u*4*l_s);
           end

           %% CDF of Capacity C_0 = log2(1 + SNR_s)         
            
           %% Emperical CDF
            C_0 = log2(ones(1,M) +  snr_s_hat);
            [CDF_C0_temp, C0_pts(k,:),~,~,eid] = cdfcalc(C_0);axis([0 7 0 1]);

            CDF_C0_sim(k,:) = CDF_C0_temp(1:M);

           %% Analytic CDF
            parfor j=1:M
                %func_C0 = @(t) 2.^(t) * log(2) .* 1/b_gamma_s^(a_gamma_s)...
                %    .* (2.^t - 1).^(a_gamma_s - 1)/gamma(a_gamma_s)...
                %    .* exp(- (2.^t - 1)/b_gamma_s);
                warning off;
                func_C0 = @(t, u) 2.^(t) * log(2) .* exp(((1*N_s + u*l_s).^2./(2*N_s + 4*u*l_s)) .* log(1./(1/N_s * (2*N_s + u*4*l_s)./(1*N_s + u*l_s)))...
                    + ((1*N_s + u*l_s).^2./(2*N_s + 4*u*l_s) - 1) .* log(2.^t - 1) - gammaln((1*N_s + u*l_s).^2./(2*N_s + 4*u*l_s))...
                    -(2.^t - 1)./(1/N_s * (2*N_s+ u*4*l_s)./(1*N_s + u*l_s))) .* ...
                    1/gamma(m_s(k)) * m_s(k)^(m_s(k)) .* u.^(m_s(k) - 1) .* exp(-m_s(k) * u);
                %CDF_C0_ana(k,j) = integral(func_C0, 0, C0_pts(k,j));

                CDF_C0_ana(k,j) = integral2(func_C0, 0, C0_pts(k,j), 0, 100);
            end
        end
        save('results_CDF_C0_diff_m_fading.mat');
        quit;
end
load('results_CDF_C0_diff_m_fading.mat');
        
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

k = 4;
hold on,
plot(C0_pts(k, :), CDF_C0_ana(k, :), 'c-', 'Linewidth',1);
hold on,
plot(C0_pts(k, sim), CDF_C0_sim(k, sim), 'x', 'Linewidth',1);


k = 5;
hold on,
plot(C0_pts(k, :), CDF_C0_ana(k, :), 'c-', 'Linewidth',1);
hold on,
plot(C0_pts(k, sim), CDF_C0_sim(k, sim), 'x', 'Linewidth',1);

Fontsize = 8;
axis([0 7 0 1]);
grid on;
xlabel('$\text{C}_0$ [bits/sec/Hz]','FontSize',Fontsize);
ylabel('CDF','FontSize',Fontsize);
hl = legend('Theoretical', 'Simulated');
set(hl, 'position',[0.62 0.12 0.28 0.14]);
set(gca,'FontSize',Fontsize);
laprint(1, '../figures/fig_CDF_C0_diff_m_fading', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');
