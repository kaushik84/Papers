%% This matlab script validates analytical expression for pdf C_0 charterized using
%% Gamma approximation for different value of Nakgama fading parameter m

clear all;
close all;
analysis = 1;                                                                           % Enable it for performing a complete new analysis
                                                                                        % otherwise a mat file is read 

if analysis
        % System Parameters
        M = 1e4;                                                                        % Number of realizations
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
    
            
           %% Emperical CDF
            C_0 = log2(ones(1,M) +  P_s * g_s_true * random('gam', m_s(k), 1/m_s(k), 1, M)/np);
            [CDF_C0_temp, C0_pts(k,:),~,~,eid] = cdfcalc(C_0);axis([0 7 0 1]);

            CDF_C0_sim(k,:) = CDF_C0_temp(1:M);

           %% Analytic CDF
            CDF_C0_ana(k,:) = 1 - gammainc(m_s(k) * (2.^C0_pts(k,:) - 1) * np/(g_s_true * P_s),...
                m_s(k), 'upper');
        end
        save('results_CDF_C0_id_diff_m_fading.mat');
        %quit;
end
load('results_CDF_C0_id_diff_m_fading.mat');
        
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


load('results_CDF_C0_diff_m_fading.mat');
        
% Plotting Curves
figure(1);
diff = 500;
sim = 1:diff:M;

k = 1;
plot(C0_pts(k, :), CDF_C0_ana(k, :), 'g-', 'Linewidth',1);
hold on,
plot(C0_pts(k, sim), CDF_C0_sim(k, sim), 'x', 'Linewidth',1);


k = 2;
hold on,
plot(C0_pts(k, :), CDF_C0_ana(k, :), 'g-', 'Linewidth',1);
hold on,
plot(C0_pts(k, sim), CDF_C0_sim(k, sim), 'x', 'Linewidth',1);

k = 3;
hold on,
plot(C0_pts(k, :), CDF_C0_ana(k, :), 'g-', 'Linewidth',1);
hold on,
plot(C0_pts(k, sim), CDF_C0_sim(k, sim), 'x', 'Linewidth',1);

k = 4;
hold on,
plot(C0_pts(k, :), CDF_C0_ana(k, :), 'g-', 'Linewidth',1);
hold on,
plot(C0_pts(k, sim), CDF_C0_sim(k, sim), 'x', 'Linewidth',1);


k = 5;
hold on,
plot(C0_pts(k, :), CDF_C0_ana(k, :), 'g-', 'Linewidth',1);
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



%laprint(1, '../figures/fig_CDF_C0_id_diff_m_fading', 'options', 'factory', 'width', 8, 'scalefonts',...
%    'on', 'factor',0.5, 'keepfontprops', 'on');
