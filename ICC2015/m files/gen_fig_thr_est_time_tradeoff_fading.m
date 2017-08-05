%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: This m-file considers the Throughput-Estimation time 
%        tradeoff for the Rayliegh fading channel. The fading is fast        
%
%        Clearly the regulated power is distributed according non central
%        chi-2 distribution. However to simplify the analysis we
%        approximate the non central distribution against the Gamma
%        distribution.
%        
%        The simulation is performed to show the following:
%        1) Analyse the throughput vs estimation curves
%        2) Confirm the theoretical analysis vs simulation results
%
%        Created on: 01.09.14
%        Last modified: 01.09.14
%        Revision History: 01.09.14 --> File generated   
%       
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

sim = 1;                                  % Enable( = 1) to perform simulation, theoretical analysis
th = 1;                                   % disable ( = 0) to plot curves, data is read from a file
                                          % Theoretical anaylsis is also included of simulation as
                                          % numerical integration is involved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noise_power = -100;                       % noise power -100 dBm
I_T = -110;                               % Interference temperature -80 dBm
f_s = 1e6;                                % 1 MHz one band
SNR_rcvd = -00;                           % SNR received at SR 0 dB
K = 0.1 * f_s;                            % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha_true = 90;                          % True Path loss between ST and PR  to achieve P_p = IT  
P_rcvd = SNR_rcvd + noise_power;          % Power received at PR
mu_P = 0.025;                             % Accuracy at PR
mu_R = 0.025;                             % Accuracy at SR:

m_p = 1;                                  % Number of realizations of channel for the link ST-PR
m_s = 1;                                  % Number of realizations of channel for the link ST-SR 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 1e3;                                  % Number of realizations 
N_temp = 4;
N_sim = ceil([10:4:50] * (1e-3 * f_s));   % N = Total of samples used for estimation = tau * f_s 
P_reg_sim = zeros(length(N_sim), M);      % Regulated transmitted power at ST
P_p_sim = zeros(length(N_sim), M);        % Power received at the PR
P_s_sim = zeros(length(N_sim), M);        % Power received at the SR
R_sim = zeros(length(N_sim), M);          % Throughput at the SR
Normalized_dist = zeros(1, M);            % Normalized distribution retrieved fro the received power

Exp_R_sim = zeros(1,length(N_sim));       % Expected throughput  
Exp_P_reg_sim= zeros(1,length(N_sim));    % Expected power regulated by ST
Exp_P_p_sim= zeros(1,length(N_sim));      % Expected power received at PR
PC_P_p_sim = zeros(1,length(N_sim));      % Confidence probability power at PR  
PC_R_sim = zeros(1,length(N_sim));        % Condidence probability rate at SR

if sim
    for i=1:length(N_sim)        
        disp(strcat('N = ',num2str(i)));        
        g_p = ones(1, M) .* random('gam',m_p, 1, 1, M);
        g_s = random('gam',m_s, 1, 1, M); 
        for j=1:M
            % Estimate the received power P_reg with additive noise amplitude (zero mean and noise variance) 
            
            samples = sqrt(g_p(j) * 10^(-alpha_true/10)) * ones(1, N_sim(i)) + random('norm', 0, sqrt(10^(noise_power/10)), 1, N_sim(i));
            
            % Checkpoint 
            % Here the samples must be Gaussian distributed with mean = 10^(P_reg/10), 
            % examine by ploting the histogram, hist(samples, 100)
            Normalized_dist(j) =  10^(-alpha_true/10)/(10^(-alpha_true/10) + 10^(noise_power/10)) * mean(samples.^2);  % 1/(10^(P_rcvd/10) + 10^(noise_power/10)) *        
        end       
            
            % Add noise margin in the path Loss 
            channel_sim(i,:) = Normalized_dist(:)'; % 10^(-alpha_true/10) .* 
                   
            
            % Perform power control over the received samples, 
            % P_reg_sim should be Inverse gamma distributed, 
            % hist(P_reg_sim(i,:), 50)
            P_reg_sim(i,:) =  10^(I_T/10) ./ (channel_sim(i,:));           
            
            % Determine P_p, P_s and R from P_reg after imposing path loss
            % hist(P_p_sim(i,:), 50)
            P_p_sim(i,:) = g_p * 10^(-alpha_true/10) .* P_reg_sim(i,:);
            % hist(P_s_sim(i,:), 50)
            P_s_sim(i,:) = g_s * 10^(-alpha_true/10) .* P_reg_sim(i,:);
            % hist(R_sim(i,:), 50)
            R_sim(i,:) = (K- N_sim(i))/K * log2(1 + P_s_sim(i,:) / 10^(noise_power/10));  
            
           %% Determining the performance meterics --> Exp_R, CP_power, CP_rate
            
            % Expected values
            Exp_P_reg_sim(i) = mean(P_reg_sim(i,:));
            Exp_R_sim(i) = mean(R_sim(i,:));
            Exp_P_p_sim(i) = mean(P_p_sim(i,:));
            
            % PC power
            upper_limit = Exp_P_p_sim(i) * (1 + mu_P);
            lower_limit = Exp_P_p_sim(i) * (1 - mu_P);
            PC_P_p_sim(i) = 1/M * length(find((P_p_sim(i,:) < upper_limit) & (P_p_sim(i,:) > lower_limit)));
            
            % PC rate
            upper_limit = Exp_R_sim(i) * (1 + mu_R);
            lower_limit = Exp_R_sim(i) * (1 - mu_R);
            PC_R_sim(i) = 1/M * length(find((R_sim(i,:) < upper_limit) & (R_sim(i,:) > lower_limit)));
            
            disp(strcat('Exp_R_sim(i) = ',num2str(Exp_R_sim(i)))); 
            disp(strcat('PC_power_sim(i) = ',num2str(PC_P_p_sim(i)))); 
            disp(strcat('PC_rate_sim(i) = ',num2str(PC_R_sim(i)))); 

    end
    save('results_thr_est_time_tradeoff_fading_sim.mat');
    if ~th      % If theoretical analysis is already simulated the computation ends here.
        quit;
    end
end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Theoretical analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th
    % Parameters of Received Power
    N_th = ceil([N_temp:0.1:N_temp] * (1e-3 * f_s));                 % N = Total of samples used for estimation = tau * f_s 
    Exp_R_th = zeros(1, length(N_th));

    % Gamma approximation to non-central distribution for the normalized with fading parameter
    a_norm_dis = N_th * (1 + m_p * 10^(SNR_rcvd/10))^2/(2 + 4 * m_p * 10^(SNR_rcvd/10) + N_th * m_p * (10^(SNR_rcvd/10))^2) ;
    b_norm_dis = 1/N_th * (2 + 4 * m_p * 10^(SNR_rcvd/10) + N_th * m_p * (10^(SNR_rcvd/10)^2) ) * (10^(noise_power/10))/(1 + m_p * 10^(SNR_rcvd/10)) / (10^(P_rcvd/10) + 10^(noise_power/10));
    b_channel = b_norm_dis*10^(-alpha_true/10);              % Gamma distributed
    %b_P_reg = b_alpha_true/10^(I_T/10);                     % Inverse Gamma distributed
    %b_P_p = b_P_reg/10^(-alpha_true/10);                    % Inverse Gamma distributed
    %b_P_s = b_P_reg/10^(-alpha_true/10);                    % Distribution is not known

    bins = 100;
    [f x_pts] = hist(P_p_sim, bins);
    
      
    
    %% Determining the performance meterics --> Exp_R, CP_power, CP_rate
    % Expected values
%     Exp_alpha_true = a .* b_alpha_true;
%     Exp_P_reg_th = 1./(b_P_reg .* (a - 1));

%     % Numerical integration
%     s = 10^(noise_power/10) * b_P_s;
%     for i=1:length(N_th)
%         %format long
%         %func = @(t) (K- N_th(i))/K * exp( - gammaln(a(i)) - t + log(log2(1 + 1./(s(i) * t))) + (a(i) - 1) * log(t)); 
%         %func = @(t) 1/log(2) * (K- N(i))/K * gamma(a(i)) .* exp(-t) .* log(1 + 1./(s(i) * t)) .* t.^(a(i) - 1); 1/log(2) * (K- N_th(i))/K  *
%         %Exp_R_th(i) = integral(func, 0, Inf);
%         error = 0.01;
%         t = ceil(0.01*a(i)):error:ceil(10*a(i));
%         y = (K- N_th(i))/K * exp( - gammaln(a(i)) - t + log(log2(1 + 1./(s(i) * t))) + (a(i) - 1) * log(t));
%         Exp_R_th(i) = sum(y) * error;
%         disp(strcat('Exp_R_th(i) = ',num2str(Exp_R_th(i)))); 
%     end
    
    %% Expected P_p_th is computed numerically
    % Exp_P_p_th = 1./(b_P_p .* (a - 1));
    pl = 10^((-alpha_true)/10);
    bins = 100;
    np = 10^(-noise_power/10);
    [f x_pts] = hist(P_p_sim, bins);
    %x_pts = linspace(min(P_p_sim), max(P_p_sim),1e5);
    % Confidence probability power at PR
    a = (N_th - 2)/4 * log(1./(N_th * (x_pts/np)))  + (N_th + 2)/4 * log( N_th * pl^2 * (x_pts/np) ./ ...
                (pl + 2 * (x_pts/np) + N_th * pl * (x_pts/np)).^2);
            b_d = 4 * N_th * pl^2 * (x_pts/np)./(pl + 2 * (x_pts/np) + N_th * pl * (x_pts/np)).^2;
            b = log(abs(hypergeom([(2 + N_th)/4, (4 + N_th) /4], N_th/2, b_d)));
            pdf_scaled_inv_ncx2 = (x_pts(3) - x_pts(2)) * 1/np * 1./(pl * (x_pts/np)) .* exp(a + b);
    Exp_P_p_th = sum(pdf_scaled_inv_ncx2 .* x_pts);
    figure(1);
    bar(x_pts,f/sum(f));
        hold on,
    plot(x_pts, pdf_scaled_inv_ncx2, 'g', 'Linewidth',2.5);     
        
    Exp_P_p_th = sum(x_pts .* pdf_scaled_inv_ncx2);   
        
    PC_P_p_th = gammainc( 1./((1 + mu_P) .* Exp_P_p_th .* b_P_p), a, 'upper')...
        - gammainc( 1./((1 - mu_P) .* Exp_P_p_th .* b_P_p), a, 'upper');

    % PC rate
    PC_R_th = gammainc( 1./((2.^(K./(K- N_th) * (1 + mu_R).* Exp_R_th) - 1) .* 10^(noise_power/10) .* b_P_s), a, 'upper')...
        -  gammainc( 1./((2.^(K./(K- N_th) * (1 - mu_R).* Exp_R_th) - 1) .* 10^(noise_power/10) .* b_P_s), a, 'upper');


%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %   Finding optimum sensing time
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     N_opt = ceil([1:0.05:30] * (1e-3 * f_s));
%     a = (1 + 10^(SNR_rcvd/10))^2 * N_opt/(2 + 4 * 10^(SNR_rcvd/10));
%     b = 1./a;      
%     b_alpha_true = b*10^(-alpha_true/10);                   % Gamma distributed
%     b_P_reg = b_alpha_true/10^(I_T/10);                     % Inverse Gamma distributed
%     b_P_p = b_P_reg/10^(-alpha_true/10);                    % Inverse Gamma distributed
%     b_P_s = b_P_reg/10^(-alpha_true/10);                    % Inverse Gamma distributed
%     epsilon_R = 0.95;                                       % Probability constraint
%     epsilon_P = 0.95;                                       % Probability constraint
%     N_optimum = 0;
% 
% 
%     for i = length(N_opt):-1:1
%         Exp_P_p_opt = 1./(b_P_p(i) .* (a(i) - 1));
% 
%         % Confidence probability power at PR
%         PC_P_p_opt = gammainc( 1./((1 + gamma_a) .* Exp_P_p_opt .* b_P_p(i)), a(i), 'upper')...
%             - gammainc( 1./((1 - gamma_a) .* Exp_P_p_opt .* b_P_p(i)), a(i), 'upper');
% 
%         % Numerical integration
%         s = 10^(noise_power/10) * b_P_s;
%         error = 0.01;
%         t = ceil(0.01*a(i)):error:ceil(10*a(i));
%         y = (K- N_opt(i))/K * exp( - gammaln(a(i)) - t + log(log2(1 + 1./(s(i) * t))) + (a(i) - 1) * log(t));
%         Exp_R_opt = sum(y) * error;
% 
%         % PC rate
%         PC_R_opt = gammainc( 1./((2.^(K./(K- N_opt(i)) * (1 + gamma_a).* Exp_R_opt) - 1) .* 10^(noise_power/10) .* b_P_s(i)), a(i), 'upper')...
%             -  gammainc( 1./((2.^(K./(K- N_opt(i)) * (1 - gamma_a).* Exp_R_opt) - 1) .* 10^(noise_power/10) .* b_P_s(i)), a(i), 'upper');
% 
%         if (PC_P_p_opt < epsilon_P) | (PC_R_opt < epsilon_R)
%             N_optimum = N_opt(i);        
%             break;
%         end
%         disp(strcat('Exp_R_opt = ',num2str(Exp_R_opt))); 
%         disp(strcat('PC_P_p_opt = ',num2str(PC_P_p_opt))); 
%         disp(strcat('PC_R_opt = ',num2str(PC_R_opt))); 
%     end 
    save('results_thr_est_time_tradeoff_fading_th.mat');
    %quit;
end
load('results_thr_est_time_tradeoff_fading_th.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fontsize = 18;
% PC_opt = min(PC_P_p_opt, PC_R_opt);
% [ax h1 h2] = plotyy(N_th * 1e-3, Exp_R_th, N_th * 1e-3, PC_P_p_th);
% set(h1,'LineWidth',1.5);
% set(h2,'LineWidth',1.5);
% set(ax,'NextPlot','add');
% h3 = plot(ax(2), N_th * 1e-3, PC_R_th,'LineWidth',1.5,'LineStyle','-.');
% set(ax,'NextPlot','add');
% h4 = plot(ax(2), [N_optimum,max(N_opt)] * 1e-3, PC_opt * ones(1, 2),'LineStyle','--','LineWidth',1.5);
% set(ax,'NextPlot','add');
% plot(ax(2), N_optimum * ones(1, 2) * 1e-3, [PC_opt, min(min(PC_P_p_th, min(PC_R_th)))],...
%     'LineWidth',1.5,'LineStyle','--', 'HandleVisibility','off');
% set(ax,'NextPlot','add');
% h5 = plot(ax(1), N_optimum * 1e-3, Exp_R_opt, 'o','LineWidth',1.5);

load('results_thr_est_time_tradeoff_fading_sim.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ax h1 h2] = plotyy(N_sim * 1e-3, Exp_R_sim, N_sim * 1e-3, PC_P_p_sim);
% set(ax,'NextPlot','add');
% h6 = plot(ax(1), N_sim * 1e-3, Exp_R_sim, 'r*', 'HandleVisibility','off');
% set(ax,'NextPlot','add');
% plot(ax(2), N_sim * 1e-3, PC_P_p_sim, 'r*', 'HandleVisibility','off');
% set(ax,'NextPlot','add');
% plot(ax(2), N_sim * 1e-3, PC_R_sim,'r*', 'HandleVisibility','off');

grid on;
axis(ax(1), [min(N_sim * 1e-3) max(N_sim * 1e-3) min(Exp_R_sim) max(Exp_R_sim)]);
axis(ax(2), [min(N_sim * 1e-3) max(N_sim * 1e-3) min(PC_P_p_sim) 1]);
ylabel(ax(1),'R_s = [bits/sec/Hz]','FontSize',Fontsize+2);
ylabel(ax(2),'Probability','FontSize',Fontsize+2);
xlabel(ax(1),'\tau = [ms]','FontSize',Fontsize+2);
xlabel(ax(2), '\tau = [ms]','FontSize',Fontsize+2);
hl = legend([h1;h2;h3;h4;h5;h6],'R', 'PC_P', 'PC_R', 'min(PC_P, PC_R) > \epsilon', 'max \tau', 'sim');
set(hl, 'Position', [.7, .5, .1, 0.2], 'FontSize', Fontsize + 2);
set(ax(1),'FontSize',Fontsize);
set(ax(2),'FontSize',Fontsize);
%Ytick_R =  round(linspace(min(Exp_R_th), max(Exp_R_th), 10)*1000)/1000;
%Ytick_PC = round(linspace(min(min(PC_P_p_th, min(PC_R_th))), 1, 10)*100)/100;
%set(ax(1), 'YTick', Ytick_R);
%set(ax(2), 'YTick', Ytick_PC);
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize',[10 7.5]);
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
print(gcf,'-dpdf','../figures/fig_thr_est_time_tradeoff_fading');
print(gcf,'-depsc','../figures/fig_thr_est_time_tradeoff_fading');
