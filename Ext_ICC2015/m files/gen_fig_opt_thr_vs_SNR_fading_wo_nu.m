%                                
%        Description: This m-file considers the variation of optimum throughput
%        determined from the Estimation-Throughput tradeoff for the fading channel.  
%        (When an outage constraint is established over the power received at PR)      
%
%        Clearly the regulated power is distributed according non central
%        chi-2 distribution. However to simplify the analysis we
%        approximate the non central distribution against the Gamma
%        distribution.
%        
%        The simulation is performed to show the following:
%        1) Analyse the optimum throughput vs received SNR
%
%        Created on: 04.08.15
%        Revision History: 04.08.15 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

sim = 0;                                      % Enable( = 1) to perform simulation, theoretical analysis
th = 0;                                       % disable ( = 0) to plot curves, data is read from a file
                                              % Theoretical anaylsis is also included of simulation as
                                              % numerical integration is involved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_tran = -00;                                 % Power transmitted by PR, the SNR received at ST can be 
                                              % varried using this parameter
noise_power = -100;                           % noise power -100 dBm
I_T = -110;                                   % Interference temperature -80 dBm
f_s = 1e6;                                    % 1 MHz one band
K = 0.1 * f_s;                                % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha_true_p = [115:-1:90];                   % True Path loss between ST and PR  to achieve P_p = IT  
alpha_true_s = 080;                           % True Path loss between ST and SR  to achieve P_p = IT  
epsilon = [0.01 0.1];                         % Outage Probability constraint at PR

Exp_opt_R_th = zeros(length(epsilon),...
    length(alpha_true_p));                    % Optimum Expected throughput  
Opt_N_th = zeros(length(epsilon),...
    length(alpha_true_p));                    % Optimum Expected throughput  
Exp_opt_R_ga_th = zeros(length(epsilon),...
    length(alpha_true_p));                    % Optimum Expected throughput  
Opt_N_ga_th = zeros(length(epsilon),...
    length(alpha_true_p));                    % Optimum Expected throughput  
snr = 10^(P_tran/10) * 10.^(-alpha_true_p/10)...  % snr received at the PR
    / 10^(noise_power/10);
m_p = 1;                                      % Number of realizations of channel for the link ST-PR
m_s = 1;                                      % Number of realizations of channel for the link ST-SR 

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th  
    %% Theoretical parameters
    N_th = [1:1:39,40:30:2980,3000:500:10000];  % N = Total of samples used for estimation = tau * f_s (ms)
    P_reg_th = zeros(1, length(N_th));          % Expected power regulated by ST
    Exp_R_th = zeros(1, length(N_th));          % Expected Throughput
    num_tp = 5000;                              % Number of test points
    P_reg_ga_th = zeros(1,length(N_th));        % Expected power regulated by ST (applying Gamma Approximation)
    Exp_R_ga_th = zeros(1,length(N_th));        % Expected throughput (applying Gamma Approximation) 
    
    for i = 1:length(epsilon)
       disp(strcat('epsilon = ',num2str(epsilon(i))));             

       for j = 1:length(snr)
           disp(strcat('snr = ',num2str(10*log10(snr(j)))));      
           for k=1:length(N_th)       
               if 1 %% Method 2  -- Working But takes to long       
                   M = 1e3;
                   g_p = ones(1, M) .* random('gam',m_p, 1, 1, M);
                   P_rcvd_sim = zeros(1, M);
                   parfor l=1:M
                        % Estimate the received power P_reg with additive noise amplitude (zero mean and noise variance) 
                        samples = sqrt(g_p(l) * 10^((P_tran)/10) * 10^(-alpha_true_p(j)/10)) * ones(1, N_th(k)) + random('norm', 0, sqrt(10^(noise_power/10)), 1, N_th(k));
                        % Checkpoint 
                        % Here the samples must be Gaussian distributed with mean = 10^(P_reg/10), 
                        % examine by ploting the histogram, hist(samples, 100)
                        P_rcvd_sim(l) = mean(samples.^2);                      

                    end 
                    test_points = linspace(10^(noise_power/10) * (1 + m_p * snr(j)), max(P_rcvd_sim), num_tp);
                    for l = num_tp:-1:1                      
                        func_prcvd = @(t) gammainc(test_points(l)./(10^(noise_power/10)/N_th(k)...
                            * (2 + 4 * t * snr(j))./(1 + t * snr(j))),...
                            N_th(k) * (1 + t * snr(j)).^2./(2 + 4 * t * snr(j))) .* exp(-t);
                        result = integral(func_prcvd, 0, 100);                
                        if result <= (1 - epsilon(i))
                            test_points(i);
                            P_reg_th(k) = 10^(I_T/10) * 10^(P_tran/10)/...
                                (test_points(l) - 10^(noise_power/10));                    
                            break;
                        end
                    end
                   %% Expected Rate
                   func_r_em = @(u) (K- N_th(k))/K * log2(1 + u * P_reg_th(k) * 10^(-alpha_true_s/10) /...
                        10^(noise_power/10)) .* exp(-u);
                    Exp_R_th(k) = integral(func_r_em, 0, 100);
                    %disp(strcat('Exp_R_th(k) = ',num2str(Exp_R_th(k)))); 
                    %disp(strcat('P_reg_th(k) = ',num2str(P_reg_th(k))));
               end

               if 1  %% Method 3 -- Approximating the first two moments of the distribution for Prcvd with 
                     %% the moments of Gamma distribution -- Not so Accurate for low SNR and for N > 1000

                    %% The mean and variance are determined from the moment generating function 
                    mean_nc = 10^(noise_power/10) * (1 + m_p * snr(j));
                    var_nc =  (10^(noise_power/10))^2/N_th(k) * (2 * + 4 * m_p *...
                        snr(j) + N_th(k) * m_s * (snr(j))^2);
                    % Parameters for Gamma distribution 
                    b = var_nc / mean_nc;
                    a = mean_nc / b;
                    P_reg_ga_th(k) = 10^(I_T/10) * 10^(P_tran/10) / (b * gammaincinv((1 - epsilon(i)),a) - 10^(noise_power/10));       

                    %% Expected Rate
                    func_r_em = @(u) (K- N_th(k))/K * log2(1 + u * P_reg_ga_th(k) * 10^(-alpha_true_s/10) /...
                        10^(noise_power/10)) .* exp(-u);
                    Exp_R_ga_th(k) = integral(func_r_em, 0, 100);
                    %disp(strcat('Exp_R_ga_th(k) = ',num2str(Exp_R_ga_th(k)))); 
                    %disp(strcat('P_reg_ga_th(k) = ',num2str(P_reg_ga_th(k))));
               end 
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Finding optimum estimation time
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Optimum throughput for the ideal model
            [Exp_opt_R_th(i,j) index] = max(Exp_R_th);
            Opt_N_th(i,j) = N_th(index);            
            [Exp_opt_R_ga_th(i,j) index] = max(Exp_R_ga_th);
            Opt_N_ga_th(i,j) = N_th(index);
            
            disp(strcat('Exp_R_th = ',num2str(Exp_opt_R_th(i,j)))); 
            disp(strcat('Opt_N_th = ',num2str(Opt_N_th(i,j))));
            disp(strcat('Exp_opt_R_ga_th = ',num2str(Exp_opt_R_ga_th(i,j)))); 
            disp(strcat('Opt_N_ga_th = ',num2str(Opt_N_ga_th(i,j))));
       end
   end

   save('results_opt_thr_vs_SNR_fading_e01_th.mat');
   quit;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_opt_thr_vs_SNR_fading_e01_th.mat');
SNR = 10*log10(snr);

%pbaspect([1 0.3 1]);
Fontsize = 9;
%% Ideal Model
if 0    %% Average constraint on the controlled power
    for i = 1:length(SNR)
        func_r_id = @(t,u) log2(1 + t * 10^(I_T/10) * 10^(P_tran/10) * 10^(-alpha_true_s/10)...
            ./ (u  * 10^(-alpha_true_p(i)/10) *  10^(noise_power/10))) .* exp(-t) .* exp(-u);
        Exp_R_id(i) = integral2(func_r_id, 0, 100, 0,100);
    end
end

if 1    %% Outage constraint on the contolled power
    Exp_R_id = zeros(length(epsilon), length(SNR));
    for i = 1:length(epsilon) 
        for j = 1:length(SNR)

               P_reg_id = 10^(I_T/10) / 10^(-alpha_true_p(j)/10) / log(1/(epsilon(i)));

                func_r_id = @(t) log2(1 + P_reg_id * t *  10^(-alpha_true_s/10)...
                    ./ (10^(noise_power/10))) .* exp(-t);
                Exp_R_id(i,j) = integral(func_r_id, 0, 100); 
        end
    end
    hold on,
    h1 = plot(SNR, Exp_R_id(1,:), 'c', 'LineWidth', 2);
    hold on,    
    plot(SNR, Exp_R_id(2,:), 'c', 'LineWidth', 2);    
end


h2 = plot(SNR, Exp_opt_R_th(1,:), 'k', 'LineWidth', 1);
SNR = round(10*log10(snr));
hold on,
h3 = plot(SNR, Exp_opt_R_ga_th(1,:), 'k--', 'LineWidth', 1);
hold on,
h4 = plot(SNR, Exp_opt_R_th(2,:), 'k', 'LineWidth', 1);
%hold on,
%h5 = plot(SNR, Exp_opt_R_th(3,:), 'k', 'LineWidth', 1);
hold on,
h6 = plot(SNR, Exp_opt_R_ga_th(2,:), 'k--', 'LineWidth', 1);
%hold on,
%h7 = plot(SNR, Exp_opt_R_ga_th(3,:), 'k--', 'LineWidth', 1);




grid on;
axis([min(SNR) max(SNR) 0 max(Exp_R_id(2,:)) * 1.01]);
ylabel('$\rs(\ttau)$ [bits/sec/Hz]','FontSize',Fontsize);
xlabel('$\gamma$ [dB]','FontSize',Fontsize);
hl = legend([h1 h2 h3],'IM', 'EM', 'Lemma 2');
%set(hl, 'position',[0.725 0.12 0.15 0.31]);
set(hl, 'position',[0.14 0.12 0.3 0.2]);
set(gca,'FontSize',Fontsize);

if 1
    %%%% Zommed Image %%%%%
    % create a new pair of axes inside current figure
    ax = axes('Position',[.52 .52 .32 .32]);
    box on; % put box around new pair of axes
    child_gca = get(gca,'Children');

    [value indexOfInterest] = find(SNR <= 1 & SNR >= 0);
    index = 2;
    plot(SNR(indexOfInterest), Exp_opt_R_th(index,indexOfInterest), 'k', 'LineWidth', 1); % plot on new axes
    hold on,
    plot(SNR(indexOfInterest), Exp_opt_R_ga_th(index,indexOfInterest), 'k--', 'LineWidth', 1); 
    hold on,
    plot(SNR(indexOfInterest), Exp_R_id(index,indexOfInterest), 'c', 'LineWidth', 1); 
    ylabel(gca, '$\rs(\ttau)$','FontSize',Fontsize);
    xlabel(gca, '$\gamma$','FontSize',Fontsize);
    set(ax, 'XTick', [0 0.5 1]);
    set(ax, 'YTick', [1.8 1.9 2.0]);
    axis(ax, [min(SNR(indexOfInterest)) max(SNR(indexOfInterest))...
        min(Exp_opt_R_ga_th(index,(indexOfInterest)))  max(Exp_R_id(index, indexOfInterest)) * 1]);
    title('Zoom');
end
laprint(1, '../figures/fig_opt_thr_vs_SNR_fading', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Optium estimatio time vs SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
%load('results_opt_thr_vs_SNR_fading_th.mat');
points = [1:2:length(snr),26];
SNR = 10*log10(snr);
h2 = plot(SNR(points), Opt_N_th(1,(points)) * 1e-3, 'k', 'LineWidth', 1);
hold on,
h3 = plot(SNR(points), Opt_N_ga_th(1,points) * 1e-3, 'k--', 'LineWidth', 1);
hold on,
h4 = plot(SNR(points), Opt_N_th(2,points) * 1e-3 , 'k', 'LineWidth', 1);
hold on,
h5 = plot(SNR(points), Opt_N_ga_th(2,points) * 1e-3, 'k--', 'LineWidth', 1);
%hold on,
%h5 = plot(SNR, Opt_N_th(3,:) * 1e-3, 'k', 'LineWidth', 1);
grid on;
axis([min(SNR) max(SNR) 0  max(Opt_N_th(1,(points))) * 1e-3]);
ylabel('$\ttau$ [ms]','FontSize',Fontsize);
xlabel('$\gamma$ [dB]','FontSize',Fontsize);
hl = legend([h1 h2 h3],'IM', 'EM', 'Lemma');
%set(hl, 'position',[0.725 0.12 0.15 0.31]);
set(hl, 'position',[0.68 0.71 0.22 0.2]);
set(gca,'FontSize',Fontsize);
laprint(2, '../figures/fig_opt_est_time_vs_SNR_fading', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');
