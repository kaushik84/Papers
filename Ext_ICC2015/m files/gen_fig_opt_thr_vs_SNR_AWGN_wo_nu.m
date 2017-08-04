%                                
%        Description: This m-file considers the variation of optimum throughput
%        determined from the Estimation-Throughput tradeoff for the AWGN channel.  
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
P_tran = -00;                                     % Power transmitted by PR, the SNR received at ST can be 
                                                  % varried using this parameter
noise_power = -100;                               % noise power -100 dBm
I_T = -110;                                       % Interference temperature -80 dBm
f_s = 1e6;                                        % 1 MHz one band
K = 0.1 * f_s;                                    % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha_true_p = [120:-0.1:90];                     % True Path loss between ST and PR  to achieve P_p = IT  
alpha_true_s = 080;                               % True Path loss between ST and SR  to achieve P_p = IT  
epsilon = [0.01, 0.1];                            % Outage Probability constraint at PR

P_reg_max = 10.^([-10 -00]/10);                   % Maximum Trasmit Power Constraint


Exp_opt_R_th = zeros(length(epsilon),...
    length(alpha_true_p), length(P_reg_max));     % Optimum Expected throughput  
Opt_N_th = zeros(length(epsilon),...
    length(alpha_true_p), length(P_reg_max));     % Optimum Estimation time 
snr = 10^(P_tran/10) * 10.^(-alpha_true_p/10)...  % snr received at the PR
    / 10^(noise_power/10);


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th  
    %% Theoretical parameters
    N_th = 1:1:7000;                               % N = Total of samples used for estimation = tau * f_s 
    P_reg_th = zeros(1, length(N_th));              % Expected power regulated by ST
    Exp_R_th = zeros(1, length(N_th));              % Expected Throughput
    for i = 1:length(epsilon)
       disp(strcat('epsilon = ',num2str(epsilon(i))));             

       for j = 1:length(snr)
           disp(strcat('snr = ',num2str(10*log10(snr(j)))));  
           for k = 1:length(P_reg_max)
               for l=1:length(N_th)       
                   %disp(strcat('N = ',num2str(N_th(k))));             

                   %% Determining the performance meterics --> Exp_R
                   %% Expected values           

                   %% Gamma Approximation to the non-central chi-squared distribution
                   mean = 10^(noise_power/10) * (1 + snr(j));       
                   var = (10^(noise_power/10))^2/N_th(l) * ( 2 + 4 * snr(j)); 

                   b = var/mean;
                   a = mean/b;

                   %% Determine the controlled power
                   P_reg_th(l) = min(P_reg_max(k), 10^(I_T/10) * 10^(P_tran/10) /...
                       (b * gammaincinv((1 - epsilon(i)),a) - 10^(noise_power/10)));       

                   %% Expected Rate
                    Exp_R_th(l) = (K- N_th(l))/K * log2(1 + P_reg_th(l) * 10^(-alpha_true_s/10) /...
                        10^(noise_power/10));
               end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   Finding optimum estimation time
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Optimum throughput for the ideal model
                [Exp_opt_R_th(i,j,k) index] = max(Exp_R_th);
                Opt_N_th(i,j, k) = N_th(index);
                disp(strcat('Exp_R_th = ',num2str(Exp_opt_R_th(i,j, k)))); 
                disp(strcat('N_opt_th = ',num2str(Opt_N_th(i,j, k))));
           end
       end
   end

   save('results_opt_thr_vs_SNR_AWGN_e01_th.mat');
   %quit;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Optium Throughput vs SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
load('results_opt_thr_vs_SNR_AWGN_e01_th.mat');
%pbaspect([1 0.3 1]);
Fontsize = 8.5;
SNR = 10*log10(snr);
h2 = plot(SNR, reshape(Exp_opt_R_th(1,:,1), 1, length(SNR)), 'k-', 'LineWidth', 1);
hold on,
h3 = plot(SNR, reshape(Exp_opt_R_th(2,:,1), 1, length(SNR)), 'k-.', 'LineWidth', 1);
hold on,
plot(SNR, reshape(Exp_opt_R_th(1,:,2), 1, length(SNR)), 'k-', 'LineWidth', 1);
hold on,
plot(SNR, reshape(Exp_opt_R_th(2,:,2), 1, length(SNR)), 'k-.', 'LineWidth', 1);
%hold on,
%h4 = plot(SNR, Exp_opt_R_th(3,:), 'k', 'LineWidth', 1);

%% Ideal Model

Exp_R_id = log2(1 + 10^(I_T/10) * 10^(P_tran/10) * 10^(-alpha_true_s/10)...
    ./ (10.^(-alpha_true_p/10) *  10^(noise_power/10)));
hold on,
h1 = plot(SNR, Exp_R_id, 'c', 'LineWidth', 3);


grid on;
axis([min(SNR) max(SNR) 0  max(Exp_R_id) * 1.01]);
ylabel('$\rs(\ttau)$ [bits/sec/Hz]','FontSize',Fontsize);
xlabel('$\gamma$ [dB]','FontSize',Fontsize);
hl = legend([h1 h2 h3],'IM', 'EM, $\opc = 0.01$', 'EM, $\opc = 0.1$');
%set(hl, 'position',[0.725 0.12 0.15 0.31]);
set(hl, 'position',[0.14 0.12 0.4 0.24]);
set(gca,'FontSize',Fontsize);

if 1
    %%%% Zommed Image %%%%%
    % create a new pair of axes inside current figure
    ax = axes('Position',[.56 .51 .33 .33]);
    box on; % put box around new pair of axes
    child_gca = get(gca,'Children');

    [value indexOfInterest] = find(SNR <= -0.4 & SNR >= -1.4);
    plot(SNR(indexOfInterest), reshape(Exp_opt_R_th(1,indexOfInterest,1), 1, length(SNR(indexOfInterest))),...
    'k', 'LineWidth', 1); % plot on new axes
    hold on,
    plot(SNR(indexOfInterest), reshape(Exp_opt_R_th(2,indexOfInterest,1), 1, length(SNR(indexOfInterest))),...
    'k-.', 'LineWidth', 1); 
    hold on,
    plot(SNR(indexOfInterest), reshape(Exp_opt_R_th(1,indexOfInterest,2), 1, length(SNR(indexOfInterest))),...
    'k-', 'LineWidth', 1); 
    hold on,
    plot(SNR(indexOfInterest), reshape(Exp_opt_R_th(2,indexOfInterest,2), 1, length(SNR(indexOfInterest))),...
    'k-.', 'LineWidth', 1); 
    plot(SNR(indexOfInterest), Exp_R_id(indexOfInterest), 'c', 'LineWidth', 3);

    % plot on new axes
    %plot(SNR(indexOfInterest), Exp_opt_R_th(3,indexOfInterest), 'k', 'LineWidth', 1); 
    ylabel(gca, '$\rs(\ttau)$','FontSize',Fontsize);
    xlabel(gca, '$\gamma$','FontSize',Fontsize);
    %set(ax, 'XTick', [-1.3 -0.9 -0.5]);
    %set(ax, 'YTick', [3.2 3.4 3.7]);
    axis(ax, [min(SNR(indexOfInterest)) max(SNR(indexOfInterest))...
        min(reshape(Exp_opt_R_th(1,indexOfInterest,1), 1, length(SNR(indexOfInterest))))...
        max(Exp_R_id(indexOfInterest)) * 1]);
    set(gca,'FontSize',Fontsize);
    handle = title(ax, 'Zoom');
end
laprint(1, '../figures/fig_opt_thr_vs_SNR_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');

if 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Optium estimatio time vs SNR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2)
    load('results_opt_thr_vs_SNR_AWGN_e01_th.mat');
    SNR = 10*log10(snr);
    h2 = plot(SNR, Opt_N_th(1,:) * 1e-3, 'k', 'LineWidth', 1);
    hold on,
    h3 = plot(SNR, Opt_N_th(2,:) * 1e-3 , 'k', 'LineWidth', 1);
    %hold on,
    %h4 = plot(SNR, Opt_N_th(3,:) * 1e-3, 'k', 'LineWidth', 1);
    grid on;
    axis([min(SNR) max(SNR) min(Opt_N_th(3,:)) * 1e-3  max(Opt_N_th(1,:)) * 1e-3]);
    ylabel('$\ttau$ [ms]','FontSize',Fontsize);
    xlabel('$\gamma$ [dB]','FontSize',Fontsize);
    %hl = legend([h2 h3 h4],'IM','EM');
    %set(hl, 'position',[0.725 0.12 0.15 0.31]);
    set(hl, 'position',[0.68 0.71 0.22 0.2]);
    set(gca,'FontSize',Fontsize);
    laprint(2, '../figures/fig_opt_est_time_vs_SNR_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'off');
end