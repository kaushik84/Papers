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
%        1) Analyse the optimum Controlled Power vs Outage Probability
%
%        Created on: 04.08.15
%        Revision History: 04.08.15 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

sim = 0;                                          % Enable( = 1) to perform simulation, theoretical analysis
th = 1;                                           % disable ( = 0) to plot curves, data is read from a file
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
alpha_true_p = [90 100 110];                         % True Path loss between ST and PR  to achieve P_p = IT  
alpha_true_s = 080;                               % True Path loss between ST and SR  to achieve P_p = IT  
epsilon = [0.01:0.04:.4];                         % Outage Probability constraint at PR

Opt_cont_Power_th = zeros(length(epsilon),...
    length(alpha_true_p));                        % Optimum Contolled Power at the ST 

Exp_opt_R_th = zeros(length(epsilon),...
    length(alpha_true_p));                        % Optimum Expected throughput

Opt_N_th = zeros(length(epsilon),...
    length(alpha_true_p));                        % Optimum Expected throughput  
snr = 10^(P_tran/10) * 10.^(-alpha_true_p/10)...  % snr received at the PR
    / 10^(noise_power/10);
P_reg_max = 10^(-00/10);                          % Maximum Trasmit Power Constraint

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th  
    %% Theoretical parameters
    N_th = ceil([0.001:0.001:5] * (1e-3 * f_s));    % N = Total of samples used for estimation = tau * f_s 
    P_reg_th = zeros(1, length(N_th));              % Expected power regulated by ST
    Exp_R_th = zeros(1, length(N_th));              % Expected Throughput
    for i = 1:length(epsilon)
       disp(strcat('epsilon = ',num2str(epsilon(i))));             

       for j = 1:length(snr)
           disp(strcat('snr = ',num2str(10*log10(snr(j)))));      
           for k=1:length(N_th)       
               %disp(strcat('N = ',num2str(N_th(k))));             

               %% Determining the performance meterics --> Exp_R
               %% Expected values           

               %% Gamma Approximation to the non-central chi-squared distribution
               mean = 10^(noise_power/10) * (1 + snr(j));       
               var = (10^(noise_power/10))^2/N_th(k) * ( 2 + 4 * snr(j)); 

               b = var/mean;
               a = mean/b;

               %% Determine the controlled power
               P_reg_th(k) = min(P_reg_max, 10^(I_T/10) * 10^(P_tran/10) /...
                   (b * gammaincinv((1 - epsilon(i)),a) - 10^(noise_power/10)));       

               %% Expected Rate
                Exp_R_th(k) = (K- N_th(k))/K * log2(1 + P_reg_th(k) * 10^(-alpha_true_s/10) /...
                    10^(noise_power/10));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Finding optimum estimation time
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Optimum throughput for the ideal model
            [Exp_opt_R_th(i,j) index] = max(Exp_R_th);
            Opt_cont_Power_th(i,j) = P_reg_th(index);            
            Opt_N_th(i,j) = N_th(index);
            disp(strcat('Opt_cont_Power_th = ',num2str(Opt_cont_Power_th(i,j)))); 
            disp(strcat('Exp_opt_R_th = ',num2str(Exp_opt_R_th(i,j)))); 
            disp(strcat('N_opt_th = ',num2str(Opt_N_th(i,j))));
       end
   end

   save('results_opt_contPower_vs_pout_th.mat');
   %quit;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Optium Throughput vs SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
load('results_opt_contPower_vs_pout_th.mat');
%pbaspect([1 0.3 1]);
Fontsize = 9;
SNR = 10*log10(snr);
h2 = plot(epsilon, 10* log10(Opt_cont_Power_th(:,1)), 'k', 'LineWidth', 1);
hold on,
h3 = plot(epsilon, 10* log10(Opt_cont_Power_th(:,2)), 'k', 'LineWidth', 1);
hold on,
h4 = plot(epsilon, 10* log10(Opt_cont_Power_th(:,3)), 'k', 'LineWidth', 1);
axis tight;
%% Ideal Model

Opt_cont_Power_id = 10^(I_T/10) * 10^(P_tran/10) ./ (10.^(-alpha_true_p/10));
hold on,
h1 = plot(epsilon, Opt_cont_Power_id * ones(1, length(epsilon)), 'c', 'LineWidth', 2);


grid on;
axis([min(SNR) max(SNR) 0  max(Exp_R_id) * 1.01]);
ylabel('$\rs(\ttau)$ [bits/sec/Hz]','FontSize',Fontsize);
xlabel('$\gamma$ [dB]','FontSize',Fontsize);
hl = legend([h1 h2 h3],'IM', 'EM');
%set(hl, 'position',[0.725 0.12 0.15 0.31]);
set(hl, 'position',[0.68 0.71 0.22 0.2]);
set(gca,'FontSize',Fontsize);

laprint(1, '../figures/fig_opt_contPower_vs_pout_th', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');