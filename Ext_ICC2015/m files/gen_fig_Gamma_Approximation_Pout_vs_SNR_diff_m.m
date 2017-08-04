%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: In this m-script  
%        the long term analysis (channel is average over different channel realization) 
%        is investigted that is, considering individual 
%        frame. Hence, the SNR is considered constant for the time interval 
%        (tau_est). 
%
%        Clearly the recieved power is distributed according to non central
%        chi-2 distribution where the SNR follows an exponentially or Gamma distribution,
%        Analytical form is of the obtained distribution does not exist or
%        difficult to evaluate.
%        We use the moment generating function to evaluate the first two moments. 
%        Finally we approximate the obtained distribuion with as Gamma distribution.
%        Moment matching is used to obtain the parameters of the Gamma
%        distribution.
%        
%        The simulation is performed to show the following:
%        1) Confirm the theoretical analysis (approximation) vs simulation
%        results, and state whether the approximation is good. For
%        performance analysis we consider variation of P_out against SNR
%
%        Created on: 22.12.14
%        Last modified: 22.12.14
%        Revision History: 22.12.14 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

sim = 1;                                         % Enable( = 1) to perform simulation, theoretical analysis
                                                  % disable ( = 0) to plot curves, data is read from a file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 1e4;                                          % Number of realizations 
P = 10.^(-05/10);%(-10:0.5:-10)/10);                        % Power transmitted by SR, the SNR received at SR 
                                                  % varried using this parameter
noise_power = 10^(-100/10);                       % noise power -100 dBm
f_s = 1e6;                                        % 1 MHz one band
K = 0.1 * f_s;                                    % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha = 10^(-100/10);                             % True Path loss between ST and PR   
m = [2];                                    % Number of realizations of the channel, Nakagami-m parameter 
P_rcvd_sim =  zeros(1,M);                         % Received Power  
N = 1000;

threshold_sim = zeros(length(m), length(P));      % Threshold attained for a certain P_d_d simulated 
threshold_th = zeros(length(m), length(P));       % Threshold attained for a certain P_d_d theretical 

P_out_sim = zeros(length(m), length(P));          % False alarm determined from the simulated threshold
P_out_th = zeros(length(m), length(P));           % False alarm determined from the theoretical threshold 
P_d_d = 0.9;                                      % Target Probability of Detection
epsilon = 0.1;                                   % Epsilon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sim
     %% Theoretical parameters
     for i=1:length(m)
         disp(strcat('m = ',num2str(m(i)))); 
         for j=1:length(P)
            g = ones(1, M) .* random('gam', m(i), 1, 1, M); 
            disp(strcat('P = ',num2str(P(j)))); 
            for k=1:M
                % Estimate the received power P_reg with additive noise amplitude (zero mean and noise variance) 
                samples = sqrt(g(k) * P(j) * alpha) * ones(1, N) + random('norm', 0, sqrt(noise_power), 1, N);
                % Checkpoint 
                % Here the samples must be Gaussian distributed with mean = 10^(P_reg/10), 
                % examine by ploting the histogram, hist(samples, 100)
                P_rcvd_sim(k) = mean(samples.^2);                                  
            end 

            %% Gamma Approxmiamtion
            mean_nc = (noise_power + m(i) * alpha * P(j));
            var_nc =  1/N * (2 * noise_power^2 + 4 * m(i) *...
                alpha * P(j) * noise_power + N * m(i) * (alpha * P(j))^2);
            % Parameters for Gamma distribution 
            b = var_nc / mean_nc;
            a = mean_nc / b;

            %% Simulated CDF
            [CDF_Power_sim, Power,~,~,eid] = cdfcalc(P_rcvd_sim);
            CDF_Power_sim = CDF_Power_sim(1:length(Power));
            
            % Plot CDF of the simulated and analytical Expressions
            figure(1);
            plot(Power, CDF_Power_sim);
            hold on
            plot(Power,gammainc(Power/b, a), 'r');
            
            % Calculate the threshold where detection probability = P_d_d
            [index value] = find(CDF_Power_sim >=  (1 - epsilon));
            threshold_sim(i,j) = Power(min(index));
            
            % Thoeretical Analysis 
            %CDF_Power_th = 1 - gammainc(a, Power/b);     
            threshold_th(i,j) = b * gammaincinv(1- epsilon, a); 
            
            %% From the two cases the P_fa is evaluated 
            P_out_sim(i,j) = gammainc(threshold_sim(i,j) * N /noise_power/2, N/2);
            P_out_th(i,j) = gammainc(threshold_th(i,j) * N /noise_power/2, N/2);

         end
     end
    save('results_Gamma_Approximation_Pout_vs_SNR_diff_m_N10.mat');
    %quit;
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_Gamma_Approximation_Pout_vs_SNR_diff_m_N10.mat');
Fontsize = 18;
%% semilog scale P_fa
if 1
    semilogy(10 *log10(P), P_fa_sim(1,:),'rs', 'Linewidth',2);
    hold on,
    semilogy(10 *log10(P), P_fa_th(1,:),'--', 'Linewidth',2);
    hold on,
    semilogy(10 *log10(P), P_fa_sim(2,:),'rs', 'Linewidth',2);
    hold on,
    semilogy(10 *log10(P), P_fa_th(2,:),'--', 'Linewidth',2);
    hold on,
    semilogy(10 *log10(P), P_fa_sim(3,:),'rs', 'Linewidth',2);
    hold on,
    semilogy(10 *log10(P), P_fa_th(3,:),'--', 'Linewidth',2);
    hold on,
    semilogy(10 *log10(P), P_fa_sim(4,:),'rs', 'Linewidth',2);
    hold on,
    semilogy(10 *log10(P), P_fa_th(4,:),'--', 'Linewidth',2);
    %hold on,
    %semilogy(10 *log10(P), P_fa_sim(5,:),'k', 'Linewidth',2);
    %hold on,
    %semilogy(10 *log10(P), P_fa_th(5,:),'k--', 'Linewidth',2);
    axis([min(10 * log10(P)) 10 1e-3 1e0]);
    ylabel('P_{fa}','FontSize',Fontsize+2);
    xlabel('SNR = [dB]','FontSize',Fontsize+2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves cosmetic makeover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid on;
%axis tight;

hl = legend('Simulated', 'Approximation');
set(hl, 'Location', 'SouthWest', 'FontSize', Fontsize + 2);
set(gca,'FontSize',Fontsize);
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize',[10 7.5]);
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
print(gcf,'-dpdf','../figures/fig_Gamma_Approximation_Pfa_vs_SNR_diff_m');
print(gcf,'-depsc','../figures/fig_Gamma_Approximation_Pfa_vs_SNR_diff_m');