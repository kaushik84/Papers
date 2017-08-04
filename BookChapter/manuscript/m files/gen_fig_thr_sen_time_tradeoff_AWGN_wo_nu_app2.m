%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: This m-file considers the Throughput-Sensing time 
%        tradeoff with SNR estimation for the AWGN channel. In this m-script  
%        the short term analysis is investigted that is, considering invidual 
%        frame. Hence, the SNR is considered constant for the time interval 
%        (tau_est + tau_sen). 
%        
%        The simulation is performed to demenstrate the performance of the 
%        second approach to the optimization problem as stated by the
%        reviewers that the considers the fluctuations in the optimal
%        tau_sen and determine its average value for the OC. Unlike other
%        pervious approach where a average value of the throughput is used
%        for finding optimum tau.
%        
%        The simulation is performed to show the following:
%        1) Analyse the throughput vs sensing for different cases 
%           - Ideal case, when SNR is known at the ST (i)
%           - Average Constarint
%           - Outage Constraint   
%        2) Confirm the theoretical analysis vs simulation results
%
%        Created on: 30.10.15
%        Last modified: 30.10.15
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
%warning off;


sim = 00;                                   % Enable( = 1) to perform simulation, theoretical analysis
th = 00;                                    % disable ( = 0) to plot curves, data is read from a file
                                            % Theoretical anaylsis is also included of simulation as
                                            % numerical integration is involved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_s = 10^(-10/10);                          % Power transmitted by ST, the SNR received at SR 
P_p = 10^(-10/10);                          % Power transmitted by PT, the SNR received at SR 
                                            % varried using this parameter
noise_power = 10^(-100/10);                 % noise power -100 dBm
f_s = 1e6;                                  % 1 MHz one band
K = 0.1 * f_s;                              % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha_p_1 = 10^(-100/10);                   % True Path loss between ST and PR   
alpha_p_2 = 10^(-100/10);                   % True Path loss between PR and SR   
alpha_s = 10^(-080/10);                     % True Path loss between ST and SR 
P_H0 = 0.8;                                 % Probability of Hypothesis 0
P_d_d = 0.90;                               % Constraint on Probability of detection P_d
mu = 0.05;                                  % Outage Probability on probability of detection

num_test_points = 10000;                    % Test points for evaluating threshold 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sim
     M = 5e4;                                           % Number of realizations 
     %% Simulations parameters
     N_est_sim = 5000;                                  % N_est_th = number of samples used for estimation = tau * f_s 

     N_sen_sim = [N_est_sim:10:12000];                 % N = Total of samples used for sensing = tau * f_s 

     epsilon_oc_sim = zeros(1,length(M));               % Inaccurate threshold due to incorrect estimation of received power

     P_f_oc_sim = zeros(1,length(N_sen_sim));           % Probability of false alarm at ST due to inaccurate estimation of threshold

     P_d_oc_sim = zeros(1,length(N_sen_sim));           % Probability of detection at ST due to accurate estimation of threshold
    

     P_rcvd_est_sen_sim = zeros(1,M);                    % Rcvd power samples under Hypothesis 0 (sensing)
     P_rcvd_est_int_sim = zeros(1,M);                    % Rcvd power samples under Hypothesis 1 (estimation) 
     
     N_s = 10;                                           % Number of Pilot Symobols 
     sigma_est_error = noise_power/N_s;                  % Variance of the Estimation error for h_s 
     dof_s = 1;                                          % Degree of freedom = 1
     g_s_hat = zeros(1,M);                               % Estimated power gain for the channel g_s        
        

     opt_tau_sen = zeros(1,M);                          % A array of optimum tau corresponsing to different variations

     %% Estimation of the sensing channel 
     parfor i=1:M
         P_rcvd_est_sen_sim(i) = mean(random('norm',...
         random('norm',0, sqrt(P_p * alpha_p_1), 1, N_est_sim),...
         sqrt(noise_power), 1, N_est_sim).^2);                               
     end  

     %% Estimation of the access channel
     g_s = P_s * alpha_s;
     parfor i=1:M
         g_s_hat(i) = sum(normrnd(sqrt(g_s), sqrt(sigma_est_error), 1, dof_s).^2);
     end

     %% Estimation of the interfererence channel
     parfor i=1:M
         P_rcvd_est_int_sim(i) = mean(random('norm',...
         random('norm',0, sqrt(P_p * alpha_p_2), 1, N_est_sim),...
         sqrt(noise_power), 1, N_est_sim).^2);                               
     end  
     
     
     Acc_energy = (P_p * alpha_p_1 + noise_power);
     
     %% Determine the variations in the performance parameters C_0, C_1, beacuse the performance parameters 
     %% P_d, P_fa depends on sensing time hence they are calculated inside the loop
     C_0 =  log2(1 + g_s_hat/noise_power);
     C_1 =  log2(1 + g_s_hat./P_rcvd_est_int_sim);
     parfor i=1:M        
        P_d_oc_sim = zeros(1, length(N_sen_sim));           % Probability of detection at ST due to inaccurate estimation of threshold
        %R_oc_sim = zeros(1,length(N_sen_sim));             % Throughput at SR with upper limit of the energy wall
            
        disp(strcat('iteration = ',num2str(i))); 

        
        %% Outage constraint            
        %% Threshold                  
        epsilon_oc_th = 4 * Acc_energy * gammaincinv(1 - mu, N_est_sim/2, 'upper') * ...
            gammaincinv(P_d_d, N_sen_sim/2, 'upper') ./ (N_est_sim * N_sen_sim); 

        %% Probability of false alarm and probability of detection
        P_f_oc_sim = gammainc(N_sen_sim/2 .* epsilon_oc_th/noise_power ,N_sen_sim/2, 'upper');
        P_d_oc_sim = gammainc(N_sen_sim/2 .* epsilon_oc_th/P_rcvd_est_sen_sim(i) ,N_sen_sim/2, 'upper');

        [values index] = max((K - N_sen_sim)/K .* (P_H0 * (1 -  P_f_oc_sim) *...
            C_0(i) +  (1 - P_H0) * (1 - P_d_oc_sim) * C_1(i)));   

        opt_tau_sen(i) = N_sen_sim(index);
        
        disp(strcat('opt_tau_sen(i) = ',num2str(opt_tau_sen(i) * 1e-3)));   
        
     end
    avg_opt_tau_sen = (opt_tau_sen);
    save('results_thr_sen_time_tradeoff_AWGN_wo_nu_snr_m10_sim2_app2.mat');
    if ~th      % If theoretical analysis is already simulated the computation ends here.
        quit;
    end
end   

