%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: This m-file considers the Throughput-Sesning time 
%        tradeoff with SNR estimation for the AWGN channel. In this m-script  
%        the short term analysis is investigted that is, considering invidual 
%        frame. Hence, the SNR is considered constant for the time interval 
%        (tau_est + tau_sen). 
%
%        Clearly the recieved power is distributed according to non central
%        chi-2 distribution, which can be approximated as Gaussian 
%        distribution for large sample count (This again needs an analysis, 
%        on how good is the approximation is?)
%        
%        The theortical analysis is performed to investigate the influence of 
%        SNR on the optimum throughput for different cases 
%           - Ideal case, when SNR is known at the ST (i)
%           - Average Constarint
%           - Outage Constraint   
%
%        Created on: 11.03.15
%        Last modified: 11.03.15
%        Revision History: 11.03.15 --> File generated   
%                          01.04.15 --> Optimum P_fa and Optimum P_d for the 
%                                       corresponding constraints added
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
warning off;


sim = 00;                                 % Enable( = 1) to perform simulation, theoretical analysis
th = 00;                                  % disable ( = 0) to plot curves, data is read from a file
                                          % Theoretical anaylsis is also included of simulation as
                                          % numerical integration is involved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_s = 10^(-10/10);                        % Power transmitted by PR, the SNR received at ST 
P_p = 10.^((-13.5:0.5:10)/10);            % Power transmitted by ST, the SNR received at SR 
                                          % varried using this parameter
noise_power = 10^(-100/10);               % noise power -100 dBm
f_s = 1e6;                                % 1 MHz one band
K = 0.1 * f_s;                            % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha_p_1 = 10^(-100/10);                 % True Path loss between ST and PR   
alpha_p_2 = 10^(-100/10);                 % True Path loss between PR and SR   
alpha_s = 10^(-080/10);                   % True Path loss between ST and SR 
P_H0 = 0.8;                               % Probability of Hypothesis 0
P_d_d = 0.90;                             % Constraint on Probability of detection P_d
mu = 0.05;                                % Outage Probability on probability of detection

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th
   %% Theoretical parameters
    N_est_th = 5 * 1000;                               % N_est_th = number of samples used for estimation = tau * f_s 
    N_th = [0.1:0.05:50] * 1000;                       % N = includes the time in samples, for the ideal case their is no
                                                       % estimation time hence sensing time = N_th and the throughput will 
                                                       % attain a nonzero value. For the accurate, the energy walls sensing 
                                                       % time begin after the estimation of the rcvd energy is performed.
                                                       % N_sen = N_th - N_est
    P_f_id_th = zeros(length(P_p), length(N_th));      % Prob. of false alarm at ST @ optimum sensing time for ideal case 
    P_f_ac_th = zeros(length(P_p), length(N_th));      % Prob. of false alarm at ST @ optimum sensing time with average constraint 
    P_f_oc_th = zeros(length(P_p), length(N_th));      % Prob. of false alarm at ST @ optimum sensing time with outage constraint 
    
    P_d_id_th = zeros(length(P_p), length(N_th));      % Prob. of detection at ST @ optimum sensing time for ideal case 
    P_d_ac_th = zeros(length(P_p), length(N_th));      % Prob. of detection at ST @ optimum sensing time with average constraint 
    P_d_oc_th = zeros(length(P_p), length(N_th));      % Prob. of detection at ST @ optimum sensing time with outage constraint 
                                                       
    
    R_id_th = zeros(length(P_p),length(N_th));         % Throughput at SR for ideal case
    R_ac_th = zeros(length(P_p),length(N_th));         % Throughput at SR with average constraint
    R_oc_th = zeros(length(P_p),length(N_th));         % Throughput at SR with outage constraint
    
    R_id_opt_th = zeros(1, length(P_p));               % Throughput at SR @ optimum sensing time for ideal case
    R_ac_opt_th = zeros(1, length(P_p));               % Throughput at SR @ optimum sensing time with average constraint
    R_oc_opt_th = zeros(1, length(P_p));               % Throughput at SR @ optimum sensing time with outage constraint
    
    P_f_id_opt_th = zeros(1, length(P_p));             % Prob. of false alarm at ST @ optimum sensing time for ideal case 
    P_f_ac_opt_th = zeros(1, length(P_p));             % Prob. of false alarm at ST @ optimum sensing time with average constraint 
    P_f_oc_opt_th = zeros(1, length(P_p));             % Prob. of false alarm at ST @ optimum sensing time with outage constraint 
    
    P_d_id_opt_th = zeros(1, length(P_p));             % Prob. of detection at ST @ optimum sensing time for ideal case 
    P_d_ac_opt_th = zeros(1, length(P_p));             % Prob. of detection at ST @ optimum sensing time with average constraint 
    P_d_oc_opt_th = zeros(1, length(P_p));             % Prob. of detection at ST @ optimum sensing time with outage constraint 

    for i = 1:length(P_p) 
        
        disp(strcat('P_p(i) = ',num2str(P_p(i))));
        
        % Accurarte energy  
        Acc_energy = (noise_power + alpha_p_1 * P_p(i));      

        for j=1:length(N_th)
            warning off;
            %disp(strcat('N = ',num2str(N_th(j))));     

            %% Determining the performance meteric R for all the different cases
            %% Ideal         
            %% Threshold  
            epsilon_id_th = qfuncinv(P_d_d) * sqrt(2/N_th(j)) * Acc_energy +...
                Acc_energy;                  

            %% Probability of false alarm and probability of dectection
            P_f_id_th(i,j) = qfunc((epsilon_id_th - noise_power)/...
                (sqrt(2/N_th(j)) * noise_power));        

            %% Expected Rate
            C_0_th = log2(1 + alpha_s * P_s / noise_power); 
            C_1_th = log2(1 + alpha_s * P_s / (P_p(i) * alpha_p_2 +...
                noise_power) );
            P_d_id_th(i,j) =  P_d_d;

            R_id_th(i,j)  = (K - N_th(j))/K * (P_H0 * (1 -  P_f_id_th(i,j)) * C_0_th +...
                (1 - P_H0) * (1 - P_d_id_th(i,j)) * C_1_th);  


            if N_th(j) > N_est_th
               N_sen_th = N_th(j) - N_est_th;
               %% Average case
               if 1
                   %% Threshold  
                   num_test_points = 1000;
                   test_points = linspace(Acc_energy* (1 + sqrt(2/N_sen_th)* qfuncinv(0.75)), ...
                       Acc_energy* (1 + sqrt(2/N_sen_th)* qfuncinv(0.9999)),num_test_points);
                   Exp_P_d_ac = zeros(1,num_test_points);
                   for k=1:num_test_points
                        ts = sqrt(2/N_sen_th);
                        te = sqrt(2/N_est_th);
                        func_exp_P_d = @(t)  t .*  test_points(k) .*...
                            exp( (qfuncinv(t)/sqrt(2)).^2 - (test_points(k)./(1 + ts * qfuncinv(t)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
                            * ts ./ (Acc_energy * te * ( 1 + ts * qfuncinv(t)).^2);
                                           Exp_P_d_ac(k) =  integral(func_exp_P_d,0,1);
                        if Exp_P_d_ac(k) >= P_d_d
                           epsilon_ac_th = test_points(k); 
                           break;
                        else
                           epsilon_ac_th = test_points(k); 
                        end
                   end    

                   %% Probability of false alarm and probability of dectection
                   P_f_ac_th(i,j) = qfunc((epsilon_ac_th - noise_power)/...
                       (sqrt(2/N_sen_th) * noise_power));        

                   P_d_ac_th(i,j) = P_d_d;


                   %% Expected Rate
                   C_0_th = log2(1 + alpha_s * P_s / noise_power); 
                   C_1_th = log2(1 + alpha_s * P_s / (P_p(i) * alpha_p_2 +...
                       noise_power) );

                   R_ac_th(i,j) = (K - N_th(j))/K * (P_H0 * (1 -  P_f_ac_th(i,j)) * C_0_th +...
                       (1 - P_H0) * (1 - P_d_ac_th(i,j)) * C_1_th); 
               end


                %% Outage case
                if 1
                    %% Threshold  
                    epsilon_oc_th = Acc_energy * (1 + qfuncinv(1 - mu) * (sqrt(2/N_est_th)))...
                         * (1 + qfuncinv(P_d_d) * (sqrt(2/N_sen_th)));            

                    %% Probability of false alarm and  expected probability of dectection
                    P_f_oc_th(i,j) = qfunc((epsilon_oc_th - noise_power)/...
                        (sqrt(2/N_sen_th) * noise_power));   

                    ts = sqrt(2/N_sen_th);
                    te = sqrt(2/N_est_th);
                    func_exp_P_d = @(t)  t .*  epsilon_oc_th .*...
                        exp( (qfuncinv(t)/sqrt(2)).^2 - (epsilon_oc_th./(1 + ts * qfuncinv(t)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
                        * ts ./ (Acc_energy * te * ( 1 + ts * qfuncinv(t)).^2);
                    P_d_oc_th(i,j) =  integral(func_exp_P_d,0,1); 

                    %% Expected Rate
                    C_0_th = log2(1 + alpha_s * P_s / noise_power); 
                    C_1_th = log2(1 + alpha_s * P_s / (P_p(i) * alpha_p_2 +...
                        noise_power) );

                    R_oc_th(i,j) = (K - N_th(j))/K * (P_H0 * (1 -  P_f_oc_th(i,j)) * C_0_th +...
                    (1 - P_H0) * (1 - P_d_oc_th(i,j)) * C_1_th);
                end
            end            
        end

        % Throughput, prob. of false alarm, prob. of detection 
        % @ optimum sensing time for the ideal model
        [temp_id index_id] = max(R_id_th(i,:));
        R_id_opt_th(i) = R_id_th(i,index_id);
        P_f_id_opt_th(i) = P_f_id_th(i,index_id);
        P_d_id_opt_th(i) = P_d_id_th(i,index_id);

        % Throughput, prob. of false alarm, prob. of detection 
        % @ optimum sensing time for the for the estimation model -- average constraint case
        [temp_ac index_ac] = max(R_ac_th(i,:));
        R_ac_opt_th(i) = R_ac_th(i,index_ac);
        P_f_ac_opt_th(i) = P_f_ac_th(i,index_ac);
        P_d_ac_opt_th(i) = P_d_ac_th(i,index_ac);
        
        % Throughput, prob. of false alarm, prob. of detection 
        % @ optimum sensing time for the for the estimation model -- average constraint case
        [temp_oc index_oc] = max(R_oc_th(i,:));
        R_oc_opt_th(i) = R_oc_th(i,index_oc);
        P_f_oc_opt_th(i) = P_f_oc_th(i,index_oc);
        P_d_oc_opt_th(i) = P_d_oc_th(i,index_oc);   
    end
    save('results_opt_thr_vs_SNR_AWGN_wo_nu_th.mat');
    quit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Throughput vs SNR @ optimum sensing time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_opt_thr_vs_SNR_AWGN_wo_nu_th.mat');
figure(1);
Fontsize = 9;
plot(10 *log10(P_p), R_id_opt_th, 'b', 'LineWidth',1.5);
hold on,
plot(10 *log10(P_p), R_ac_opt_th, 'LineStyle', '--', 'LineWidth', 1.5);
hold on,
plot(10 *log10(P_p), R_oc_opt_th, 'LineStyle', '-.', 'LineWidth', 1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves cosmetic makeover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid on;
axis([min(10 *log10(P_p)) max(10 *log10(P_p))...
    min(R_oc_opt_th)  max(R_id_opt_th) * 1.015]);
ylabel('$\trs $ [bits/sec/Hz]','FontSize',Fontsize);
xlabel('$\snrrcvd$ [dB]','FontSize',Fontsize);
hl = legend('$\trs$', '$\trsac$', '$\trsoc$');
set(hl, 'Location', 'SouthEast', 'FontSize', Fontsize);
set(gca,'FontSize',Fontsize);
laprint(1, '../figures/fig_opt_thr_vs_SNR_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Probability of false alarm vs SNR @ optimum sensing time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
Fontsize = 9;
plot(10 *log10(P_p), P_f_id_opt_th, 'b', 'LineWidth',1.5);
hold on,
plot(10 *log10(P_p), P_f_ac_opt_th, 'LineStyle', '--', 'LineWidth', 1.5);
hold on,
plot(10 *log10(P_p), P_f_oc_opt_th, 'LineStyle', '-.', 'LineWidth', 1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves cosmetic makeover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid on;
axis([min(10 *log10(P_p)) max(10 *log10(P_p))...
    min(P_f_oc_opt_th)  max(P_f_oc_opt_th) * 1.015]);
ylabel('$\pd$','FontSize',Fontsize);
xlabel('$\snrrcvd$ [dB]','FontSize',Fontsize);
hl = legend('$\trs$', '$\trsac$', '$\trsoc$');
set(hl, 'Location', 'SouthEast', 'FontSize', Fontsize);
set(gca,'FontSize',Fontsize);
laprint(2, '../figures/fig_opt_P_f_vs_SNR_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');
