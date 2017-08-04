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
%        Clearly the recieved power is distributed according to central
%        chi-2 distribution, 
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
P_p = 10.^((-15:0.3:10)/10);             % Power transmitted by ST, the SNR received at SR 
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
num_test_points = 5000;                  % Test points for evaluating threshold 

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th
   %% Theoretical parameters
    N_est_th = 5 * 1000;                                   % N_est_th = number of samples used for estimation = tau * f_s 
    N_sen_th = [0.1:0.05:75] * 1000;                       % N_sen_th = includes the time in samples, for the ideal case their is no
                                                           % estimation time hence sensing time = N_th and the throughput will 
                                                           % attain a nonzero value. For the accurate, the energy walls sensing 
                                                           % time begin after the estimation of the rcvd energy is performed.
                                                           % N_sen = N_th - N_est
    P_f_id_th = zeros(length(P_p), length(N_sen_th));      % Prob. of false alarm at ST @ optimum sensing time for ideal case 
    P_f_ac_th = zeros(length(P_p), length(N_sen_th));      % Prob. of false alarm at ST @ optimum sensing time with average constraint 
    P_f_oc_th = zeros(length(P_p), length(N_sen_th));      % Prob. of false alarm at ST @ optimum sensing time with outage constraint 
    
    P_d_id_th = zeros(length(P_p), length(N_sen_th));      % Prob. of detection at ST @ optimum sensing time for ideal case 
    P_d_ac_th = zeros(length(P_p), length(N_sen_th));      % Prob. of detection at ST @ optimum sensing time with average constraint 
    P_d_oc_th = zeros(length(P_p), length(N_sen_th));      % Prob. of detection at ST @ optimum sensing time with outage constraint 
                                                       
    
    R_id_th = zeros(length(P_p),length(N_sen_th));         % Throughput at SR for ideal case
    R_ac_th = zeros(length(P_p),length(N_sen_th));         % Throughput at SR with average constraint
    R_oc_th = zeros(length(P_p),length(N_sen_th));         % Throughput at SR with outage constraint
    
    R_id_opt_th = zeros(1, length(P_p));                   % Throughput at SR @ optimum sensing time for ideal case
    R_ac_opt_th = zeros(1, length(P_p));                   % Throughput at SR @ optimum sensing time with average constraint
    R_oc_opt_th = zeros(1, length(P_p));                   % Throughput at SR @ optimum sensing time with outage constraint
    
    P_f_id_opt_th = zeros(1, length(P_p));                 % Prob. of false alarm at ST @ optimum sensing time for ideal case 
    P_f_ac_opt_th = zeros(1, length(P_p));                 % Prob. of false alarm at ST @ optimum sensing time with average constraint 
    P_f_oc_opt_th = zeros(1, length(P_p));                 % Prob. of false alarm at ST @ optimum sensing time with outage constraint 
    
    P_d_id_opt_th = zeros(1, length(P_p));                 % Prob. of detection at ST @ optimum sensing time for ideal case 
    P_d_ac_opt_th = zeros(1, length(P_p));                 % Prob. of detection at ST @ optimum sensing time with average constraint 
    P_d_oc_opt_th = zeros(1, length(P_p));                 % Prob. of detection at ST @ optimum sensing time with outage constraint 

    
    for i = 1:length(P_p) 
        
        disp(strcat('P_p(i) = ',num2str(10 * log10(P_p(i)))));
        
        [C_0, C_1] = calc_capacities(alpha_p_2 * P_p(i) / noise_power,... 
            alpha_s * P_s / noise_power,....
            noise_power);                               % Capacities at SR, with and without intereference from PT  
         
        % Accurarte energy  
        Acc_energy = (noise_power + alpha_p_1 * P_p(i));   

        parfor j=1:length(N_sen_th)
            %disp(strcat('N = ',num2str(N_sen_th(j))));     

            %% Determining the performance meteric R for all the different cases
            %% Ideal         
            %% Threshold   
            epsilon_id_th = 2 * Acc_energy/N_sen_th(j) *...
                gammaincinv(P_d_d, N_sen_th(j)/2, 'upper');

            %% Probability of false alarm and probability of dectection
            P_f_id_th(i,j) = gammainc( N_sen_th(j) * epsilon_id_th/(2 * noise_power), N_sen_th(j)/2, 'upper');
            P_d_id_th(i,j) =  P_d_d;

            R_id_th(i,j)  = (K - N_sen_th(j))/K * (P_H0 * (1 -  P_f_id_th(i,j)) * C_0 +...
                (1 - P_H0) * (1 - P_d_id_th(i,j)) * C_1);  


            if N_sen_th(j) > N_est_th
               N_sen_ni_th = N_sen_th(j)% - N_est_th;
               %% Average case
               if 1
                   %% Threshold  
                   test_points = linspace(epsilon_id_th, 0.1 * epsilon_id_th, num_test_points);
                   Exp_P_d_ac = zeros(1,num_test_points);
                   for k=1:num_test_points
                        func_exp_P_d = @(t) gammainc( N_sen_ni_th/2 * test_points(k) ./ t, N_sen_ni_th/2, 'upper') *...
                            N_est_th/Acc_energy .* exp(-N_est_th/2 * log(2) - gammaln(N_est_th/2)  +...
                            (N_est_th/2 - 1) * log(N_est_th * t/Acc_energy) +...
                            (-t * N_est_th/(2 * Acc_energy)));                    
                        Exp_P_d_ac(k) =  integral(func_exp_P_d, 0, test_points(k) * 100);
                        if Exp_P_d_ac(k) >= P_d_d
                           epsilon_ac_th = test_points(k);
                           %% Expected P_d 
                           P_d_ac_th(i,j) = Exp_P_d_ac(k);
                           break;
                        else
                           epsilon_ac_th = test_points(k); 
                        end
                   end    

                   %% Probability of false alarm and probability of dectection   
                   P_f_ac_th(i,j) = gammainc( N_sen_ni_th * epsilon_ac_th/...
                       (2 * noise_power), N_sen_ni_th/2, 'upper');  


                   %% Expected Rate
                   R_ac_th(i,j) = (K - N_sen_th(j))/K * (P_H0 * (1 -  P_f_ac_th(i,j)) * C_0 +...
                       (1 - P_H0) * (1 - P_d_ac_th(i,j)) * C_1); 
               end


                %% Outage case
                if 1
                    %% Threshold  
                    epsilon_oc_th = 4 * Acc_energy * gammaincinv(1 - mu, N_est_th/2, 'upper') * ...
                        gammaincinv(P_d_d, N_sen_ni_th/2, 'upper') / (N_est_th * N_sen_ni_th); 

                    %% Probability of false alarm and  expected probability of dectection
                    P_f_oc_th(i,j) = gammainc(N_sen_ni_th * epsilon_oc_th/...
                        (2 * noise_power),N_sen_ni_th/2, 'upper');   

                    func_exp_P_d = @(t) gammainc( N_sen_ni_th/2 * epsilon_oc_th ./ t, N_sen_ni_th/2, 'upper') *...
                        N_est_th/Acc_energy .* exp(-N_est_th/2 * log(2) - gammaln(N_est_th/2)  +...
                        (N_est_th/2 - 1) * log(N_est_th * t/Acc_energy) +...
                        (-t * N_est_th/(2 * Acc_energy)));
                    P_d_oc_th(i,j) =  integral(func_exp_P_d, 0, epsilon_oc_th * 100);

                    %% Expected Rate
                    R_oc_th(i,j) = (K - N_sen_th(j))/K * (P_H0 * (1 -  P_f_oc_th(i,j)) * C_0 +...
                    (1 - P_H0) * (1 - P_d_oc_th(i,j)) * C_1);
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
        
        disp(strcat('P_f_id_opt_th(i) = ',num2str(P_f_id_opt_th(i)))); 
        disp(strcat('P_f_ac_opt_th(i) = ',num2str(P_f_ac_opt_th(i)))); 
        disp(strcat('P_f_oc_opt_th(i) = ',num2str(P_f_oc_opt_th(i))));
        disp(strcat('P_d_id_opt_th(i) = ',num2str(P_d_id_opt_th(i)))); 
        disp(strcat('P_d_ac_opt_th(i) = ',num2str(P_d_ac_opt_th(i)))); 
        disp(strcat('P_d_oc_opt_th(i) = ',num2str(P_d_oc_opt_th(i)))); 
        disp(strcat('R_id_opt_th(i) = ',num2str(R_id_opt_th(i)))); 
        disp(strcat('R_ac_opt_th(i) = ',num2str(R_ac_opt_th(i)))); 
        disp(strcat('R_oc_opt_th(i) = ',num2str(R_oc_opt_th(i)))); 
        save('results_opt_thr_vs_SNR_AWGN_wo_nu_th2.mat');

    end
    quit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Throughput vs SNR @ optimum sensing time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_opt_thr_vs_SNR_AWGN_wo_nu_th2.mat');
figure(1);
Fontsize = 8;
plot(10 *log10(P_p), R_id_opt_th, 'c-', 'LineWidth',2);
hold on,
plot(10 *log10(P_p), R_ac_opt_th, 'b-', 'LineWidth', 1);
hold on,
plot(10 *log10(P_p), R_oc_opt_th, 'b--', 'LineWidth', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves cosmetic makeover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid on;
axis([min(10 *log10(P_p)) max(10 *log10(P_p))...
    min(R_oc_opt_th)  max(R_id_opt_th) * 1.015]);
ylabel('$\trs(\test = \SI{5}{ms},\ttsen)$ [bits/sec/Hz]','FontSize',Fontsize);
xlabel('$\snrrcvd$ [dB]','FontSize',Fontsize);

%hl = legend('$\trs$', '$\trsac$', '$\trsoc$');
%set(hl, 'position',[0.725 0.12 0.15 0.23]);
hl = legend('IM', 'EM-AC, Problem 1', 'EM-OC, Problem 2');
set(hl, 'position',[0.55 0.12 0.35 0.21]);

set(gca,'FontSize',Fontsize);
laprint(1, '../figures/fig_opt_thr_vs_SNR_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');

if 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Detection Probability vs SNR @ optimum sensing time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves theoretical analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2);
    Fontsize = 9;
    plot(10 *log10(P_p), P_d_id_opt_th, 'c', 'LineWidth', 4);
    hold on,
    plot(10 *log10(P_p), P_d_ac_opt_th, 'LineStyle', '-', 'LineWidth', 1);
    hold on,
    plot(10 *log10(P_p), P_d_oc_opt_th, 'LineStyle', '--', 'LineWidth', 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves cosmetic makeover
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    grid on;
    axis([min(10 *log10(P_p)) max(10 *log10(P_p))...
        min(P_d_oc_opt_th)  max(P_d_oc_opt_th) * 1.015]);
    ylabel('$\pd$','FontSize',Fontsize);
    xlabel('$\snrrcvd$ [dB]','FontSize',Fontsize);
    hl = legend('$\trs$', '$\trsac$', '$\trsoc$');
    set(hl, 'Location', 'SouthEast', 'FontSize', Fontsize);
    set(gca,'FontSize',Fontsize);
    laprint(2, '../figures/fig_opt_P_d_vs_SNR_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
end