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
P_p = 10.^([-10]/10);                     % Power transmitted by ST, the SNR received at SR 
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
num_test_points = 10000;                  % Test points for evaluating threshold 

%% Flags -- The flags below can be used to enable and thereby investigate individual cases 
Adjacent = 1;
Outage = 1;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th
   %% Theoretical parameters
    N_est_th = 1000:100:10000;             % -10 dB     % N_est_th = number of samples used for estimation = tau * f_s 

    N_sen_th = 50:50:25000;               % -10 dB     % N = includes the time in samples, for the ideal case their is no
                                                       % estimation time hence sensing time = N_sen_th and the throughput will 
                                                       % attain a nonzero value. For the accurate, the energy walls sensing 
                                                       % time begin after the estimation of the rcvd energy is performed.
                                                       % N_sen_ni_th = N_sen_th - N_est
    
             

    P_f_id_th = zeros(length(N_est_th), length(N_sen_th));   % Probability of false alarm at ST for the ideal case
    P_f_ac_th = zeros(length(N_est_th), length(N_sen_th));   % Probability of false alarm at ST due to accurate estimation of threshold
    P_f_oc_th = zeros(length(N_est_th), length(N_sen_th));   % Probability of false alarm at ST due to inaccurate estimation of threshold

    P_d_id_th = zeros(length(N_est_th), length(N_sen_th));   % Probability of detection at ST for the ideal case
    P_d_ac_th = zeros(length(N_est_th), length(N_sen_th));   % Probability of detection at ST due to accurate estimation of threshold
    P_d_oc_th = zeros(length(N_est_th), length(N_sen_th));   % Probability of detection at ST due to inaccurate estimation of threshold
    
    
    R_id_th = zeros(length(N_est_th), length(N_sen_th));               % Throughput at SR ideal case
    R_ac_th = zeros(length(N_est_th), length(N_sen_th));               % Throughput at SR with accurate estimation
    R_oc_th = zeros(length(N_est_th), length(N_sen_th));               % Throughput at SR with upper limit of the energy wall
    

    for i = 1:length(P_p) 
        
        disp(strcat('P_p(i) = ',num2str(P_p(i))));
        
        [C_0, C_1] = calc_capacities(alpha_p_2 * P_p(i) / noise_power,... 
            alpha_s * P_s / noise_power,....
            noise_power);                               % Capacities at SR, with and without intereference from PT  
         
        
        % Accurarte energy  
        Acc_energy = (noise_power + alpha_p_1 * P_p(i));    
        for j=1:length(N_est_th)
            disp(strcat('N_est = ',num2str(N_est_th(j))));
            N_est_th_tmp = N_est_th(j);
            parfor k=1:length(N_sen_th)
               epsilon_id_th = 0;
               epsilon_ac_th = 0;
               epsilon_oc_th = 0;
               %disp(strcat('N_sen_th = ',num2str(N_sen_th(k))));

               %% Determining the performance meteric R for all the different cases
               %% Ideal         
               %% Threshold  
                epsilon_id_th = 2 * Acc_energy/N_sen_th(k) *...
                gammaincinv(P_d_d, N_sen_th(k)/2, 'upper');               

               %% Probability of false alarm and probability of dectection
                P_f_id_temp(k) = gammainc( N_sen_th(k) * epsilon_id_th/(2 * noise_power),N_sen_th(k)/2, 'upper');        
                P_d_id_temp(k) =  P_d_d;


               %% Expected Rate
                R_id_th(k)  = (K - N_sen_th(k))/K * (P_H0 * (1 -  P_f_id_temp(k)) * C_0 +...
                    (1 - P_H0) * (1 - P_d_id_temp(k)) * C_1);  

                if N_sen_th(k) > N_est_th_tmp

                   N_sen_ni_th = N_sen_th(k);% - N_est_th(j);
                   %% Adjacent case
                   if Adjacent
                       %% Threshold  
                       test_points = linspace(epsilon_id_th, 0.1 * epsilon_id_th, num_test_points);
                       Exp_P_d_ac = zeros(1,num_test_points);
                       for l=1:num_test_points
                            func_exp_P_d = @(t) gammainc( N_sen_ni_th/2 * test_points(l) ./ t, N_sen_ni_th/2, 'upper') *...
                            N_est_th(j)/Acc_energy .* exp(-N_est_th(j)/2 * log(2) - gammaln(N_est_th(j)/2)  +...
                            (N_est_th(j)/2 - 1) * log(N_est_th(j) * t/Acc_energy) +...
                            (-t * N_est_th(j)/(2 * Acc_energy)));                    
                            Exp_P_d_ac(l) =  integral(func_exp_P_d, 0, test_points(l) * 100);
                            if Exp_P_d_ac(l) >= P_d_d
                                epsilon_ac_th = test_points(l);
                                %% Expected P_d 
                                P_d_ac_th(j,k) = Exp_P_d_ac(l);
                                break;
                            else
                                epsilon_ac_th = test_points(l); 
                            end
                       end    
                       
                       %% Probability of false alarm and probability of dectection
                       P_f_ac_th(j,k) = gammainc( N_sen_ni_th * epsilon_ac_th/...
                       (2 * noise_power), N_sen_ni_th/2, 'upper');                               
                       
                       %% Expected Rate
                       R_ac_th(j, k) = (K - N_sen_th(k))/K * (P_H0 * (1 -  P_f_ac_th(j,k)) * C_0 +...
                           (1 - P_H0) * (1 - P_d_ac_th(j,k)) * C_1); 
                   end

                    %% Outage case
                    if Outage
                        %% Threshold  
                        epsilon_oc_th = 4 * Acc_energy * gammaincinv(1 - mu, N_est_th(j)/2, 'upper') * ...
                            gammaincinv(P_d_d, N_sen_ni_th/2, 'upper') / (N_est_th(j) * N_sen_ni_th);            

                        %% Probability of false alarm and expected probability of dectection
                        P_f_oc_th(j, k) = gammainc( N_sen_ni_th * epsilon_oc_th/...
                            (2 * noise_power), N_sen_ni_th/2, 'upper');     

                        func_exp_P_d = @(t) gammainc( N_sen_ni_th/2 * epsilon_oc_th ./ t, N_sen_ni_th/2, 'upper') *...
                            N_est_th(j)/Acc_energy .* exp(-N_est_th(j)/2 * log(2) - gammaln(N_est_th(j)/2)  +...
                            (N_est_th(j)/2 - 1) * log(N_est_th(j) * t/Acc_energy) +...
                            (-t * N_est_th(j)/(2 * Acc_energy)));
                        P_d_oc_th(j, k) =  integral(func_exp_P_d, 0, epsilon_oc_th * 100);
                        
 
                        %% Expected Rate
                        R_oc_th(j, k) = (K - N_sen_th(k))/K * (P_H0 * (1 -  P_f_oc_th(j, k)) * C_0 +...
                        (1 - P_H0) * (1 - P_d_oc_th(j, k)) * C_1);
                    end
                end            
            end
            save('results_opt_thr_vs_est_time_sen_time_diff_SNR_wo_nu_m10_mu05_oc_th2.mat');
        end
    end
    quit;
end


if 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves theoretical analysis --  Optimum throughput vs est time    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('results_opt_thr_vs_est_time_sen_time_diff_SNR_wo_nu_m10_mu05_oc_th2.mat');

    %% ac, oc, mu = 0.05   
    Fontsize = 8;
    figure(3);
    scale1 = 1;
    scale2 = 4;
    surf(N_sen_th(1:scale2:length(N_sen_th)) * 1e-3, N_est_th(1:scale1:length(N_est_th)) * 1e-3,....
    R_oc_th(1:scale1:length(N_est_th), 1:scale2:length(N_sen_th)),...
    'EdgeColor','none','LineStyle','none');    
%     plot3( N_est_th(1:scale1:length(N_est_th)) * 1e-3, N_sen_th(1:scale2:length(N_sen_th)) * 1e-3,....
%     R_oc_th(1:scale1:length(N_est_th), 1:scale2:length(N_sen_th)));

    caxis([0,3.2]);
    zlabel('$\trs(\test,\tsen)$ [bits/sec/Hz]','FontSize',Fontsize);
    xlabel('$\tsen$ [ms]','FontSize',Fontsize);
    ylabel('$\test$ [ms]','FontSize',Fontsize);
    set(gca,'FontSize',Fontsize);

    [value index] = max(R_oc_th(:));
    [row colummn] = find(R_oc_th == value);

    N_oc_opt_est_th = N_est_th(row);
    N_oc_opt_sen_th = N_sen_th(colummn);
    
    hold on,
    h1 = plot3(N_oc_opt_sen_th * 1e-3, N_oc_opt_est_th * 1e-3, value, 'ws');
 

   %legend([h1], '$\trs(\ttest,\ttsen)$');
    laprint(3, '../figures/fig_opt_thr_vs_est_time_sen_time_oc_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
    
    figure(4);
    surf(N_sen_th(1:scale2:length(N_sen_th)) * 1e-3, N_est_th(1:scale1:length(N_est_th)) * 1e-3,....
    R_ac_th(1:scale1:length(N_est_th), 1:scale2:length(N_sen_th)),...
    'EdgeColor','none','LineStyle','none');
    caxis([0,3.2]);
    
    [value index] = max(R_ac_th(:));
    [row colummn] = find(R_ac_th == value);

    N_ac_opt_est_th = N_est_th(row);
    N_ac_opt_sen_th = N_sen_th(colummn); 
    
        
    hold on,
    h1 = plot3(N_ac_opt_sen_th * 1e-3, N_ac_opt_est_th * 1e-3, value, 'ws');
    
    zlabel('$\trs(\test,\tsen)$ [bits/sec/Hz]','FontSize',Fontsize);
    xlabel('$\tsen$ [ms]','FontSize',Fontsize);
    ylabel('$\test$ [ms]','FontSize',Fontsize);
    set(gca,'FontSize',Fontsize);

    laprint(4, '../figures/fig_opt_thr_vs_est_time_sen_time_ac_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
    
    
end
