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
mu = 0.10;                                % Outage Probability on probability of detection

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th
   %% Theoretical parameters
    N_est_th = 500:100:10000;            % -10 dB     % N_est_th = number of samples used for estimation = tau * f_s 
    %N_est_th = 200:100:5000;              % -5 dB      % N_est_th = number of samples used for estimation = tau * f_s 

    N_th = 100:100:25000;          % -10 dB     % N = includes the time in samples, for the ideal case their is no

    %N_th = 20:20:7000;       % -5 dB      % N = includes the time in samples, for the ideal case their is no
                                                       % estimation time hence sensing time = N_th and the throughput will 
                                                       % attain a nonzero value. For the accurate, the energy walls sensing 
                                                       % time begin after the estimation of the rcvd energy is performed.
                                                       % N_sen = N_th - N_est
    
    % For storing intermediate values of probability of false alarm
    P_f_id_temp = zeros(1, length(N_th));                 
    P_f_ac_temp = zeros(1, length(N_th));                
    P_f_oc_temp = zeros(1, length(N_th));                

    % For storing intermediate values of probability of detection
    P_d_id_temp = zeros(1, length(N_th));                 
    P_d_ac_temp = zeros(1, length(N_th));                 
    P_d_oc_temp = zeros(1, length(N_th));                 

    P_f_id_th = zeros(length(P_p), length(N_est_th));   % Probability of false alarm at ST for the ideal case
    P_f_ac_th = zeros(length(P_p), length(N_est_th));   % Probability of false alarm at ST due to accurate estimation of threshold
    P_f_oc_th = zeros(length(P_p), length(N_est_th));   % Probability of false alarm at ST due to inaccurate estimation of threshold

    P_d_id_th = zeros(length(P_p), length(N_est_th));   % Probability of detection at ST for the ideal case
    P_d_ac_th = zeros(length(P_p), length(N_est_th));   % Probability of detection at ST due to accurate estimation of threshold
    P_d_oc_th = zeros(length(P_p), length(N_est_th));   % Probability of detection at ST due to inaccurate estimation of threshold
    
    
    R_id_th = zeros(1, length(N_th));                   % Throughput at SR ideal case
    R_ac_th = zeros(1, length(N_th));                   % Throughput at SR with accurate estimation
    R_oc_th = zeros(1, length(N_th));                   % Throughput at SR with upper limit of the energy wall
    
    R_id_opt_th = zeros(length(P_p), length(N_est_th)); % Throughput at SR ideal case
    R_ac_opt_th = zeros(length(P_p), length(N_est_th)); % Throughput at SR with accurate estimation
    R_oc_opt_th = zeros(length(P_p), length(N_est_th)); % Throughput at SR with upper limit of the energy wall
    
    for i = 1:length(P_p) 
        
        disp(strcat('P_p(i) = ',num2str(P_p(i))));
        
        % Accurarte energy  
        Acc_energy = (noise_power + alpha_p_1 * P_p(i));    
        for j=1:length(N_est_th)
            disp(strcat('N_est = ',num2str(N_est_th(j))));
            parfor k=1:length(N_th)
               warning off;
               epsilon_id_th = 0;
               epsilon_ac_th = 0;
               epsilon_oc_th = 0;
               %disp(strcat('N_sen = ',num2str(N_th(k))));


               %disp(strcat('N = ',num2str(N_th(k))));     

               %% Determining the performance meteric R for all the different cases
               %% Ideal         
               %% Threshold  
                epsilon_id_th = qfuncinv(P_d_d) * sqrt(2/N_th(k)) * Acc_energy +...
                    Acc_energy;                  

               %% Probability of false alarm and probability of dectection
                P_f_id_temp(k) = qfunc((epsilon_id_th - noise_power)/...
                    (sqrt(2/N_th(k)) * noise_power));        

               %% Expected Rate
                C_0_th = log2(1 + alpha_s * P_s / noise_power); 
                C_1_th = log2(1 + alpha_s * P_s / (P_p(i) * alpha_p_2 +...
                    noise_power) );
                P_d_id_temp(k) =  P_d_d;

                R_id_th(k)  = (K - N_th(k))/K * (P_H0 * (1 -  P_f_id_temp(k)) * C_0_th +...
                    (1 - P_H0) * (1 - P_d_id_temp(k)) * C_1_th);  

                if N_th(k) > N_est_th(j)
                   N_sen_th = N_th(k) - N_est_th(j);
                   %% Adjacent case
                   if 0
                       %% Threshold  
                       num_test_points = 1000;
                       test_points = linspace(Acc_energy* (1 + sqrt(2/N_sen_th)* qfuncinv(0.35)), ...
                           Acc_energy* (1 + sqrt(2/N_sen_th)* sqrt(2) * (-erfcinv(10^-6 * eps))),num_test_points);
                       Exp_P_d_ac = zeros(1,num_test_points);
                       for l=1:num_test_points
                            ts = sqrt(2/N_sen_th);
                            te = sqrt(2/N_est_th(j));
                            func_exp_P_d = @(t)  t .*  test_points(l) .*...
                                exp( (qfuncinv(t)/sqrt(2)).^2 - (test_points(l)./(1 + ts * qfuncinv(t)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
                                * ts ./ (Acc_energy * te * ( 1 + ts * qfuncinv(t)).^2);
                            Exp_P_d_ac(l) =  integral(func_exp_P_d,0,1);
                            if Exp_P_d_ac(l) > 0.8 && Exp_P_d_ac(l) < 0.92 
                                my_int =  my_integral_pd(Acc_energy, test_points(l),...
                                    ts, te, eps*10^-6, eps, 10^5);
                                Exp_P_d_ac(l) =  integral(func_exp_P_d, 0 , 1  - eps) + my_int;
                            end
                            if Exp_P_d_ac(l) >= P_d_d
                               epsilon_ac_th = test_points(l); 
                               break;
                            else
                               epsilon_ac_th = test_points(l); 
                            end
                       end    
                       %% Probability of false alarm and probability of dectection
                       P_f_ac_temp(k) = qfunc((epsilon_ac_th - noise_power)/...
                           (sqrt(2/N_sen_th) * noise_power));                                
                       % Test if the integral converges by determining the
                       % integral over the density to check if the value is
                       % equal to 1
                       if 1 
                            ts = sqrt(2/N_sen_th);
                            te = sqrt(2/N_est_th(j));
                            
                            P_d_ac_temp(k) =  integral(func_exp_P_d, 0, 1);                            
                            if P_d_ac_temp(k) < 0.89
                                my_int =  my_integral_pd(Acc_energy, epsilon_ac_th,...
                                    ts, te, eps*10^-8, eps, 10^7);
                                func_exp_P_d = @(t)  t .*  epsilon_ac_th .*...
                                exp( (qfuncinv(t)/sqrt(2)).^2 - (epsilon_ac_th./(1 + ts * qfuncinv(t)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
                                * ts ./ (Acc_energy * te * ( 1 + ts * qfuncinv(t)).^2);
                                P_d_ac_temp(k) =  integral(func_exp_P_d,0,1 - eps) + my_int;
                                if P_d_ac_temp(k) < 0.88 || P_d_ac_temp(k) > 0.92
                                    disp(strcat('AC: The integral failed to converge ','N_th: ',num2str(N_th(k))));
                                    R_ac_th(k) = 0;                                
                                    continue;
                                end
                            end
                       end
                       %% Expected Rate
                       C_0_th = log2(1 + alpha_s * P_s / noise_power); 
                       C_1_th = log2(1 + alpha_s * P_s / (P_p(i) * alpha_p_2 +...
                           noise_power) );

                       R_ac_th(k) = (K - N_th(k))/K * (P_H0 * (1 -  P_f_ac_temp(k)) * C_0_th +...
                           (1 - P_H0) * (1 - P_d_ac_temp(k)) * C_1_th); 
                   end


                    %% Outage case
                    if 1
                        %% Threshold  
                        epsilon_oc_th = Acc_energy * (1 + qfuncinv(1 - mu) * (sqrt(2/N_est_th(j))))...
                             * (1 + qfuncinv(P_d_d) * (sqrt(2/N_sen_th)));            

                        %% Probability of false alarm and expected probability of dectection
                        P_f_oc_temp(k) = qfunc((epsilon_oc_th - noise_power)/...
                            (sqrt(2/N_sen_th) * noise_power));   

                        ts = sqrt(2/N_sen_th);
                        te = sqrt(2/N_est_th(j));
                        func_exp_P_d = @(t)  t .*  epsilon_oc_th .*...
                            exp( (qfuncinv(t)/sqrt(2)).^2 - (epsilon_oc_th./(1 + ts * qfuncinv(t)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
                            * ts ./ (Acc_energy * te * ( 1 + ts * qfuncinv(t)).^2);
                        P_d_oc_temp(k) =  integral(func_exp_P_d,0,1); 
                        
                        % Test if the integral converges by determining the
                        % integral over the density to check if the value is
                        % equal to 1
                        if 1 
                             ts = sqrt(2/N_sen_th);
                             te = sqrt(2/N_est_th(j));
                             if P_d_oc_temp(k) < 0.99
                                 my_int =  my_integral_pd(Acc_energy, epsilon_oc_th,...
                                     ts, te, eps*10^-8, eps, 10^7);
                                 P_d_oc_temp(k) =  integral(func_exp_P_d, 0, 1) + my_int;
                                 if P_d_oc_temp(k) < 0.90
                                    disp(strcat('OC: The integral failed to converge ','N_th: ',num2str(N_th(k))));
                                    R_oc_th(k) = 0;
                                    continue;                                    
                                 end
                             end
                        end

                        %% Expected Rate
                        C_0_th = log2(1 + alpha_s * P_s / noise_power); 
                        C_1_th = log2(1 + alpha_s * P_s / (P_p(i) * alpha_p_2 +...
                            noise_power) );

                        R_oc_th(k) = (K - N_th(k))/K * (P_H0 * (1 -  P_f_oc_temp(k)) * C_0_th +...
                        (1 - P_H0) * (1 - P_d_oc_temp(k)) * C_1_th);
                    end
                end            
            end

            % The probability of detecttion and probability of false alarm at optimum sensing time
            % ideal model
            [temp_id index_id] = max(R_id_th);
            R_id_opt_th(i,j) = temp_id;
            P_f_id_th(i,j) = P_f_id_temp(index_id);
            P_d_id_th(i,j) = P_d_id_temp(index_id);

            % average constraint case
            [temp_ac index_ac] = max(R_ac_th);
            R_ac_opt_th(i,j) = temp_ac;
            P_f_ac_th(i,j) = P_f_ac_temp(index_ac);
            P_d_ac_th(i,j) = P_d_ac_temp(index_ac);

            % outage constraint case
            [temp_oc index_oc] = max(R_oc_th);
            R_oc_opt_th(i,j) = temp_oc;
            P_f_oc_th(i,j) = P_f_oc_temp(index_oc);
            P_d_oc_th(i,j) = P_d_oc_temp(index_oc);
            
            disp(strcat('P_f_id_th(i,j) = ',num2str(P_f_id_th(i,j)))); 
            disp(strcat('P_d_id_th(i,j) = ',num2str(P_d_id_th(i,j)))); 
            disp(strcat('R_id_opt_th(i,j) = ',num2str(R_id_opt_th(i,j)))); 
            
            disp(strcat('AC: N_th(index_ac) = ',num2str(N_th(index_ac)))); 
            disp(strcat('P_f_ac_th(i,j) = ',num2str(P_f_ac_th(i,j)))); 
            disp(strcat('P_d_ac_th(i,j) = ',num2str(P_d_ac_th(i,j)))); 
            disp(strcat('R_ac_opt_th(i,j) = ',num2str(R_ac_opt_th(i,j)))); 
            
            disp(strcat('OC: N_th(index_oc) = ',num2str(N_th(index_oc)))); 
            disp(strcat('P_f_oc_th(i,j) = ',num2str(P_f_oc_th(i,j)))); 
            disp(strcat('P_d_oc_th(i,j) = ',num2str(P_d_oc_th(i,j)))); 
            disp(strcat('R_oc_opt_th(i,j) = ',num2str(R_oc_opt_th(i,j)))); 

        end
    end
    save('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu10_oc_th.mat');
    quit;
end


if 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves theoretical analysis --  Proability of detection vs est
    %   time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fontsize = 9;
    plot(N_est_th * 1e-3, P_d_id_th(1,:), 'LineWidth', 1.5);
    hold on, 
    plot(N_est_th * 1e-3, P_d_ac_th(1,:), 'LineWidth', 1.5);
    hold on,
    plot(N_est_th * 1e-3, P_d_oc_th(1,:), 'LineWidth', 1.5);
    %hold on,
    %plot(N_est_th * 1e-3, P_d_oc_th(2,:), 'LineWidth', 1.5, 'HandleVisibility','off');
    %hold on,
    %plot(N_est_th * 1e-3, P_d_oc_th(3,:), 'LineWidth', 1.5, 'HandleVisibility','off');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves cosmetic makeover
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    grid on;
    axis([min(N_est_th) * 1e-3 max(N_est_th) * 1e-3 min(P_d_id_th(1,:)) * 0.995  1.02]);
    ylabel('$\pd$','FontSize',Fontsize+2);
    xlabel('$\test$ = [ms]','FontSize',Fontsize+2);
    hl = legend('$\pd$', '$\pdac$', '$\pdoc$');
    set(hl, 'Location', 'SouthEast', 'FontSize', Fontsize);
    set(gca,'FontSize',Fontsize);
    %laprint(1, '../figures/fig_P_d_vs_est_time_diff_snr_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
    %    'on', 'factor',0.5, 'keepfontprops', 'on');
end

if 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves theoretical analysis --  Proability of false alarm vs est
    %   time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fontsize = 9;
    figure(2)
    semilogy(N_est_th * 1e-3, P_f_id_th(1,:), 'LineWidth', 1.5);
    hold on,
    semilogy(N_est_th * 1e-3, P_f_ac_th(1,:), 'LineWidth', 1.5);
    hold on,
    semilogy(N_est_th * 1e-3, P_f_oc_th(1,:), 'LineWidth', 1.5);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves cosmetic makeover
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    grid on;
    axis tight;
    ylabel('$\pd$','FontSize',Fontsize+2);
    xlabel('$\test$ = [ms]','FontSize',Fontsize+2);
    hl = legend('$\pf$', '$\pfac$', '$\pfoc$');
    set(hl, 'Location', 'SouthEast', 'FontSize', Fontsize);
    set(gca,'FontSize',Fontsize);
    %laprint(1, '../figures/fig_P_f_vs_est_time_diff_snr_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
    %    'on', 'factor',0.5, 'keepfontprops', 'on');
end

if 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves theoretical analysis --  Optimum throughput vs est time    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu05_th.mat');
    Fontsize = 9;
    figure(3);
    scale = 6:length(N_est_th);
    plot(N_est_th(scale) * 1e-3, R_id_opt_th(1,(scale)), 'LineWidth', 1.5);
    hold on,
    plot(N_est_th(scale) * 1e-3, R_ac_opt_th(1,(scale)), 'LineWidth', 1.5);
    hold on,
    plot(N_est_th(scale) * 1e-3, R_oc_opt_th(1,(scale)), 'LineWidth', 1.5);
    
    
    [temp index] = max(R_ac_opt_th(scale));
    R_ac_opt_sen_opt_est = R_ac_opt_th(1,index);
    N_ac_opt_est = N_est_th(index);

    [temp index] = max(R_oc_opt_th(scale));
    R_oc_opt_sen_opt_est = R_oc_opt_th(1,index);
    N_oc_opt_est = N_est_th(index);
    
    
    hold on,
    plot(N_ac_opt_est * 1e-3, R_ac_opt_sen_opt_est, 'ks', 'LineWidth',1, 'HandleVisibility','off');
    hold on,
    plot(N_oc_opt_est * 1e-3, R_oc_opt_sen_opt_est, 'ks', 'LineWidth',1);
    
    grid on;
    axis([min(N_est_th(scale) * 1e-3) max(N_est_th(scale) * 1e-3) min(R_oc_opt_th(1,(scale))) max(R_id_opt_th(1,(scale))) * 1.01]);
    
    
    %% mu = 0.10 
    if 1
        load('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu10_oc_th.mat');

        hold on,
        plot(N_est_th(scale) * 1e-3, R_oc_opt_th(1,(scale)), 'LineWidth', 1.5, 'HandleVisibility','off');

        [temp index] = max(R_oc_opt_th(scale));
        R_oc_opt_sen_opt_est = R_oc_opt_th(1,index);
        N_oc_opt_est = N_est_th(index);
        hold on,
        plot(N_oc_opt_est * 1e-3, R_oc_opt_sen_opt_est, 'ks', 'LineWidth',1, 'HandleVisibility','off');
    end
    %% mu = 0.15 
    if 1
        load('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu15_oc_th.mat');

        hold on,
        plot(N_est_th(scale) * 1e-3, R_oc_opt_th(1,(scale)), 'LineWidth', 1.5);

        [temp index] = max(R_oc_opt_th(scale));
        R_oc_opt_sen_opt_est = R_oc_opt_th(1,index);
        N_oc_opt_est = N_est_th(index);
        hold on,
        plot(N_oc_opt_est * 1e-3, R_oc_opt_sen_opt_est, 'ks', 'LineWidth',1);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves cosmetic makeover
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ylabel('$\trs$ [bits/sec/Hz]','FontSize',Fontsize);
    xlabel('$\test$ [ms]','FontSize',Fontsize);
    hl = legend('$\trs$', '$\trsac$', '$\trsoc$','Opt. $\ttest$');
    set(hl, 'Location', 'SouthEast', 'FontSize', Fontsize);
    set(gca,'FontSize',Fontsize);
    laprint(3, '../figures/fig_opt_thr_vs_est_time_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
end