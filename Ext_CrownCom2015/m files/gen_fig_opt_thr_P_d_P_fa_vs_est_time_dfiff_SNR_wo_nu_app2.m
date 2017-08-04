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

warning off;

sim = 00;                                                     % Enable( = 1) to perform simulation, theoretical analysis
th = 00;                                                      % disable ( = 0) to plot curves, data is read from a file
                                                              % Theoretical anaylsis is also included of simulation as
                                                              % numerical integration is involved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   System Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_s = 10^(-10/10);                                            % Power transmitted by PR, the SNR received at ST 
P_p = 10.^([-10]/10);                                         % Power transmitted by ST, the SNR received at SR 
                                                   % varried using this parameter
noise_power = 10^(-100/10);                                   % noise power -100 dBm
f_s = 1e6;                                                    % 1 MHz one band
K = 0.1 * f_s;                                                % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha_p_1 = 10^(-100/10);                                     % True Path loss between ST and PR   
alpha_p_2 = 10^(-100/10);                                     % True Path loss between PR and SR   
alpha_s = 10^(-080/10);                                       % True Path loss between ST and SR 
P_H0 = 0.8;                                                   % Probability of Hypothesis 0
P_d_d = 0.90;                                                 % Constraint on Probability of detection P_d
mu = 0.05;                                                    % Outage Probability on probability of detection
num_test_points = 10000;                                      % Test points for evaluating threshold 

N_s = 10;                                                     % Number of Pilot Symobols 
sigma_est_error = noise_power/N_s;                            % Variance of the Estimation error for h_s 
dof_s = 1;                                                    % Degree of freedom = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 5e4;                                                      % Number of realizations            
N_est_sim = 1000:500:10000;                                   % Estimation time 
g_s_hat = zeros(1,M);                                         % Estimated power gain for the channel g_s        
P_f_ac_sim = zeros(length(P_p), length(N_est_sim));           % Probability of false alarm at ST due to inaccurate estimation of threshold
P_f_oc_sim = zeros(length(P_p), length(N_est_sim));           % Probability of false alarm at ST due to inaccurate estimation of threshold
P_d_oc_sim = zeros(length(P_p), length(N_est_sim));           % Probability of detection at ST due to inaccurate estimation of threshold
P_d_oc_sim = zeros(length(P_p), length(N_est_sim));           % Probability of detection at ST due to inaccurate estimation of threshold
R_ac_opt_sim = zeros(length(P_p), length(N_est_sim));         % Throughput at SR with upper limit of the energy wall
R_oc_opt_sim = zeros(length(P_p), length(N_est_sim));         % Throughput at SR with upper limit of the energy wall
N_sen_ac_opt_sim = zeros(1, M);                               % Different realizations of the optimum sensing correponding to different channel realizations
N_sen_oc_opt_sim = zeros(1, M);                               % Different realizations of the optimum sensing correponding to different channel realizations
Exp_N_sen_ac_opt_sim = zeros(length(P_p), length(N_est_sim)); % Expected Optimum sensing time for different values of estimation time
Exp_N_sen_oc_opt_sim = zeros(length(P_p), length(N_est_sim)); % Expected Optimum sensing time for different values of estimation time



%% Flags -- The flags below can be used to enable and thereby investigate individual cases 
Adjacent = 1;
Outage = 1;
    

if sim


    for i = 1:length(P_p)         
        disp(strcat('P_p(i) = ',num2str(P_p(i))));        
        % Accurarte energy  
        True_rcvd_energy_hp1 = (noise_power + alpha_p_1 * P_p(i));    
        for j=1:length(N_est_sim)
            disp(strcat('N_est = ',num2str(N_est_sim(j))));
            N_est_th_local = N_est_sim(j);
            
             %% Estimation of the sensing channel 
             parfor l=1:M
                 P_rcvd_est_sen_sim(l) = mean(random('norm',...
                 random('norm',0, sqrt(P_p * alpha_p_1), 1, N_est_th_local),...
                 sqrt(noise_power), 1, N_est_th_local).^2);                               
             end  

             %% Estimation of the access channel
             g_s = P_s * alpha_s;
             parfor l=1:M
                 g_s_hat(l) = sum(normrnd(sqrt(g_s), sqrt(sigma_est_error), 1, dof_s).^2);
             end

             %% Estimation of the interfererence channel
             parfor l=1:M
                 P_rcvd_est_int_sim(l) = mean(random('norm',...
                 random('norm',0, sqrt(P_p * alpha_p_2), 1, N_est_th_local),...
                 sqrt(noise_power), 1, N_est_th_local).^2);                               
             end  


             %% Determine the variations in the performance parameters C_0, C_1, beacuse the performance parameters 
             %% P_d, P_fa depends on sensing time hence they are calculated inside the loop
             C_0 =  log2(1 + g_s_hat/noise_power);
             C_1 =  log2(1 + g_s_hat./P_rcvd_est_int_sim);  
             
             %% Sensing time extends from the estimation time
             N_sen_sim = N_est_th_local:10:30000;              
             
             
             %% AC Threshold
             % Ideal case is required for the
             epsilon_ac_sim = zeros(1, length(N_sen_sim));
             epsilon_id_th = 2 * True_rcvd_energy_hp1./N_sen_sim .*...
                gammaincinv(P_d_d, N_sen_sim/2, 'upper'); 
             parfor m = 1:length(N_sen_sim)                 
                 Exp_P_d_ac = 0;
                 disp(strcat('AC iteration = ',num2str(m)));  
                 test_points = linspace(epsilon_id_th(m), 0.1 * epsilon_id_th(m), num_test_points);
                 for l = 1:num_test_points
                     func_exp_P_d = @(t) gammainc( N_sen_sim(m)/2 * test_points(l) ./ t, N_sen_sim(m)/2, 'upper') *...
                     N_est_sim(j)/True_rcvd_energy_hp1 .* exp(-N_est_sim(j)/2 * log(2) - gammaln(N_est_sim(j)/2)  +...
                     (N_est_sim(j)/2 - 1) * log(N_est_sim(j) * t/True_rcvd_energy_hp1) +...
                     (-t * N_est_sim(j)/(2 * True_rcvd_energy_hp1)));                    
                     Exp_P_d_ac =  integral(func_exp_P_d, 0, test_points(l) * 100);
                     if Exp_P_d_ac >= P_d_d
                         epsilon_ac_sim(m) = test_points(l);
                            break;
                     else
                         epsilon_ac_sim(m) = test_points(l); 
                     end
                 end
             end
             
             %% OC Threshold     
             epsilon_oc_sim = 4 * True_rcvd_energy_hp1 * gammaincinv(1 - mu, N_est_th_local/2, 'upper') * ...
                   gammaincinv(P_d_d, N_sen_sim/2, 'upper') ./ (N_est_th_local * N_sen_sim); 

             
             
             parfor k=1:M
                 
                   %disp(strcat('iteration = ',num2str(k)));    
                   if Adjacent
                       %% Average Constraint                      

                       %% Probability of false alarm and probability of detection
                       P_f_ac_local = gammainc(N_sen_sim/2 .* epsilon_ac_sim/noise_power, N_sen_sim/2, 'upper');
                       P_d_ac_local = gammainc(N_sen_sim/2 .* epsilon_ac_sim/P_rcvd_est_sen_sim(k), N_sen_sim/2, 'upper');

                       [values index] = max((K - N_sen_sim)/K .* (P_H0 * (1 -  P_f_ac_local) *...
                           C_0(k) +  (1 - P_H0) * (1 - P_d_ac_local) * C_1(k)));   

                        N_sen_ac_opt_sim(k) = N_sen_sim(index);
                   end
                   
                   if Outage
                       %% Outage constraint            
                       
                       %% Probability of false alarm and probability of detection
                       P_f_oc_local = gammainc(N_sen_sim/2 .* epsilon_oc_sim/noise_power ,N_sen_sim/2, 'upper');
                       P_d_oc_local = gammainc(N_sen_sim/2 .* epsilon_oc_sim/P_rcvd_est_sen_sim(k) ,N_sen_sim/2, 'upper');

                       [values index] = max((K - N_sen_sim)/K .* (P_H0 * (1 -  P_f_oc_local) *...
                           C_0(k) +  (1 - P_H0) * (1 - P_d_oc_local) * C_1(k)));   

                        N_sen_oc_opt_sim(k) = N_sen_sim(index);
                        %disp(strcat('N_sen_oc_opt_sim(k) = ',num2str(N_sen_oc_opt_sim(k))));
                   end

            end

            % The probability of detection and probability of false alarm at optimum sensing time
            
            %% Average constraint case   
            Exp_N_sen_ac_opt_sim(i,j) =  round(mean(N_sen_ac_opt_sim)/10) * 10; 
            
            %% Determine new threshold at the optimum expected sensing time 
            epsilon_ac_opt_sim = 0;
            epsilon_id_th = 0;
            epsilon_id_th = 2 * True_rcvd_energy_hp1./Exp_N_sen_ac_opt_sim(i,j) .*...
                gammaincinv(P_d_d, Exp_N_sen_ac_opt_sim(i,j)/2, 'upper');
            test_points = linspace(epsilon_id_th, 0.1 * epsilon_id_th, num_test_points);
            Exp_P_d_ac = 0;
            for l = 1:num_test_points
                func_exp_P_d = @(t) gammainc( Exp_N_sen_ac_opt_sim(i,j)/2 * test_points(l) ./ t, Exp_N_sen_ac_opt_sim(i,j)/2, 'upper') *...
                N_est_sim(j)/True_rcvd_energy_hp1 .* exp(-N_est_sim(j)/2 * log(2) - gammaln(N_est_sim(j)/2)  +...
                (N_est_sim(j)/2 - 1) * log(N_est_sim(j) * t/True_rcvd_energy_hp1) +...
                (-t * N_est_sim(j)/(2 * True_rcvd_energy_hp1)));                    
                Exp_P_d_ac =  integral(func_exp_P_d, 0, test_points(l) * 100);
                if Exp_P_d_ac >= P_d_d
                    epsilon_ac_opt_sim = test_points(l);
                       break;
                else
                    epsilon_ac_opt_sim = test_points(l); 
                end
            end
            
           
            P_f_ac_sim(i,j) = gammainc( Exp_N_sen_ac_opt_sim(i,j)/2 .* epsilon_ac_opt_sim/noise_power , Exp_N_sen_ac_opt_sim(i,j)/2, 'upper');
            P_d_ac_sim(i,j) = mean(gammainc(Exp_N_sen_ac_opt_sim(i,j)/2 .* epsilon_ac_opt_sim./P_rcvd_est_sen_sim, Exp_N_sen_ac_opt_sim(i,j)/2, 'upper'));

            R_ac_opt_sim(i,j) = (K - Exp_N_sen_ac_opt_sim(i,j))/K .* (P_H0 * (1 -  P_f_ac_sim(i,j)) *...
                   mean(C_0) +  (1 - P_H0) * (1 - P_d_ac_sim(i,j)) * mean(C_1));
            
               
            %% Outage constraint case            
            Exp_N_sen_oc_opt_sim(i,j) =  round(mean(N_sen_oc_opt_sim)/10) * 10; 
            
            %% Determine new threshold at the optimum expected sensing time 
            epsilon_oc_opt_sim = 4 * True_rcvd_energy_hp1 * gammaincinv(1 - mu, N_est_th_local/2, 'upper') * ...
                   gammaincinv(P_d_d, Exp_N_sen_oc_opt_sim(i,j)/2, 'upper') ./ (N_est_th_local * Exp_N_sen_oc_opt_sim(i,j)); 

           
            P_f_oc_sim(i,j) = gammainc( Exp_N_sen_oc_opt_sim(i,j)/2 .* epsilon_oc_opt_sim/noise_power , Exp_N_sen_oc_opt_sim(i,j)/2, 'upper');
            P_d_oc_sim(i,j) = mean(gammainc(Exp_N_sen_oc_opt_sim(i,j)/2 .* epsilon_oc_opt_sim./P_rcvd_est_sen_sim, Exp_N_sen_oc_opt_sim(i,j)/2, 'upper'));

            R_oc_opt_sim(i,j) = (K - Exp_N_sen_oc_opt_sim(i,j))/K .* (P_H0 * (1 -  P_f_oc_sim(i,j)) *...
                   mean(C_0) +  (1 - P_H0) * (1 - P_d_oc_sim(i,j)) * mean(C_1));
            
            disp(strcat('AC: Exp_N_sen_ac_opt_sim(index_oc) = ',num2str(Exp_N_sen_ac_opt_sim(i,j)))); 
            disp(strcat('P_f_ac_sim(i,j) = ',num2str(P_f_ac_sim(i,j)))); 
            disp(strcat('P_d_ac_sim(i,j) = ',num2str(P_d_ac_sim(i,j)))); 
            disp(strcat('R_ac_opt_sim(i,j) = ',num2str(R_ac_opt_sim(i,j)))); 
            
            disp(strcat('OC: Exp_N_sen_oc_opt_sim(index_oc) = ',num2str(Exp_N_sen_oc_opt_sim(i,j)))); 
            disp(strcat('P_f_oc_sim(i,j) = ',num2str(P_f_oc_sim(i,j)))); 
            disp(strcat('P_d_oc_sim(i,j) = ',num2str(P_d_oc_sim(i,j)))); 
            disp(strcat('R_oc_opt_sim(i,j) = ',num2str(R_oc_opt_sim(i,j)))); 
            save('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu05_oc_sim2_app2.mat');
        end
    end
    quit;
end

if 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves theoretical analysis --  Proability of detection vs est
    %   time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu05_oc_sim2_app2.mat');
    
    %% ac, oc, mu = 0.05    
    Fontsize = 9;
    diff = 5;
    pts = [10:diff:length(P_d_id_th(1,:)), length(P_d_id_th(1,:))];
    figure(1);
    plot(N_est_th(pts) * 1e-3, P_d_id_th(1,pts),'c-', 'LineWidth', 4);
    hold on, 
    plot(N_est_th(pts) * 1e-3, P_d_ac_th(1,pts),'-', 'LineWidth', 2);
    hold on,
    plot(N_est_th(pts) * 1e-3, P_d_oc_th(1,pts),'--', 'LineWidth', 2);
    
    %% oc, mu = 0.10 
    if 1
        load('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu10_oc_th2.mat');

        hold on,
        plot(N_est_th(pts) * 1e-3, P_d_oc_th(1,pts),'--', 'LineWidth', 2);
    end
    
    %% oc, mu = 0.15 
    if 1
        load('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu15_oc_th2.mat');

        hold on,
        plot(N_est_th(pts) * 1e-3, P_d_oc_th(1,pts),'--', 'LineWidth', 2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves cosmetic makeover
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    grid on;
    axis([min(N_est_th((pts))) * 1e-3 max(N_est_th((pts))) * 1e-3 min(P_d_id_th(1,(pts))) * 0.93  1.00]);
    ylabel('$\e{\pd}{\pd}$','FontSize',Fontsize);
    xlabel('$\test$ [ms]','FontSize',Fontsize);
    %hl = legend('$\pd$', '$\pdac$', '$\pdoc$');
    %set(hl, 'position',[0.758 0.118 0.09 0.23]);
    hl = legend('IM', 'EM-AC', 'EM-OC');
    set(hl, 'position',[0.615 0.12 0.28 0.23]);
    set(gca,'FontSize',Fontsize);
    laprint(1, '../figures/fig_P_d_vs_est_time_diff_mu_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
end

if 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves theoretical analysis --  Proability of false alarm vs est
    %   time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% ac, oc, mu = 0.05    
    load('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu05_oc_sim2_app2.mat');    
    Fontsize = 9;
    figure(2);
    diff1 = 5;
    pts1 = [10:diff1:length(P_f_id_th(1,:)), length(P_f_id_th(1,:))]
    diff2 = 10;
    pts2 = [7:diff1:length(P_f_id_th(1,:)), length(P_f_id_th(1,:))];
    semilogy(N_est_th(pts1) * 1e-3, P_f_id_th(1,pts1),'c-', 'LineWidth', 4);
    hold on,
    semilogy(N_est_th(pts2) * 1e-3, P_f_ac_th(1,pts2),'-', 'LineWidth', 2);
    hold on,
    semilogy(N_est_th(pts1) * 1e-3, P_f_oc_th(1,pts1),'--', 'LineWidth', 2);    
    axis([min(N_est_th(pts1)) * 1e-3 max(N_est_th(pts1))* 1e-3...
        1e-4 * .5 1e-1 * 2]);
    
    %% oc, mu = 0.10 
    if 1
        load('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu10_oc_th2.mat');

        hold on,
        semilogy(N_est_th(pts1) * 1e-3, P_f_oc_th(1,pts1),'--', 'LineWidth', 2);           
    end
    
    %% oc, mu = 0.15 
    if 1
        load('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu15_oc_th2.mat');

        hold on,
        semilogy(N_est_th(pts1) * 1e-3, P_f_oc_th(1,pts1),'--', 'LineWidth', 2); 
    end
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves cosmetic makeover
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    grid on;

    ylabel('$\pfa$','FontSize',Fontsize);
    xlabel('$\test$ = [ms]','FontSize',Fontsize);
    %hl = legend('$\pfa$', '$\pfaac$', '$\pfaoc$');
    %set(hl, 'position',[0.752 0.685 0.09 0.23]);
    hl = legend('IM', 'EM-AC', 'EM-OC');
    set(hl, 'position',[0.615 0.685 0.28 0.23]);
    set(gca,'FontSize',Fontsize);
    laprint(2, '../figures/fig_P_f_vs_est_time_diff_mu_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
end

if 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Second Approach Plot curves simulation analysis --  Optimum throughput vs est time    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu05_oc_sim2_app2.mat');

    %% ac, oc, mu = 0.05   
    Fontsize = 9;
    figure(3);
    scale = 1:1:length(N_est_sim);
    if 1
        %plot(N_est_sim(scale) * 1e-3, R_id_opt_sim(1,(scale)),'c-','LineWidth', 4);
        %hold on,
        plot(N_est_sim(scale) * 1e-3, R_ac_opt_sim(1,(scale)),'k+', 'LineWidth', 1);
        hold on,
    end
    plot(N_est_sim(scale) * 1e-3, R_oc_opt_sim(1,(scale)),'kx', 'LineWidth', 1);
    
    if 0  
        [temp index] = max(R_ac_opt_th);
        R_ac_opt_sen_opt_est = R_ac_opt_th(1,index);
        N_ac_opt_est = N_est_th(index);

    end
    [temp index] = max(R_oc_opt_sim);
    R_oc_opt_sen_opt_est = R_oc_opt_sim(1,index);
    N_oc_opt_est = N_est_sim(index);
    
    if 0
        hold on,
        plot(N_ac_opt_est * 1e-3, R_ac_opt_sen_opt_est, 'rs', 'LineWidth',2, 'HandleVisibility','off');
    end
    %hold on,
    %plot(N_oc_opt_est * 1e-3, R_oc_opt_sen_opt_est, 'rs', 'LineWidth',2);
    

    
    %% mu = 0.10 
    if 1
        load('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu10_oc_sim2_app2.mat');

        hold on,
        plot(N_est_sim(scale) * 1e-3, R_oc_opt_sim(1,(scale)),'kx', 'LineWidth', 1);
        
        [temp index] = max(R_oc_opt_sim);
        R_oc_opt_sen_opt_est = R_oc_opt_sim(1,index);
        N_oc_opt_est = N_est_sim(index);
        %hold on,
        %plot(N_oc_opt_est * 1e-3, R_oc_opt_sen_opt_est, 'rs', 'LineWidth',2, 'HandleVisibility','off');
    end
    
    %% mu = 0.15
    if 1
        load('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu15_oc_sim2_app2.mat');

        hold on,
        plot(N_est_sim(scale) * 1e-3, R_oc_opt_sim(1,(scale)),'kx', 'LineWidth', 1);
        
        [temp index] = max(R_oc_opt_sim);
        R_oc_opt_sen_opt_est = R_oc_opt_sim(1,index);
        N_oc_opt_est = N_est_sim(index);
        %hold on,
        %plot(N_oc_opt_est * 1e-3, R_oc_opt_sen_opt_est, 'rs', 'LineWidth',2);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   First Approach Plot curves simulation analysis --  Optimum throughput vs est time    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu05_oc_th2.mat');
    %% ac, oc, mu = 0.05   
    Fontsize = 9;
    figure(3);
    [value scale] = find(N_est_th >= 1000 & N_est_th <= 10000);
    if 1
        plot(N_est_th(scale) * 1e-3, R_id_opt_th(1,(scale)),'c-','LineWidth', 4);
        hold on,
        plot(N_est_th(scale) * 1e-3, R_ac_opt_th(1,(scale)),'-', 'LineWidth', 1);
        hold on,
    end
    hold on,
    plot(N_est_th(scale) * 1e-3, R_oc_opt_th(1,(scale)),'--', 'LineWidth', 1);
    
    if 1  
        [temp index] = max(R_ac_opt_th);
        R_ac_opt_sen_opt_est = R_ac_opt_th(1,index);
        N_ac_opt_est = N_est_th(index);

    end
    [temp index] = max(R_oc_opt_th);
    R_oc_opt_sen_opt_est = R_oc_opt_th(1,index);
    N_oc_opt_est = N_est_th(index);
    
    if 1
        hold on,
        plot(N_ac_opt_est * 1e-3, R_ac_opt_sen_opt_est, 'rs', 'LineWidth',1, 'HandleVisibility','off');
    end
    
    hold on,
    plot(N_oc_opt_est * 1e-3, R_oc_opt_sen_opt_est, 'rs', 'LineWidth',1);
    
    grid on;
    axis([min(N_est_th(scale) * 1e-3) max(N_est_th(scale) * 1e-3) min(R_oc_opt_th(1,(scale))) max(R_id_opt_th(1,(scale))) * 1.02]);
    
    
    
    %% mu = 0.10 
    if 1
        load('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu10_oc_th2.mat');

        hold on,
        plot(N_est_th(scale) * 1e-3, R_oc_opt_th(1,(scale)),'--', 'LineWidth', 1, 'HandleVisibility','off');

        [temp index] = max(R_oc_opt_th);
        R_oc_opt_sen_opt_est = R_oc_opt_th(1,index);
        N_oc_opt_est = N_est_th(index);
        hold on,
        plot(N_oc_opt_est * 1e-3, R_oc_opt_sen_opt_est, 'rs', 'LineWidth',1, 'HandleVisibility','off');
    end
    
    %% mu = 0.15
    if 1
        load('results_opt_thr_P_d_P_f_vs_est_time_diff_SNR_wo_nu_m10_mu15_oc_th2.mat');

        hold on,
        plot(N_est_th(scale) * 1e-3, R_oc_opt_th(1,(scale)),'--', 'LineWidth', 1);

        [temp index] = max(R_oc_opt_th);
        R_oc_opt_sen_opt_est = R_oc_opt_th(1,index);
        N_oc_opt_est = N_est_th(index);
        hold on,
        plot(N_oc_opt_est * 1e-3, R_oc_opt_sen_opt_est, 'rs', 'LineWidth',1);
    end
    
    
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves cosmetic makeover
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ylabel('$\trs(\test,\ttsen)$ [bits/sec/Hz]','FontSize',Fontsize);
    xlabel('$\test$ [ms]','FontSize',Fontsize);
    %hl = legend('$\trs$', '$\trsac$', '$\trsoc$','$\ttest$');
    %set(hl, 'position',[0.725 0.12 0.15 0.28]);
    hl = legend('IM', 'EM-AC', 'EM-OC','$\trs(\ttest,\ttsen)$');
    set(hl, 'position',[0.61 0.12 0.26 0.28]);
    set(gca,'FontSize',Fontsize);
    laprint(3, '../figures/fig_opt_thr_vs_est_time_diff_mu_AWGN_2app', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
end
