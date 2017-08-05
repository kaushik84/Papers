%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: This m-file considers the Throughput-Sesning time 
%        tradeoff for the AWGN channel. In this m-script the short term analysis 
%        is investigated that is, considering invidual frame. 
%        Hence, the SNR is considered constant for the time interval 
%        (tau_est + tau_sen). 
%        This assumption however will become weak for lower SNR, because low
%        SNR regimes presents high tau_sen, leading to a low probability
%        that the channel remains constant.
%
%
%        Created on: 24.09.15
%        Revision History: 24.09.15 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

th = 00;                                  % disable ( = 0) to plot curves, data is read from a file
                                          % Theoretical anaylsis is also included of simulation as
                                          % numerical integration is involved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_full = 10^(+00/10);                     % Max Transmit power transmitted by ST 
P_p = 10.^(+10/10);                       % Power transmitted by PT, the SNR received at ST 
                                          % varried using this parameter
noise_power = 10^(-100/10);               % noise power -100 dBm
f_s = 1e6;                                % 1 MHz one band
K = 0.1 * f_s;                            % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha_p_1 = 10^(-120/10);                 % True Path loss between ST and PR   
alpha_p_2 = 10^(-120/10);                 % True Path loss between PT and SR   
alpha_p_3 = 10.^([-(110:-0.3:92),-91 -90]/10);        % True Path loss between ST and PR   
alpha_s = 10^(-090/10);                   % True Path loss between ST and SR 
P_H0 = 0.8;                               % Probability of Hypothesis 0
P_d_d = 0.90;                             % Constraint on Probability of detection P_d
theta_it = 10^(-110/10);                  % Interference temperature 
rho_pd = 0.05;                             % Outage contraint on the detection probability
rho_cont = 0.1;                           % Outage contraint on the controlled power


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th
    
 
   %% Theoretical parameters
    N_est_hp1_th = 1000;                                      % N_est_th = number of samples used for estimation = tau_est_hp1 * f_s 
    N_est_hp3_th = 1000;                                      % N_est_th = number of samples used for estimation = tau_est_hp3 * f_s 
    N_sen_th = [1:1:9,10:10:90,100:20:980,1000:50:6000];      % N = includes the time in samples, for the ideal case their is no
                                                              % estimation time hence sensing time = N_th and the throughput will 
                                                              % attain a nonzero value. For the accurate, the energy walls sensing 
                                                              % time begin after the estimation of the rcvd energy is performed.
                                                              % N_sen = N_th - N_est


    
    P_f_id_th = zeros(length(alpha_p_3),length(N_sen_th));             % False Alarm Probability at ST for the ideal case
    P_d_id_th = zeros(length(alpha_p_3),length(N_sen_th));             % Detection Probability at ST for the ideal case
    P_f_em_th = zeros(length(alpha_p_3),length(N_sen_th));             % False Alarm Probability at ST for the EM case
    P_d_em_th = zeros(length(alpha_p_3),length(N_sen_th));             % Detection Probability at ST for the EM case

    P_reg_id_th = zeros(length(alpha_p_3),length(N_sen_th));           % Regulated power at the ST for the ideal case
    P_reg_em_th = zeros(length(alpha_p_3),length(N_sen_th));           % Regulated power at the ST for the EM case

    
    R_id_th = zeros(length(alpha_p_3),length(N_sen_th));               % Throughput at SR ideal case
    
    R_opt_id_th = zeros(1, length(alpha_p_3));               % Throughput at SR @ optimum sensing time for ideal case
    R_opt_em_th = zeros(1, length(alpha_p_3));               % Throughput at SR @ optimum sensing time for EM case
        
    
    P_f_opt_id_th = zeros(1, length(alpha_p_3));             % Prob. of false alarm at ST @ optimum sensing time for ideal case 
    P_f_opt_em_th = zeros(1, length(alpha_p_3));             % Prob. of false alarm at ST @ optimum sensing time for EM case
    
    P_d_opt_id_th = zeros(1, length(alpha_p_3));             % Prob. of detection at ST @ optimum sensing time for ideal case 
    P_d_opt_em_th = zeros(1, length(alpha_p_3));             % Prob. of detection at ST @ optimum sensing time with average constraint 

    P_reg_opt_id_th = zeros(1, length(alpha_p_3));           % Optimum Controlled Power at ST @ optimum sensing time for ideal case 
    P_reg_opt_em_th = zeros(1, length(alpha_p_3));           % Optimum Controlled Power at ST @ optimum sensing time with average constraint 
        
        
    N_sen_opt_id_th = zeros(1, length(alpha_p_3));            % Optimum sensing time for ideal case
    N_sen_opt_em_th = zeros(1, length(alpha_p_3));            % Optimum sensing time for EM case
    
    
    num_tp = 1000;                                     % Number of test points 
    

    
    for j=1:length(alpha_p_3)          
  
        disp(strcat('alpha_p_3 = ', num2str(10 * log10(alpha_p_3(j)))));     
        %alpha_p_3_local = alpha_p_3(j);
        parfor i=1:length(N_sen_th)            
            warning off;
            % Received power  
            rcvd_energy_pt = (noise_power + alpha_p_1 * P_p);      
            rcvd_energy_pr = (noise_power + alpha_p_3(j) * P_p);
            %disp(strcat('N_sen_th = ', num2str(N_sen_th(i))));     

           %% Determining the performance meteric R for all the different cases
           %% Ideal Model        
           %% Threshold  
            epsilon_id_th = gammaincinv(P_d_d, N_sen_th(i)/2, 'upper') * 2 * rcvd_energy_pt/ N_sen_th(i); 


           %% Probability false alarm and probability of dectection
            P_f_id_th(j,i) = gammainc(N_sen_th(i)/2 * epsilon_id_th/ noise_power', N_sen_th(i)/2, 'upper');           
            P_d_id_th(j,i) = gammainc(N_sen_th(i)/2 * epsilon_id_th/ rcvd_energy_pt, N_sen_th(i)/2, 'upper');        
   

           %% Determine the regulated power
            P_reg_id_th(j,i) = min(theta_it / ((1 - P_H0) * alpha_p_3(j)), P_full);


           %% Expected Rate
            C_0_th = log2(1 + alpha_s * P_full / noise_power); 
            C_1_th = log2(1 + alpha_s * P_full / (P_p * alpha_p_2 +...
                noise_power) );

            C_2_th = log2(1 + alpha_s * P_reg_id_th(j,i) / (P_p * alpha_p_2 +...
                noise_power) );
            C_3_th = log2(1 + alpha_s * P_reg_id_th(j,i) / (P_p * alpha_p_2 +...
                noise_power) );

            R_id_th(j,i) = (K - N_sen_th(i))/K * (P_H0 * (1 -  P_f_id_th(j,i)) * C_0_th +...
                (1 - P_H0) * (1 - P_d_id_th(j,i)) * C_1_th + P_H0 * P_f_id_th(j,i) * C_2_th +...
                (1 - P_H0) * P_d_id_th(j,i) * C_3_th);            

          %% Estimation Model 
           if N_sen_th(i) >= N_est_hp1_th
              %% Threshold              
                epsilon_em_th = 4 * rcvd_energy_pt * gammaincinv(1 - rho_pd, N_est_hp1_th/2, 'upper') *...
                    gammaincinv(P_d_d, N_sen_th(i)/2, 'upper')/ (N_est_hp1_th * N_sen_th(i));  



                P_f_em_th(j,i) = gammainc(N_sen_th(i)/2 * epsilon_em_th/ noise_power', N_sen_th(i)/2, 'upper');           
            
               %% Evaluate the intergral using the distribution function                
                if 1
                    x_pts = linspace(0,1,num_tp);
                    CDF_pts = 1 - gammainc(N_est_hp1_th * N_sen_th(i) * epsilon_em_th./...
                        (4 * rcvd_energy_pt * gammaincinv(x_pts, N_sen_th(i)/2, 'upper')),...
                        N_est_hp1_th/2, 'upper'); 
                        PDF_pts = diff(CDF_pts);
                    P_d_em_th(j,i) = sum(PDF_pts .* x_pts(2:end));
              %% Earlier Appraoch
                else                  
                    func_exp_pd = @(t) t .*...
                                exp(gammaln(N_sen_th(i)/2) - gammaln(N_est_hp1_th/2) - ...
                                 N_est_hp1_th * N_sen_th(i) * epsilon_em_th./...
                                 (4 * rcvd_energy_pt * gammaincinv(t, N_sen_th(i)/2, 'upper')) +...
                                 gammaincinv(t, N_sen_th(i)/2, 'upper') + N_est_hp1_th/2 * log(N_est_hp1_th * N_sen_th(i) *...
                                 epsilon_em_th./(4 * rcvd_energy_pt * gammaincinv(t, N_sen_th(i)/2, 'upper'))) -...
                                 N_sen_th(i)/2 * log(gammaincinv(t, N_sen_th(i)/2, 'upper')));
                    P_d_em_th(j,i) = integral(func_exp_pd, 0 ,1);  
                end

               %% Determine the regulated power           
                result = zeros(1, num_tp);
                if 1
                    %gammainc(theta_it./((1 - t) * P_full  + t * test_points(k)) + noise_power,...
                             %N_est_hp3_th/2) .*;
                             % 
                    P_reg_tp = linspace(P_full, P_full * 1e-2, num_tp);
                    for k = 1:num_tp
                        func_dis_preg = @(t) (gammainc((theta_it * P_p ./((1 - P_H0) *...
                             ((1  - t) * P_full + t * P_reg_tp(k)))...
                            + noise_power) * N_est_hp3_th/ (2 * rcvd_energy_pr), N_est_hp3_th/2, 'upper')) .*...
                            exp(gammaln(N_sen_th(i)/2) - gammaln(N_est_hp1_th/2) - ...
                             N_est_hp1_th * N_sen_th(i) * epsilon_em_th./...
                             (4 * rcvd_energy_pt * gammaincinv(t, N_sen_th(i)/2, 'upper')) +...
                             gammaincinv(t, N_sen_th(i)/2, 'upper') + N_est_hp1_th/2 * log(N_est_hp1_th * N_sen_th(i) *...
                             epsilon_em_th./(4 * rcvd_energy_pt * gammaincinv(t, N_sen_th(i)/2, 'upper'))) -...
                             N_sen_th(i)/2 * log(gammaincinv(t, N_sen_th(i)/2, 'upper')));
                        result(k) = integral(func_dis_preg, 0, 1);
                        if result(k) <= rho_cont
                           P_reg_em_th(j,i) = P_reg_tp(k);
                            break;             
                        end
                    end
                end

               %% Expected Rate
                [C_0_th, C_1_th] = calc_capacities(alpha_p_2 * P_p / noise_power,... 
                    alpha_s * P_full / noise_power,...
                    noise_power, N_est_hp3_th);      

                if (P_reg_em_th(j,i) ~= 0)
                % with power regulation at the ST
                    [C_2_th, C_3_th] = calc_capacities(alpha_p_2 * P_p/ noise_power,... 
                    alpha_s * P_reg_em_th(j,i)  / noise_power,...
                    noise_power, N_est_hp3_th);
                else
                    C_2_th = 0;
                    C_3_th = 0;
                end


                R_em_th(j,i) = (K - N_sen_th(i) -( N_est_hp3_th))/K * (P_H0 * (1 -  P_f_em_th(j,i)) * C_0_th +...
                    (1 - P_H0) * (1 - P_d_em_th(j,i)) * C_1_th + P_H0 * P_f_em_th(j,i) * C_2_th +...
                    (1 - P_H0) * P_d_em_th(j,i) * C_3_th);   
           else
           
                 R_em_th(j,i) = 0;
           end

            %disp(strcat('P_reg_id_th = ',num2str(P_reg_id_th(j,i))));                 
            %disp(strcat('R_id_th = ',num2str(R_id_th(j,i))));         
            %disp(strcat('P_reg_em_th = ',num2str(P_reg_em_th(j,i))));         
            %disp(strcat('R_em_th = ',num2str(R_em_th(j,i))));

        end
        % Throughput, prob. of false alarm, prob. of detection 
        % @ optimum sensing time for the IM
        [temp_id index_id] = max(R_id_th(j,:));
        R_opt_id_th(j) = R_id_th(j,index_id);
        P_d_opt_id_th(j) = P_d_id_th(j,index_id);
        P_f_opt_id_th(j) = P_f_id_th(j,index_id);
        P_reg_opt_id_th(j) = P_reg_id_th(j,index_id);
        N_sen_opt_id_th(j) = N_sen_th(index_id);

        % Throughput, prob. of false alarm, prob. of detection 
        % @ optimum sensing time for the EM
        [temp_em index_em] = max(R_em_th(j,:));
        R_opt_em_th(j) = R_em_th(j,index_em);
        P_d_opt_em_th(j) = P_d_em_th(j,index_em);
        P_f_opt_em_th(j) = P_f_em_th(j,index_em);
        P_reg_opt_em_th(j) = P_reg_em_th(j,index_em);
        N_sen_opt_em_th(j) = N_sen_th(index_em);
        
        disp(strcat('N_sen_opt_id_th = ',num2str(N_sen_opt_id_th(j))));                 
        disp(strcat('P_reg_opt_id_th = ',num2str(P_reg_opt_id_th(j))));                 
        disp(strcat('R_opt_id_th = ',num2str(R_opt_id_th(j))));         
        disp(strcat('N_sen_opt_em_th = ',num2str(N_sen_opt_em_th(j))));                 
        disp(strcat('P_reg_opt_em_th = ',num2str(P_reg_opt_em_th(j))));         
        disp(strcat('R_opt_em_th = ',num2str(R_opt_em_th(j))));
	
	save('results_opt_thr_vs_SNR_AWGN_P_full_p00_P_p_pd_05_wo_nu_th2.mat');

    end
    quit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Throughput vs SNR @ optimum sensing time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
%% P_d = 0.10
load('results_opt_thr_vs_SNR_AWGN_P_full_p00_P_p_pd_10_wo_nu_th2.mat');
Fontsize = 8;
% * P_p/noise_power
plot(10 *log10(alpha_p_3), R_opt_id_th, 'c-', 'LineWidth',2);
hold on,
plot(10 *log10(alpha_p_3), R_opt_em_th, 'b-', 'LineWidth', 1);


%% P_d = 0.05
if 0
    load('results_opt_thr_vs_SNR_AWGN_P_full_p00_P_p_pd_05_wo_nu_th2.mat');
    Fontsize = 8;
    % * P_p/noise_power
    hold on,
    plot(10 *log10(alpha_p_3), R_opt_em_th, 'b-', 'LineWidth', 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves cosmetic makeover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid on;
axis([min(10 *log10(alpha_p_3))... 
    max(10 *log10(alpha_p_3))...
    min(R_opt_em_th) * 0.87  max(R_opt_id_th) * 1.02]);
ylabel('$\rs(\testpt, \testptsr, \testpr, \ttsen)$ [bits/sec/Hz]','FontSize',Fontsize);
xlabel('$\phpth$ [dB]','FontSize',Fontsize);

%hl = legend('$\trs$', '$\trsac$', '$\trsoc$');
%set(hl, 'position',[0.725 0.12 0.15 0.23]);
hl = legend('IM', 'EM');
set(hl, 'Location', 'SouthEast', 'FontSize', Fontsize);
set(hl, 'position',[0.68 0.71 0.22 0.21]);
set(gca,'FontSize',Fontsize);
laprint(1, '../figures/fig_opt_thr_vs_SNR_AWGN', 'options', 'factory', 'width', 6, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Power control vs Channel PT-PR @ optimum sensing time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% P_d = 0.05
if 0
    load('results_opt_thr_vs_SNR_AWGN_P_full_p00_P_p_pd_05_wo_nu_th2.mat');
    figure(2);
    Fontsize = 8;
    % * P_p/noise_power

    plot(10 *log10(alpha_p_3), 10 * log10(P_reg_opt_em_th), 'b-', 'LineWidth', 1);
end

%% P_d = 0.10
load('results_opt_thr_vs_SNR_AWGN_P_full_p00_P_p_pd_10_wo_nu_th2.mat');
figure(2);
Fontsize = 8;
% * P_p/noise_power
[value index ]= find(10 * log10(P_reg_opt_em_th) == -Inf) 
P_reg_opt_em_th(index) = 0.01 * ones(1,length(index));
%hold on,
plot(10 *log10(alpha_p_3), 10 * log10(P_reg_opt_id_th), 'c-', 'LineWidth',2);
hold on,
plot(10 *log10(alpha_p_3), 10 * log10(P_reg_opt_em_th), 'b-', 'LineWidth', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves cosmetic makeover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid on;
axis([min(10 *log10(alpha_p_3))... 
    max(10 *log10(alpha_p_3))...
    -24  1]);
ylabel('$\preg$ [dBm]','FontSize',Fontsize);
xlabel('$\phpth$ [dB]','FontSize',Fontsize);

%hl = legend('$\trs$', '$\trsac$', '$\trsoc$');
%set(hl, 'position',[0.725 0.12 0.15 0.23]);
hl = legend('IM', 'EM');
set(hl, 'Location', 'SouthEast', 'FontSize', Fontsize);
set(hl, 'position',[0.68 0.71 0.22 0.21]);
set(gca,'FontSize',Fontsize);
laprint(2, '../figures/fig_opt_cont_power_vs_SNR_AWGN', 'options', 'factory', 'width', 6, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');