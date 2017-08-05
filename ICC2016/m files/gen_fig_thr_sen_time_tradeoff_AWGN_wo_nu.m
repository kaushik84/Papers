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
%        Created on: 10.09.15
%        Revision History: 10.09.15 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

sim = 00;                                  % Enable( = 1) to perform simulation, theoretical analysis
th = 00;                                  % disable ( = 0) to plot curves, data is read from a file
                                          % Theoretical anaylsis is also included of simulation as
                                          % numerical integration is involved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_full = 10^(+00/10);                     % Max Transmit power transmitted by ST 
P_p = 10^(+10/10);                        % Power transmitted by PT, the SNR received at ST 
                                          % varried using this parameter
noise_power = 10^(-100/10);               % noise power -100 dBm
f_s = 1e6;                                % 1 MHz one band
K = 0.1 * f_s;                            % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha_p_1 = 10^(-120/10);                 % True Path loss between ST and PR   
alpha_p_2 = 10^(-120/10);                 % True Path loss between PT and SR   
alpha_p_3 = 10^(-100/10);                 % True Path loss between ST and PR   
alpha_s = 10^(-090/10);                   % True Path loss between ST and SR 
P_H0 = 0.8;                               % Probability of Hypothesis 0
P_d_d = 0.90;                             % Constraint on Probability of detection P_d
theta_it = 10^(-110/10);                  % Interference temperature 
rho_pd = 0.05;                             % Outage contraint on the detection probability
rho_cont = 0.1;                           % Outage contraint on the controlled power


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sim
     M = 1e4;                                       % Number of realizations 
     %% Simulations parameters
     N_est_hp1_sim = 1000;                          % N_est_th = number of samples used for estimation = tau * f_s
     N_est_hp3_sim = 1000;                          % N_est_th = number of samples used for estimation = tau * f_s 
     N_sen_sim = 600:2000:20000;                    % N = Total of samples used for estimation = tau * f_s 

     %epsilon_id_sim = zeros(1,length(N_sen_sim));   % Ideal threshold 

     P_f_id_sim = zeros(1,length(N_sen_sim));       % Probability of false alarm at ST for the ideal case
     P_d_id_sim = zeros(1,length(N_sen_sim));       % Probability of detection at ST for the ideal case

     P_f_em_sim = zeros(1,length(N_sen_sim));       % Probability of false alarm at ST for the EM case
     P_d_em_sim = zeros(1,length(N_sen_sim));       % Probability of detection at ST for the EM case
     
     
     P_reg_id_sim = zeros(1,length(N_sen_sim));     % Regulated power at the ST for the ideal case
     P_reg_em_sim = zeros(1,length(N_sen_sim));     % Regulated power at the ST for the EM case

       
     R_id_sim = zeros(1,length(N_sen_sim));         % Throughput at SR ideal case
     R_em_sim = zeros(1,length(N_sen_sim));         % Throughput at SR em case

     num_tp = 5000;                                 % Test point for evaluating P_d, P_f
     
     for i=1:length(N_sen_sim)
        disp(strcat('N_sen_sim = ',num2str(N_sen_sim(i)))); 
        P_rcvd_hp1_est_sim = zeros(1,M);                % Rcvd power samples under Hypothesis 1 
        P_rcvd_hp1_sen_H1_sim = zeros(1,M);             % Rcvd power samples under Hypothesis 1 
        P_rcvd_hp1_sen_H0_sim = zeros(1,M);             % Rcvd power samples under Hypothesis 1 
        P_rcvd_hp3_est_sim = zeros(1,M);                % Rcvd power samples under Hypothesis 1 
        
        epsilon_id_sim = 0;
        epsilon_em_sim = 0;
        
        
       %% Sensing Phase for EM model
       %% Sensing of channel h_p1
        parfor k=1:M 
            % Listens the received power from PT without additive noise amplitude (zero mean and noise variance) 
            samples = random('norm', 0, sqrt(P_p * alpha_p_1), 1, N_sen_sim(i)) + random('norm', 0, sqrt(noise_power), 1, N_sen_sim(i));
            P_rcvd_hp1_sen_H1_sim(k) = mean(samples.^2);                                  
        end 
        
        parfor k=1:M
            % Listens the received power P_reg without additive noise amplitude (zero mean and noise variance) 
            samples = random('norm', 0, sqrt(noise_power), 1, N_sen_sim(i));
            P_rcvd_hp1_sen_H0_sim(k) = mean(samples.^2);                                  
        end 

       %% Determining the performance metric R for ideal case and EM case
       %% Threshold     
        threshold_tp = linspace(min(P_rcvd_hp1_sen_H1_sim), mean(P_rcvd_hp1_sen_H1_sim),...
            num_tp);
        for k = length(threshold_tp):-1:1
            if (length(find(P_rcvd_hp1_sen_H1_sim > threshold_tp(k)))/M) >= P_d_d   
                epsilon_id_sim = threshold_tp(k);
                break;
            end
        end
        
       %% False alram probability and probability of detection
        P_f_id_sim(i) = length(find(P_rcvd_hp1_sen_H0_sim > epsilon_id_sim))/M;
        P_d_id_sim(i) = length(find(P_rcvd_hp1_sen_H1_sim > epsilon_id_sim))/M;
 
        
       %% Determine the regulated power
        P_reg_id_sim(i) = min(theta_it / ((1 - P_H0) * alpha_p_3), P_full);
        
                
       %% Expected Rate
        C_0_sim = log2(1 + alpha_s * P_full / noise_power); 
        C_1_sim = log2(1 + alpha_s * P_full / (P_p * alpha_p_2 +...
            noise_power) );
        C_2_sim = log2(1 + alpha_s * P_reg_id_sim(i) / noise_power); 
        C_3_sim = log2(1 + alpha_s * P_reg_id_sim(i) / (P_p * alpha_p_2 +...
            noise_power) );        

        R_id_sim(i) = (K - N_sen_sim(i))/K * (P_H0 * (1 -  P_f_id_sim(i)) *...
            C_0_sim +  (1 - P_H0) * (1 - P_d_id_sim(i)) * C_1_sim +... 
            P_H0 * P_f_id_sim(i) * C_2_sim + (1 - P_H0) * P_d_id_sim(i) * C_3_sim);  

       %% Estimation Phase for EM model
       if N_sen_sim(i) >= N_est_hp1_sim 
           %% Estimation of channel h_p1
            parfor k=1:M
                % Estimate the received power from PT with additive noise amplitude (zero mean and noise variance) 
                samples = random('norm', 0, sqrt(P_p * alpha_p_1), 1, N_est_hp1_sim) + random('norm', 0, sqrt(noise_power), 1, N_est_hp1_sim);
                P_rcvd_hp1_est_sim(k) = mean(samples.^2);                                  
            end 
            %% Estimation of channel h_p3
            parfor k=1:M
                % Estimate the received power from PT with additive noise amplitude (zero mean and noise variance) 
                samples = random('norm', 0, sqrt(P_p * alpha_p_3), 1, N_est_hp3_sim)  + random('norm', 0, sqrt(noise_power), 1, N_est_hp3_sim);
                P_rcvd_hp3_est_sim(k) = mean(samples.^2);                                  
            end      


           %% Threshold 
           threshold_tp = linspace(epsilon_id_sim * .1, epsilon_id_sim,...
            num_tp);
            P_d_variations = zeros(1,M);
            for k = 1:num_tp
                P_d_variations(:) = gammainc(N_sen_sim(i) * threshold_tp(k)./...
                    (2 * P_rcvd_hp1_est_sim), N_sen_sim(i)/2, 'upper');
                if (length(find(P_d_variations <= P_d_d))/M >= rho_pd)
                    epsilon_em_sim = threshold_tp(k);

                    break;
                end
            end
            
            rcvd_energy_pt = (noise_power + alpha_p_1 * P_p);
            epsilon_em_sim = 4 * rcvd_energy_pt * gammaincinv(1 - rho_pd, N_est_hp1_sim/2, 'upper') *...
                gammaincinv(P_d_d, N_sen_sim(i)/2, 'upper')/ (N_est_hp1_sim * N_sen_sim(i)); 

           %% Determine the regulated power                
            P_reg_tp = linspace(P_full, P_full * 1e-2, num_tp);        
            for k=1:num_tp
                P_rcvd_variations = (1 - P_H0) * (P_rcvd_hp3_est_sim - noise_power)/P_p .*...
                    ((1 - P_d_variations) * P_full + P_d_variations * P_reg_tp(k));  
                if length(find(P_rcvd_variations >= theta_it))/M <= rho_cont
                    P_reg_em_sim(i) = P_reg_tp(k); 
                    break;
                end                
            end     


           %% Determine the regulated power
            P_f_em_sim(i) = length(find(P_rcvd_hp1_sen_H0_sim > epsilon_em_sim))/M;
           %% Expected probability of detection is calculated
           %% as P_d is random variable;
            P_d_em_sim(i) = mean(P_d_variations);    

           %% Expected Rate
            [C_0_sim, C_1_sim] = calc_capacities(alpha_p_2 * P_p / noise_power,... 
                alpha_s * P_full / noise_power,...
                noise_power, N_est_hp3_sim);      

            if (P_reg_em_sim(i) ~= 0)
                % with power regulation at the ST
                [C_2_sim, C_3_sim] = calc_capacities(alpha_p_2 * P_p / noise_power,... 
                    alpha_s * P_reg_em_sim(i) / noise_power,...
                    noise_power, N_est_hp3_sim);
            else
                C_2_sim = 0;
                C_3_sim = 0;
            end

            R_em_sim(i) = (K - (N_est_hp3_sim) - N_sen_sim(i))/K *...
                (P_H0 * (1 -  P_f_em_sim(i)) * C_0_sim +  (1 - P_H0) * (1 - P_d_em_sim(i)) * C_1_sim +... 
                P_H0 * P_f_em_sim(i) * C_2_sim + (1 - P_H0) * P_d_em_sim(i) * C_3_sim);    
	else
	    R_em_sim(i) = 0;	
	end		
            
        disp(strcat('epsilon_id_sim = ',num2str(epsilon_id_sim)));   
        disp(strcat('P_reg_id_sim = ',num2str(P_reg_id_sim(i))));         
        disp(strcat('R_id_sim = ',num2str(R_id_sim(i)))); 
        disp(strcat('epsilon_em_sim = ',num2str(epsilon_em_sim))); 
        disp(strcat('P_reg_em_sim = ',num2str(P_reg_em_sim(i))));         
        disp(strcat('R_em_sim = ',num2str(R_em_sim(i)))); 
    	save('results_thr_sen_time_tradeoff_AWGN_P_full_p00_P_p_p05_pd_10_wo_nu_sim2.mat');
    end
    if ~th      % If theoretical analysis is already simulated the computation ends here.
        quit;
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th
    
 
   %% Theoretical parameters
    N_est_hp1_th = 1000;                               % N_est_th = number of samples used for estimation = tau_est_hp1 * f_s 
    N_est_hp3_th = 1000;                               % N_est_th = number of samples used for estimation = tau_est_hp3 * f_s 
    N_sen_th = [5:100:N_est_hp1_th-1,N_est_hp1_th:100:20000];                           % N = includes the time in samples, for the ideal case their is no
                                                       % estimation time hence sensing time = N_th and the throughput will 
                                                       % attain a nonzero value. For the accurate, the energy walls sensing 
                                                       % time begin after the estimation of the rcvd energy is performed.
                                                       % N_sen = N_th - N_est
    %epsilon_id_th = zeros(1,length(N_sen_th));         % Ideal threshold 
    %epsilon_em_th = zeros(1,length(N_sen_th));         % Estimation Model 

    
    P_f_id_th = zeros(1,length(N_sen_th));             % False Alarm Probability at ST for the ideal case
    P_d_id_th = zeros(1,length(N_sen_th));             % Detection Probability at ST for the ideal case
    P_f_em_th = zeros(1,length(N_sen_th));             % False Alarm Probability at ST for the EM case
    P_d_em_th = zeros(1,length(N_sen_th));             % Detection Probability at ST for the EM case

    P_reg_id_th = zeros(1,length(N_sen_th));           % Regulated power at the ST for the ideal case
    P_reg_em_th = zeros(1,length(N_sen_th));           % Regulated power at the ST for the EM case

    
    R_id_th = zeros(1,length(N_sen_th));               % Throughput at SR ideal case
    
    num_tp = 5000;                                    % Number of test points 
    
    %% Determine the energy walls
    % Accurarte energy  
    rcvd_energy_pt = (noise_power + alpha_p_1 * P_p);      
    rcvd_energy_pr = (noise_power + alpha_p_3 * P_p);      
    parfor i=1:length(N_sen_th)
        epsilon_id_th = 0;
        epsilon_em_th = 0;
        
    	warning off;
        disp(strcat('N_sen_th = ',num2str(N_sen_th(i))));     
        
       %% Determining the performance meteric R for all the different cases
       %% Ideal Model        
       %% Threshold  
       epsilon_id_th = gammaincinv(P_d_d, N_sen_th(i)/2, 'upper') * 2 * rcvd_energy_pt/ N_sen_th(i); 
              
       %% Probability of false alarm and probability of dectection
       P_f_id_th(i) = gammainc(N_sen_th(i)/2 * epsilon_id_th/ noise_power', N_sen_th(i)/2, 'upper');           
       P_d_id_th(i) = gammainc(N_sen_th(i)/2 * epsilon_id_th/ rcvd_energy_pt, N_sen_th(i)/2, 'upper'); 
  
      %% Determine the regulated power
       P_reg_id_th(i) = min(theta_it / ((1 - P_H0) * alpha_p_3), P_full);
       
       %% Expected Rate
        C_0_th = log2(1 + alpha_s * P_full / noise_power); 
        C_1_th = log2(1 + alpha_s * P_full / (P_p * alpha_p_2 +...
            noise_power) );
        C_2_th = log2(1 + alpha_s * P_reg_id_th(i) / (P_p * alpha_p_2 +...
            noise_power) );
        C_3_th = log2(1 + alpha_s * P_reg_id_th(i) / (P_p * alpha_p_2 +...
            noise_power) );
        
        R_id_th(i) = (K - N_sen_th(i))/K * (P_H0 * (1 -  P_f_id_th(i)) * C_0_th +...
            (1 - P_H0) * (1 - P_d_id_th(i)) * C_1_th + P_H0 * P_f_id_th(i) * C_2_th +...
            (1 - P_H0) * P_d_id_th(i) * C_3_th);   
        
        
       %% Estimation Model        
        if N_sen_th(i) >= N_est_hp3_th        
            
           %% Threshold 
            epsilon_em_th = 4 * rcvd_energy_pt * gammaincinv(1 - rho_pd, N_est_hp1_th/2, 'upper') *...
                gammaincinv(P_d_d, N_sen_th(i)/2, 'upper')/ (N_est_hp1_th * N_sen_th(i)); 

           %% Probability of false alarm and probability of dectection            
            P_f_em_th(i) = gammainc(N_sen_th(i)/2 * epsilon_em_th/ noise_power', N_sen_th(i)/2, 'upper');           
           
            %% Evaluate the intergral using the distribution function
            if 1
                x_pts = linspace(0,1,num_tp);
                CDF_pts = 1 - gammainc(N_est_hp1_th * N_sen_th(i) * epsilon_em_th./...
                    (4 * rcvd_energy_pt * gammaincinv(x_pts, N_sen_th(i)/2, 'upper')),...
                    N_est_hp1_th/2, 'upper'); 
                PDF_pts = diff(CDF_pts);
                P_d_em_th(i) = sum(PDF_pts .* x_pts(2:end));
            
            %% Earlier Appraoch
            else          
                func_exp_pd = @(t) t .*...
                            exp(gammaln(N_sen_th(i)/2) - gammaln(N_est_hp1_th/2) - ...
                             N_est_hp1_th * N_sen_th(i) * epsilon_em_th./...
                             (4 * rcvd_energy_pt * gammaincinv(t, N_sen_th(i)/2, 'upper')) +...
                             gammaincinv(t, N_sen_th(i)/2, 'upper') + N_est_hp1_th/2 * log(N_est_hp1_th * N_sen_th(i) *...
                             epsilon_em_th./(4 * rcvd_energy_pt * gammaincinv(t, N_sen_th(i)/2, 'upper'))) -...
                             N_sen_th(i)/2 * log(gammaincinv(t, N_sen_th(i)/2, 'upper')));
                P_d_em_th(i) = integral(func_exp_pd, 0 ,1);    
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
                       P_reg_em_th(i) = P_reg_tp(k);
                        break;             
                    end
                end
            end 

           %% Expected Rate
            [C_0_th, C_1_th] = calc_capacities(alpha_p_2 * P_p / noise_power,... 
                alpha_s * P_full / noise_power,...
                noise_power, N_est_hp1_th);   


            if (P_reg_em_th(i) ~= 0)
                % with power regulation at the ST
                [C_2_th, C_3_th] = calc_capacities(alpha_p_2 * P_p/ noise_power,... 
                    alpha_s * P_reg_em_th(i)  / noise_power,...
                    noise_power, N_est_hp1_th);
            else
                C_2_th = 0;
                C_3_th = 0;
            end


            R_em_th(i) = (K - N_sen_th(i) -(N_est_hp3_th))/K * (P_H0 * (1 -  P_f_em_th(i)) * C_0_th +...
                (1 - P_H0) * (1 - P_d_em_th(i)) * C_1_th + P_H0 * P_f_em_th(i) * C_2_th +...
                (1 - P_H0) * P_d_em_th(i) * C_3_th);     
        else
            R_em_th(i) = 0;
        end
        
        disp(strcat('epsilon_id_th = ',num2str(epsilon_id_th)));         
        disp(strcat('P_reg_id_th = ',num2str(P_reg_id_th(i))));                 
        disp(strcat('R_id_th = ',num2str(R_id_th(i))));         
        disp(strcat('epsilon_em_th = ',num2str(epsilon_em_th))); 
        disp(strcat('P_reg_em_th = ',num2str(P_reg_em_th(i))));         
        disp(strcat('R_em_th = ',num2str(R_em_th(i))));
        
    end
    % Optimum throughput for the ideal model
    [R_id_opt index] = max(R_id_th);
    N_id_opt = N_sen_th(index);
    
    % Optimum throughput for the estimation model
    [R_em_opt index] = max(R_em_th);
    N_em_opt = N_sen_th(index);
    
    save('results_thr_sen_time_tradeoff_AWGN_P_full_p00_P_p_p10_pd_05_wo_nu_th2.mat');
    quit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Sensing throughput tradeoff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if 1
    load('results_thr_sen_time_tradeoff_AWGN_P_full_p00_P_p_p10_pd_10_wo_nu_th2.mat');
    load('results_thr_sen_time_tradeoff_AWGN_P_full_p00_P_p_p10_pd_10_wo_nu_sim2.mat');
    
    index_sim = 1:length(R_id_sim);
    index_th = 1:length(R_id_th) - 50;
    
    h1 = plot(N_sen_th * 1e-3, R_id_th, 'c', 'LineWidth',2);
    hold on,
    plot(N_id_opt * 1e-3, R_id_opt, 'rs', 'HandleVisibility','off');


    hold on,
    h2 = plot(N_sen_th * 1e-3, R_em_th, 'b', 'LineWidth',1);
    hold on,
    h3 = plot(N_em_opt * 1e-3, R_em_opt, 'rs', 'HandleVisibility','off');

    if 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Plot curves simulation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hold on,
        h4 = plot(N_sen_sim * 1e-3, R_id_sim, 'ko', 'LineWidth',1);
        hold on,
        plot(N_sen_sim * 1e-3, R_em_sim, 'ko', 'LineWidth',1);
    end 
end

load('results_thr_sen_time_tradeoff_AWGN_P_full_p00_P_p_p10_pd_05_wo_nu_th2.mat');
load('results_thr_sen_time_tradeoff_AWGN_P_full_p00_P_p_p10_pd_05_wo_nu_sim2.mat');



hold on,
plot(N_sen_th * 1e-3, R_em_th, 'b', 'LineWidth', 1);
hold on,
plot(N_em_opt * 1e-3, R_em_opt, 'rs', 'HandleVisibility','off');

if 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves simulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold on,
    plot(N_sen_sim * 1e-3, R_em_sim, 'ko', 'LineWidth',1);
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves cosmetic makeover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fontsize = 8;
grid on;
axis([0 max(N_sen_th(index_th)) * 1e-3 1.8  max(R_id_opt) * 1.01]);
ylabel('$\rs$ [bits/sec/Hz]','FontSize',Fontsize);
xlabel('$\tsen$ [ms]','FontSize',Fontsize);
hl = legend([h1 h2 h3 h4],'$\rs$, IM', '$\rs$, EM',...
    '$\rs(\ttsen)$', 'Simulated');
set(hl, 'Location', 'SouthEast', 'FontSize', Fontsize);
set(hl, 'position',[0.58 0.12 0.31 0.36]);
set(gca,'FontSize',Fontsize);
laprint(1, '../figures/fig_thr_sen_time_tradeoff_AWGN', 'options', 'factory', 'width', 6, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');
