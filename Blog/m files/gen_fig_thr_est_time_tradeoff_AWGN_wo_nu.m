%                                
%        Description: This m-file considers the Estimation-Throughput tradeoff 
%        for the AWGN channel.  
%        (When an outage constraint is established over the power received at PR)        
%
%        Clearly the regulated power is distributed according non central
%        chi-2 distribution. However to simplify the analysis we
%        approximate the non central distribution against the Gamma
%        distribution.
%        
%        The simulation is performed to show the following:
%        1) Analyse the throughput vs estimation curves
%        2) Confirm the theoretical analysis vs simulation results
%
%        Created on: 04.08.15
%        Revision History: 04.08.15 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

sim = 00;                                            % Enable( = 1) to perform simulation, theoretical analysis
th = 00;                                             % disable ( = 0) to plot curves, data is read from a file
                                                     % Theoretical anaylsis is also included of simulation as
                                                     % numerical integration is involved
al = 00;                                             % An alternative approach that doesn't consider the power control
                                                     % and transmits a constant power, uses estimation samples to satisfy
                                                     % the outage constaint on the primary interference power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_tran = -00;                                        % Power transmitted by PR, the SNR received at ST can be 
                                                     % varried using this parameter
noise_power = -100;                                  % noise power -100 dBm
I_T = -110;                                          % Interference temperature -80 dBm
f_s = 1e6;                                           % 1 MHz one band
K = 0.1 * f_s;                                       % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha_true_p1 = 100;                                 % True Path loss between ST and PR  to achieve P_p = IT  
alpha_true_p2 = 100;                                 % True Path loss between PT and SR  
alpha_true_s = 080;                                  % True Path loss between ST and SR  to achieve P_p = IT  
epsilon = 0.01;                                      % Outage Probability constraint at PR
snr = 10^(P_tran/10) * 10^(-alpha_true_p1/10)...     % snr received at the PR
    / 10^(noise_power/10);
P_reg_max = 10^(-00/10);                             % Maximum Trasmit Power Constraint

E_s = 1;                                             % Pilot Energy 
N_s = 10;                                            % Number of Pilot Symobols 
sigma_est_error = 10^(noise_power/10)/N_s;           % Variance of the Estimation error for h_s 
dof_s = 1;                                           % Degree of freedom = 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if sim
    M = 1e4;                                  % Number of realizations 
    N_sim = ceil([0.5:0.5:20] * (1e-3 * f_s));    % N = Total of samples used for estimation = tau * f_s 
    P_reg_sim = zeros(1, length(N_sim));      % Regulated transmitted power at ST
    P_rcvd_sim = zeros(1, M);                 % Received Power at ST

    Exp_R_sim = zeros(1,length(N_sim));       % Expected throughput  
    num_tp = 10000;                            % Number of test points
    for i=1:length(N_sim)
        disp(strcat('N = ',num2str(N_sim(i)))); 
        parfor j=1:M
            % Estimate the received power P_reg with additive noise amplitude (zero mean and noise variance) 
            samples_st = sqrt(10^((P_tran)/10) * 10^(-alpha_true_p1/10)) * ones(1, N_sim(i)) + random('norm', 0, sqrt(10^(noise_power/10)), 1, N_sim(i));
            samples_sr = sqrt(10^((P_tran)/10) * 10^(-alpha_true_p2/10)) * ones(1, N_sim(i)) + random('norm', 0, sqrt(10^(noise_power/10)), 1, N_sim(i));
            P_rcvd_st_sim(j) = mean(samples_st.^2);               
            P_rcvd_sr_sim(j) = mean(samples_sr.^2);    
            alpha_s_est(j) = mean(normrnd(sqrt(10^(-alpha_true_s/10)), sqrt(10^(noise_power/10)), 1, N_s).^2);
        end            
            % Perform power control over the received samples, 
            % hist(P_reg_sim(i,:), 200)
            % print(gcf,'-dpdf','../figures/hist_P_reg')
            test_points = linspace(min(P_rcvd_st_sim), max(P_rcvd_st_sim), num_tp);
            for k=1:num_tp
                %(length(find(P_rcvd_sim < test_points(k)))/M)
                if (length(find(P_rcvd_st_sim < test_points(k)))/M) >= (1 - epsilon)   
                    P_reg_sim(i) = 10^(I_T/10) * 10^(P_tran/10) /...
                        (test_points(k) - 10^(noise_power/10));       
                    break;
                end
            end                   
            
            if 0
            
                C_1_ = calc_capacities(10^(-alpha_true_p2/10)/10^(noise_power/10),...
                    10^(-alpha_true_s/10), P_reg_sim(i), 10^(noise_power/10), N_sim(i), N_s);
            end
            C_1 = mean(log2(1  + P_reg_sim(i) * alpha_s_est ./P_rcvd_sr_sim));
            Exp_R_sim(i) = (K- N_sim(i))/K * C_1;     

            disp(strcat('Exp_R_sim(i) = ',num2str(Exp_R_sim(i)))); 
            disp(strcat('P_reg_sim(i) = ',num2str(P_reg_sim(i)))); 
    end
    save('results_thr_est_time_tradeoff_AWGN_wo_nu_snr_m00_e01_sim.mat');
    if ~th      % If theoretical analysis is already simulated the computation ends here.
        quit;
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th
   %% Theoretical parameters
    N_th = [1:2:25,28:4:20000];                           % N = Total of samples used for estimation = tau * f_s 
    P_reg_th = zeros(1,length(N_th));             % Expected power regulated by ST
    Exp_R_th = zeros(1,length(N_th));             % Expected throughput  

    
    for i=1:length(N_th)       
        disp(strcat('N_th = ',num2str(N_th(i))));             
           
       %% Determining the performance meterics --> Exp_R
       %% Expected values           
       
       %% Gamma Approximation to the non-central chi-squared distribution
       mean = 10^(noise_power/10) * (1 + snr);       
       var = (10^(noise_power/10))^2/N_th(i) * ( 2 + 4 * snr); 
       
       b = var/mean;
       a = mean/b;
       
       %% Determine the controlled power
       P_reg_th(i) =  min(P_reg_max, 10^(I_T/10) * 10^(P_tran/10) / (b * gammaincinv(epsilon, a,'upper')...
           - 10^(noise_power/10)));       
      
       C_1 = calc_capacities(10^(-alpha_true_p2/10)/10^(noise_power/10),...
                10^(-alpha_true_s/10), P_reg_th(i), 10^(noise_power/10), N_th(i), N_s);
       
       %% Expected Rate
        Exp_R_th(i) = (K- N_th(i))/K * C_1;
        
       %%
        disp(strcat('Exp_R_th(i) = ',num2str(Exp_R_th(i)))); 
        disp(strcat('P_reg_th(i) = ',num2str(P_reg_th(i)))); 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Finding optimum estimation time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimum throughput for the ideal model
    [Exp_R_opt index] = max(Exp_R_th);
    N_opt = N_th(index);

    save('results_thr_est_time_tradeoff_AWGN_wo_nu_snr_m00_e10_th.mat');
    %quit;
end

time = 10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_thr_est_time_tradeoff_AWGN_wo_nu_snr_m00_e01_th.mat');
%pbaspect([1 0.3 1]);
[values scope_th] = find(N_th <= time);
gap = 12;
downsampled = 1:ceil(length(scope_th)/gap):length(scope_th);
Fontsize = 8.5;
h2 = plot(N_th(scope_th) * 1e-3, Exp_R_th(scope_th), 'k-', 'LineWidth', 1);



%% Ideal Model
Exp_R_id = log2(1 + 10^(I_T/10) ./(10^(P_tran/10) * 10.^(-alpha_true_p1/10))...
    * 10^(-alpha_true_s/10)./ (10^(P_tran/10) *  (10^(-alpha_true_p2/10) + 10^(noise_power/10))));
hold on,
%h1 = plot(N_th(scope_th) * 1e-3, Exp_R_id * ones(1, length(N_th(scope_th))), 'c', 'LineWidth', 1);
%if 1
%    hold on;
%    plot(N_th(downsampled) * 1e-3, Exp_R_id * ones(1, length(N_th(downsampled))), 'cd', 'LineWidth', 1);
%    lp1 = plot(0,0,'-c','Marker','d', 'LineWidth', 1, 'visible','off');
%end
Y_min = min(Exp_R_th(scope_th));
Y_max = max(Exp_R_id);

load('results_thr_est_time_tradeoff_AWGN_wo_nu_snr_m00_e10_th.mat');
%pbaspect([1 0.3 1]);
%h5 = plot(N_th(scope_th) * 1e-3, Exp_R_th(scope_th), 'k-', 'LineWidth', 1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_thr_est_time_tradeoff_AWGN_wo_nu_snr_m00_e01_th.mat');
hold on,
h3 = plot(N_opt * 1e-3, Exp_R_opt, 'rs', 'LineWidth', 1);
load('results_thr_est_time_tradeoff_AWGN_wo_nu_snr_m00_e10_th.mat');
%hold on,
%plot(N_opt * 1e-3, Exp_R_opt, 'rs', 'LineWidth', 1);

axis([0 max(N_th(scope_th)) * 1e-3  1.5 Y_max  * 1.01]);
ylabel('$\rs(\test)$ [bits/sec/Hz]','FontSize',Fontsize);
xlabel('$\test$ [ms]','FontSize',Fontsize);
%hl = legend([h1 h2 h3 h4],'IM','EM', '$\trs(\ttau)$', 'Simulated');
%if 1
%    hl = legend([lp1 h2 h3 h4],'IM','EM', '$\trs(\ttau)$', 'Simulated');
%end
    %hl = legend([h1 h2 h3 h5 h4],'IM','EM, $\opc = 0.1$', 'EM, $\opc = 0.01$', '$\trs(\ttau)$', 'Simulated');
%set(hl, 'position',[0.725 0.12 0.15 0.31]);
%set(hl, 'position',[0.15 0.13 0.34 0.30]);
%set(hl, 'position',[0.15 0.13 0.4 0.35]);
set(gca,'FontSize',Fontsize);
laprint(1, '../figures/fig_thr_est_time_tradeoff_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');


