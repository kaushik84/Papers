%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: This m-file considers the Throughput-Estimation time 
%        tradeoff for the AWGN channel.        
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
%        Created on: 17.09.14 --> Noise uncertainty added    
%        Revision History: 17.09.14 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

sim = 1;                                  % Enable( = 1) to perform simulation, theoretical analysis
th = 0;                                   % disable ( = 0) to plot curves, data is read from a file
                                          % Theoretical anaylsis is also included of simulation as
                                          % numerical integration is involved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu = 3;                                   % Noise uncertainty upper  
P_tran = -00;                             % Power transmitted by PR, the SNR received at ST can be 
                                          % varried using this parameternoise_power = -100;  
noise_power = -100;                       % Actual noise power -100 dBm
noise_power_p = -100 + nu;                % noise power at PR with nu noise power + nu dBm
noise_power_s = noise_power;              % noise power at PR without nu noise power dBm
I_T = -110;                               % Interference temperature -80 dBm
f_s = 1e6;                                % 1 MHz one band
K = 0.1 * f_s;                            % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha_true_p = 100;                       % True Path loss between ST and PR  to achieve P_p = IT  
alpha_true_s = 080;                       % True Path loss between ST and SR  to achieve P_p = IT  
mu_P = 0.025;                             % Accuracy at PR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if sim
    M = 1e4;                                  % Number of realizations 
    N_sim = ceil([1:4:20] * (1e-3 * f_s));    % N = Total of samples used for estimation = tau * f_s 
    P_reg_sim = zeros(length(N_sim), M);      % Regulated transmitted power at ST
    P_p_sim = zeros(length(N_sim), M);        % Power received at the PR
    R_sim = zeros(length(N_sim), M);          % Throughput at the SR
    P_rcvd_sim = zeros(1, M);                 % Received Power at ST

    Exp_R_sim = zeros(1,length(N_sim));       % Expected throughput  
    Exp_P_reg_sim = zeros(1,length(N_sim));   % Expected power regulated by ST
    Exp_P_p_sim = zeros(1,length(N_sim));     % Expected power received at PR
    PC_P_p_sim = zeros(1,length(N_sim));      % Confidence probability power at PR  
    PC_R_sim = zeros(1,length(N_sim));        % Condidence probability rate at SR
    for i=1:length(N_sim)
        disp(strcat('N = ',num2str(i))); 
        parfor j=1:M
            % Estimate the received power P_reg with additive noise amplitude (zero mean and noise variance) 
            samples = sqrt(10^((P_tran)/10) * 10^(-alpha_true_p/10)) * ones(1, N_sim(i)) + random('norm', 0, sqrt(10^(noise_power_p/10)), 1, N_sim(i));
            % Checkpoint 
            % Here the samples must be Gaussian distributed with mean = 10^(P_reg/10), 
            % examine by ploting the histogram, hist(samples, 100)
            P_rcvd_sim(i,j) = mean(samples.^2);                      
            
        end
            nc_sim = zeros(1, M);
            % Determine the scaling factor, which is the Expected Received
            % power at PR is scaled to meet the interference temperature
            nc_sim = 1 / mean(10^(-alpha_true_p/10) ./ P_rcvd_sim(i,:)); 
            
            % Perform power control over the received samples, 
            % P_reg_sim  = nc 
            % hist(P_reg_sim(i,:), 50)
            P_reg_sim(i,:) =  10^(I_T/10) * nc_sim ./ P_rcvd_sim(i,:);           
            
            % Determine P_p, P_s and R from P_reg after imposing path loss
            % hist(P_p_sim(i,:), 50)
            P_p_sim(i,:) = 10^(-alpha_true_p/10) * P_reg_sim(i,:);
           
            % hist(R_sim(i,:), 50)
            R_sim(i,:) = (K- N_sim(i))/K * log2(1 + 10^(-alpha_true_s/10) * P_reg_sim(i,:) / 10^(noise_power_s/10));     
            
           %% Determining the performance meterics --> Exp_R, CP_power, CP_rate            
            % Expected values
            Exp_P_reg_sim(i) = mean(P_reg_sim(i,:));
            Exp_R_sim(i) = mean(R_sim(i,:));
            Exp_P_p_sim(i) = mean(P_p_sim(i,:));
            
            % PC power
            upper_limit = Exp_P_p_sim(i) * (1 + mu_P);
            lower_limit = Exp_P_p_sim(i) * (1 - mu_P);
            PC_P_p_sim(i) = 1/M * length(find((P_p_sim(i,:) < upper_limit) & (P_p_sim(i,:) > lower_limit)));
            
            disp(strcat('Exp_R_sim(i) = ',num2str(Exp_R_sim(i)))); 
            disp(strcat('PC_power_sim(i) = ',num2str(PC_P_p_sim(i)))); 
    end
    save('results_thr_est_time_tradeoff_AWGN_nu_u_snr_m00_sim.mat');
    if ~th      % If theoretical analysis is already simulated the computation ends here.
        quit;
    end
end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Theoretical analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th
   %% Theoretical parameters
    N_th = ceil([0.5:0.1:20] * (1e-3 * f_s));   % N = Total of samples used for estimation = tau * f_s 
    Exp_P_reg_th = zeros(1,length(N_th));     % Expected power regulated by ST
    PC_P_p_th = zeros(1,length(N_th));        % Confidence probability power at PR  
    PC_R_th = zeros(1,length(N_th));          % Condidence probability rate at SR
    nc_th = 1;

    for i=1:length(N_th)       
        disp(strcat('N = ',num2str(i)));             
     
        nc_th = (10^(noise_power_p/10)/10^(-alpha_true_p/10) +  10^(P_tran/10));
        
       %% Determining the performance meterics --> Exp_R, CP_power, CP_rate
       %% Expected values            
        Exp_P_p_th = 10^(I_T/10);
       
        PC_P_p_th(i) = marcumq(sqrt(10^(-alpha_true_p/10) * 10^(P_tran/10) * ...
            N_th(i)/(10^(noise_power_p/10))),...
            sqrt(10^(I_T/10) * nc_th * N_th(i) * 10^(-alpha_true_p/10)/...
            (10^(noise_power_p/10) * (1 + mu_P) * Exp_P_p_th)), ceil(N_th(i)/2 - 1)) -... 
        marcumq(sqrt(10^(-alpha_true_p/10) * 10^(P_tran/10) * ...
            N_th(i)/(10^(noise_power_p/10))),...
            sqrt(10^(I_T/10) * nc_th * N_th(i) * 10^(-alpha_true_p/10)/...
            (10^(noise_power_p/10) * (1 - mu_P) * Exp_P_p_th)), ceil(N_th(i)/2 - 1));

       %% Expected Rate
        Exp_R_th(i) = (K- N_th(i))/K * log2(1 + 10^(I_T/10)/(10^(P_tran/10) *  10^(-alpha_true_p/10) +...
            10^(noise_power_p/10)) * nc_th * 10^(-alpha_true_s/10)  / 10^(noise_power_s/10));
        disp(strcat('Exp_R_th(i) = ',num2str(Exp_R_th(i)))); 
        disp(strcat('PC_power_th(i) = ',num2str(PC_P_p_th(i)))); 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Finding optimum sensing time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N_opt = ([0.04:0.04:20] * (1e-3 * f_s));
    epsilon_R = 0.95;                                       % Probability constraint
    epsilon_P = 0.95;                                       % Probability constraint
    N_optimum = 0;


    for i = length(N_opt):-1:1
        Exp_P_p_opt = 10^(I_T/10);

        % Confidence probability power at PR
        PC_P_p_opt = marcumq(sqrt(10^(-alpha_true_p/10) * 10^(P_tran/10) * ...
            N_opt(i)/(10^(noise_power_p/10))),...
            sqrt(10^(I_T/10) * nc_th * N_opt(i) * 10^(-alpha_true_p/10)/...
            (10^(noise_power_p/10) * (1 + mu_P) * Exp_P_p_th)), ceil(N_opt(i)/2 - 1)) -... 
        marcumq(sqrt(10^(-alpha_true_p/10) * 10^(P_tran/10) * ...
            N_opt(i)/(10^(noise_power_p/10))),...
            sqrt(10^(I_T/10) * nc_th * N_opt(i) * 10^(-alpha_true_p/10)/...
            (10^(noise_power_p/10) * (1 - mu_P) * Exp_P_p_th)), ceil(N_opt(i)/2 - 1));

        % Optimum PC        
        if (PC_P_p_opt < epsilon_P)
            N_optimum = N_opt(i);        
            break;
        end
        disp(strcat('PC_P_p_opt = ',num2str(PC_P_p_opt))); 
    end 
    Exp_R_opt = (K- N_optimum)/K * log2(1 + 10^(I_T/10)/(10^(P_tran/10) *  10^(-alpha_true_p/10) +...
            10^(noise_power_p/10)) * nc_th * 10^(-alpha_true_s/10)  / 10^(noise_power_s/10));
    save('results_thr_est_time_tradeoff_AWGN_nu_u_snr_m00_th.mat');
%    quit;
end
load('results_thr_est_time_tradeoff_AWGN_nu_u_snr_m10_th.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fontsize = 18;
PC_opt = PC_P_p_opt;
[ax h1 h2] = plotyy(N_th * 1e-3, Exp_R_th, N_th * 1e-3, PC_P_p_th);
set(h1,'LineWidth',1.5);
set(h2,'LineWidth',1.5);
%set(ax,'NextPlot','add');
%h3 = plot(ax(2), N_th * 1e-3, PC_R_th,'LineWidth',1.5,'LineStyle','-.');
set(ax,'NextPlot','add');
h4 = plot(ax(2), [N_optimum,max(N_opt)] * 1e-3, PC_opt * ones(1, 2),'LineStyle','--','LineWidth',1.5);
set(ax,'NextPlot','add');
plot(ax(2), N_optimum * ones(1, 2) * 1e-3, [PC_opt, min(min(PC_P_p_th, min(PC_R_th)))],...
    'LineWidth',1.5,'LineStyle','--', 'HandleVisibility','off');
set(ax,'NextPlot','add');
h5 = plot(ax(1), N_optimum * 1e-3, Exp_R_opt, 'o','LineWidth',1.5);

load('results_thr_est_time_tradeoff_AWGN_nu_u_snr_m10_sim.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(ax,'NextPlot','add');
h6 = plot(ax(1), N_sim * 1e-3, Exp_R_sim, 'r*', 'HandleVisibility','off');
set(ax,'NextPlot','add');
plot(ax(2), N_sim * 1e-3, PC_P_p_sim, 'r*', 'HandleVisibility','off');

grid on;
axis(ax(1), [min(N_th * 1e-3) max(N_th * 1e-3) min(Exp_R_th) max(Exp_R_th)]);
axis(ax(2), [min(N_th * 1e-3) max(N_th * 1e-3) min(PC_P_p_th) 1]);
ylabel(ax(1),'R_s = [bits/sec/Hz]','FontSize',Fontsize+2);
ylabel(ax(2),'Probability','FontSize',Fontsize+2);
xlabel(ax(1),'\tau = [ms]','FontSize',Fontsize+2);
xlabel(ax(2), '\tau = [ms]','FontSize',Fontsize+2);
hl = legend([h1;h2;h4;h5;h6],'R', 'PC', 'PC = \epsilon', 'max \tau', 'sim');
set(hl, 'Position', [.7, .5, .1, 0.2], 'FontSize', Fontsize + 2);
set(ax(1),'FontSize',Fontsize);
set(ax(2),'FontSize',Fontsize);
Ytick_R =  round(linspace(min(Exp_R_th), max(Exp_R_th), 10)*100)/100;
Ytick_PC = round(linspace(min(PC_P_p_th), 1, 10)*100)/100;
set(ax(1), 'YTick', Ytick_R);
set(ax(2), 'YTick', Ytick_PC);
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize',[10 7.5]);
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
%print(gcf,'-dpdf','../figures/fig_thr_est_time_tradeoff_AWGN_nu_upper');
%print(gcf,'-depsc','../figures/fig_thr_est_time_tradeoff_AWGN_nu_upper');
