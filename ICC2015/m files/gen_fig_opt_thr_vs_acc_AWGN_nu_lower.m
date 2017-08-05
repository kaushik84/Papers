%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: This m-file considers the Throughput-accuracy at PR 
%        anlysis for the AWGN channel.        
%        
%        The simulation is performed to show the following:
%        1) Analyse the throughput vs accuracy curves
%
%        Created on: 17.09.14
%        Last modified: 17.09.14
%        Revision History: 17.09.14 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

sim = 0;                                  % Enable( = 1) to perform simulation, theoretical analysis
th = 1;                                   % disable ( = 0) to plot curves, data is read from a file
                                          % Theoretical anaylsis is also included of simulation as
                                          % numerical integration is involved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu = -3;                                  % Noise uncertainty upper  
P_tran = 00;                              % Power transmitted by PR  
noise_power = -100;                       % Actual noise power -100 dBm
noise_power_p = -100 + nu;                % noise power at PR with nu noise power + nu dBm
noise_power_s = noise_power;              % noise power at PR without nu noise power dBm
I_T = -110;                               % Interference temperature -80 dBm
f_s = 1e6;                                % 1 MHz one band
K = 0.1 * f_s;                            % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha_true_p = 100;                       % True Path loss between ST and PR  to achieve P_p = IT  
alpha_true_s = 080;                       % True Path loss between ST and SR  to achieve P_p = IT  
epsilon_P = [0.95];                       % Probability constraint


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if sim
    M = 1e4;                                  % Number of realizations 
    mu_P_sim = 0.7:0.05:1;                    % Accuracy at PR
    N_sim = ceil([1:3:20] * (1e-3 * f_s));    % N = Total of samples used for estimation = tau * f_s 
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
        for j=1:M
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
    save('results_opt_thr_vs_acc_AWGN_nu_l_sim.mat');
    if ~th      % If theoretical analysis is already simulated the computation ends here.
        quit;
    end
end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Theoretical analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th
    %% Theoretical parameters
    mu_P_th = 0.001:0.002:0.1;                                   % Accuracy at PR
    Exp_R_opt = zeros(length(epsilon_P),length(mu_P_th));        % Expected power regulated by ST
    N_optimum = zeros(length(epsilon_P),length(mu_P_th));        % Optimum Sensing Time ST


    nc_th = (10^(noise_power_p/10)/10^(-alpha_true_p/10) + 10^(P_tran/10));
    for e  = 1:length(epsilon_P)
        disp(strcat('epsilon_P = ',num2str(epsilon_P(e)))); 
        parfor m = 1:length(mu_P_th)    

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Finding optimum Rate by determing the optimum sensing time
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            disp(strcat('mu_P_th = ',num2str(mu_P_th(m)))); 

            N_opt = ([0.01:0.01:100] * (1e-3 * f_s));

            for i = length(N_opt):-1:1           
                Exp_P_p_opt = 10^(I_T/10);

                % Confidence probability power at PR
                PC_P_p_opt = marcumq(sqrt(10^(-alpha_true_p/10) * 10^(P_tran/10) * ...
                    N_opt(i)/(10^(noise_power_p/10))),...
                    sqrt(10^(I_T/10) * nc_th * N_opt(i) * 10^(-alpha_true_p/10)/...
                    (10^(noise_power_p/10) * (1 + mu_P_th(m)) * Exp_P_p_opt)), ceil(N_opt(i)/2 - 1)) -... 
                marcumq(sqrt(10^(-alpha_true_p/10) * 10^(P_tran/10) * ...
                    N_opt(i)/(10^(noise_power_p/10))),...
                    sqrt(10^(I_T/10) * nc_th * N_opt(i) * 10^(-alpha_true_p/10)/...
                    (10^(noise_power_p/10) * (1 - mu_P_th(m)) * Exp_P_p_opt)), ceil(N_opt(i)/2 - 1));

                % Optimum PC        
                if (PC_P_p_opt < epsilon_P(e))
                    N_optimum(e, m) = N_opt(i);        
                    break;
                end
                %disp(strcat('PC_P_p_opt = ',num2str(PC_P_p_opt))); 
            end 
            Exp_R_opt(e, m) = (K- N_optimum(e, m))/K * log2(1 + 10^(I_T/10)/(10^(P_tran/10) *  10^(-alpha_true_p/10) +...
                    10^(noise_power_p/10)) * nc_th * 10^(-alpha_true_s/10)  / 10^(noise_power_s/10));
            disp(strcat('Exp_R_opt = ',num2str(Exp_R_opt(e, m))));
        end
    end
    save('results_opt_thr_vs_acc_AWGN_nu_l_th.mat');
    quit;
end
load('results_opt_thr_vs_acc_AWGN_nu_l_th.mat');

temp = 6:1:50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fontsize = 18;
plot(mu_P_th(temp), Exp_R_opt(1,(temp)), 'LineWidth',1.5, 'HandleVisibility','off');
hold on,
plot(mu_P_th(temp), Exp_R_opt(2,(temp)), 'LineWidth',1.5, 'HandleVisibility','off');
hold on,
plot(mu_P_th(temp), Exp_R_opt(3,(temp)), 'LineWidth',1.5, 'HandleVisibility','off');
if 0 
    load('results_opt_thr_vs_acc_AWGN_nu_l_sim.mat');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves simulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(ax,'NextPlot','add');
    h6 = plot(ax(1), N_sim * 1e-3, Exp_R_sim, 'r*', 'HandleVisibility','off');
    set(ax,'NextPlot','add');
    plot(ax(2), N_sim * 1e-3, PC_P_p_sim, 'r*', 'HandleVisibility','off');
end

grid on;
axis([min(mu_P_th(temp)) max(mu_P_th(temp)) min(min(Exp_R_opt(:,temp))) max(max(Exp_R_opt(:,10:50)))]);
ylabel('R_s = [bits/sec/Hz]','FontSize',Fontsize+2);
xlabel('\mu','FontSize',Fontsize+2);
set(gca,'FontSize',Fontsize);
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize',[10 7.5]);
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
%print(gcf,'-dpdf','../figures/fig_opt_thr_vs_acc_AWGN');
%print(gcf,'-depsc','../figures/fig_opt_thr_vs_acc_AWGN');