%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: This m-file considers the Throughput-Estimation time 
%        tradeoff for the Rayliegh fading channel. The fading is fast        
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
%        Created on: 01.09.14
%        Last modified: 01.09.14
%        Revision History: 01.09.14 --> File generated   
%       
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
P_tran = -00;                             % Power transmitted by PR, the SNR received at ST can be 
                                          % varried using this parameter
noise_power = -100;                       % noise power -100 dBm
I_T = -110;                               % Interference temperature -80 dBm
f_s = 1e6;                                % 1 MHz one band
K = 0.1 * f_s;                            % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha_true_p = 100;                       % True Path loss between ST and PR  to achieve P_p = IT  
alpha_true_s = 080;                       % True Path loss between ST and SR  to achieve P_p = IT  
mu_P = 0.025;                             % Accuracy at PR
m_p = 1;                                  % Number of realizations of channel for the link ST-PR
m_s = 1;                                  % Number of realizations of channel for the link ST-SR 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sim
    M = 1e4;                                  % Number of realizations 
    % N = Total of samples used for estimation = tau * f_s 
    N_sim = round([0.001:0.01:0.09,0.1:0.1:0.9,1:1:10,20,30] * (1e-3 * f_s));     
    P_reg_sim = zeros(length(N_sim), M);      % Regulated transmitted power at ST
    P_p_sim = zeros(length(N_sim), M);        % Power received at the PR
    P_s_sim = zeros(length(N_sim), M);        % Power received at the SR
    R_sim = zeros(length(N_sim), M);          % Throughput at the SR
    P_rcvd_sim = zeros(length(N_sim), M);                 % Received Power at ST

    Exp_R_sim = zeros(1,length(N_sim));       % Expected throughput  
    Exp_P_reg_sim= zeros(1,length(N_sim));    % Expected power regulated by ST
    Exp_P_p_sim = zeros(1,length(N_sim));     % Expected power received at PR
    PC_P_p_sim = zeros(1,length(N_sim));      % Confidence probability power at PR  
    PC_R_sim = zeros(1,length(N_sim));        % Condidence probability rate at SR
   
    %% Theoretical parameters
    PC_P_p_th = zeros(1,length(N_sim));        % Confidence probability power at PR  
    nc_th = 1;    
    
    bins = 1e6;    
    
    for i=1:length(N_sim)        
        disp(strcat('N = ',num2str(i)));        
        g_p = ones(1, M) .* random('gam',m_p, 1, 1, M);
        g_s = random('gam',m_s, 1, 1, M); 
        for j=1:M
            % Estimate the received power P_reg with additive noise amplitude (zero mean and noise variance) 
            channel = g_p(j) * 10^(-alpha_true_p/10);
            samples = sqrt(10^(P_tran/10) * channel) * ones(1, N_sim(i)) + random('norm', 0, sqrt(10^(noise_power/10)), 1, N_sim(i));
            
            % Checkpoint 
            % Here the samples must be Gaussian distributed with mean = 10^(P_reg/10), 
            % examine by ploting the histogram, hist(samples, 100)
            P_rcvd_sim(i,j) = mean(samples.^2);    
        end       
            
        % Determine the scaling factor, which is the Expected Received
        % power at PR is scaled to meet the interference temperature
        %nc_sim = 1;     
        nc_sim = 1 / mean(g_p * 10^(-alpha_true_p/10) ./ P_rcvd_sim(i,:));             

        %temp = 1 / mean(g_p * 10^(-alpha_true_p/10) ./ P_rcvd_sim(i,:));
        disp(strcat('nc_sim = ',num2str(nc_sim))); 


        % Perform power control over the received samples, 
        % P_reg_sim should be Inverse gamma distributed, 
        % hist(P_reg_sim(i,:), 50)
        P_reg_sim(i,:) =  10^(I_T/10)  * nc_sim ./ P_rcvd_sim(i,:);           

        % Determine P_p, P_s and R from P_reg after imposing path loss
        % hist(P_p_sim(i,:), 50)
        P_p_sim(i,:) = g_p * 10^(-alpha_true_p/10) .* P_reg_sim(i,:);
        % hist(P_s_sim(i,:), 50)
        P_s_sim(i,:) = g_s * 10^(-alpha_true_s/10) .* P_reg_sim(i,:);
        % hist(R_sim(i,:), 50)
        R_sim(i,:) = (K- N_sim(i))/K * log2(1 + P_s_sim(i,:) / 10^(noise_power/10));  

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
        disp(strcat('PC_P_p_sim(i) = ',num2str(PC_P_p_sim(i)))); 
        if 0
            % Determine the PC_P_p theoretically         
            [f x_pts] = hist(P_p_sim(i,:), bins);    

            % To determine the nc_th theoretically the expected is computed
            % from the histogram of the receved power at PR, histogram of 
            % P_p_sim is cannot be used as it is already scaled with nc_sim
            if 0
                bins_nc = 1000;
                [f_nc x_pts_nc] = hist(g_p * 10^(-alpha_true_p/10) * 10^(I_T/10) ./ P_rcvd_sim(i,:), bins_nc);

                % nc_th is computed by dividing scaling the P_reg with the
                % interference temperature and the mean value of the power recieved
                % at the PR
                nc_th = 10^(I_T/10) / density_PR_fading(x_pts_nc,...
                N_sim(i), 10^(P_tran/10), 10^(I_T/10), 10^(-alpha_true_p/10), ...
                10^(noise_power/10), 1, x_pts_nc(2) - x_pts_nc(1), 'exp');
            else
                % or speeden up the simulation nc_sim can be used
                nc_th =  nc_sim;

            end
            disp(strcat('nc_th = ',num2str(nc_th))); 

            Exp_P_p_th = 10^(I_T/10); 
            index = find(x_pts <=  Exp_P_p_th * (1 + mu_P) & ...
                x_pts >=  Exp_P_p_th * (1 - mu_P));

            if isempty(index) 
                PC_P_p_th(i) = 0;
            else        
                PC_P_p_th(i) = density_PR_fading(x_pts(index),...
                N_sim(i), 10^(P_tran/10), 10^(I_T/10), 10^(-alpha_true_p/10), ...
                10^(noise_power/10), nc_th, x_pts(2) - x_pts(1), 'dis');       
            end
            disp(strcat('PC_P_p_th(i) = ',num2str(PC_P_p_th(i))));
        end
    end
%    save('results_thr_est_time_tradeoff_fading_wo_nu_snr_m10_sim.mat');
%    quit;
end   
