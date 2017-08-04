%                                
%        Description: This m-file considers the Estimation-Throughput tradeoff 
%        for the fading channel.  
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

sim = 0;                                      % Enable( = 1) to perform simulation, theoretical analysis
th = 0;                                       % disable ( = 0) to plot curves, data is read from a file
                                              % Theoretical anaylsis is also included of simulation as
                                              % numerical integration is involved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_tran = -00;                                  % Power transmitted by PR, the SNR received at ST can be 
                                              % varried using this parameter
noise_power = -100;                           % noise power -100 dBm
I_T = -110;                                   % Interference temperature -80 dBm
f_s = 1e6;                                    % 1 MHz one band
K = 0.1 * f_s;                                % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha_true_p = 100;                           % True Path loss between ST and PR  to achieve P_p = IT  
alpha_true_s = 080;                           % True Path loss between ST and SR  to achieve P_p = IT  
epsilon = 0.01;                               % Outage Probability constraint at PR
snr = 10^(P_tran/10) * 10^(-alpha_true_p/10)...  % snr received at the PR
    / 10^(noise_power/10);
m_p = 1;                                      % Number of realizations of channel for the link ST-PR
m_s = 1;                                      % Number of realizations of channel for the link ST-SR 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if sim
    M = 1e6;                                  % Number of realizations 
    if snr == 1
        N_sim = ceil(10.^[0.2:0.4:4.4]);          % N = Total of samples used for estimation = tau * f_s (ms) 
    elseif snr == 0.1
        N_sim = [100:2000:20000];            % N = Total of samples used for estimation = tau * f_s (ms)      
    end
    P_reg_sim = zeros(1, length(N_sim));      % Regulated transmitted power at ST
    P_rcvd_sim = zeros(1, M);                 % Received Power at ST

    Exp_R_sim = zeros(1,length(N_sim));       % Expected throughput  
    num_tp = 10000;                           % Number of test points
    for i=1:length(N_sim)
        disp(strcat('N = ',num2str(N_sim(i)))); 
        g_p = ones(1, M) .* random('gam',m_p, 1, 1, M);
        g_s = random('gam',m_s, 1, 1, M); 
        parfor j=1:M
            % Estimate the received power P_reg with additive noise amplitude (zero mean and noise variance) 
            samples = sqrt(g_p(j) * 10^((P_tran)/10) * 10^(-alpha_true_p/10)) * ones(1, N_sim(i)) + random('norm', 0, sqrt(10^(noise_power/10)), 1, N_sim(i));
            % Checkpoint 
            % Here the samples must be Gaussian distributed with mean = 10^(P_reg/10), 
            % examine by ploting the histogram, hist(samples, 100)
            P_rcvd_sim(j) = mean(samples.^2);                      
            
        end            
            % Perform power control over the received samples, 
            % hist(P_reg_sim(i,:), 200)
            % print(gcf,'-dpdf','../figures/hist_P_reg')
            test_points = linspace(10^(noise_power/10) * (1 + m_p * snr), max(P_rcvd_sim), num_tp);
            for k=1:num_tp
                %(length(find(P_rcvd_sim < test_points(k)))/M)
                if (length(find(P_rcvd_sim < test_points(k)))/M) >= (1 - epsilon)   
                    P_reg_sim(i) = 10^(I_T/10) * 10^(P_tran/10) /...
                        (test_points(k) - 10^(noise_power/10));       
                    break;
                end
            end            
       
          
            Exp_R_sim(i) = mean((K- N_sim(i))/K * log2(1 + g_s * 10^(-alpha_true_s/10) * P_reg_sim(i) / 10^(noise_power/10)));     
            
            disp(strcat('Exp_R_sim(i) = ',num2str(Exp_R_sim(i)))); 
            disp(strcat('P_reg_sim(i) = ',num2str(P_reg_sim(i)))); 
    end
    save('results_thr_est_time_tradeoff_fading_wo_nu_snr_m00_log_e01_sim.mat');
    if ~th      % If theoretical analysis is already simulated the computation ends here.
        quit;
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th
   %% Theoretical parameters    
    if snr == 1
        N_th = [1:1:39,40:30:2980,3000:500:10000,10500:0500:20000];  % N = Total of samples used for estimation = tau * f_s (ms)
    elseif snr == 0.1
        N_th = [100:200:20000];  % N = Total of samples used for estimation = tau * f_s (ms)
    end
   
   
    P_reg_th = zeros(1,length(N_th));          % Expected power regulated by ST
    Exp_R_th = zeros(1,length(N_th));          % Expected throughput  
    num_tp = 5000;                             % Number of test points
    P_reg_ga_th = zeros(1,length(N_th));       % Expected power regulated by ST (applying Gamma Approximation)
    Exp_R_ga_th = zeros(1,length(N_th));       % Expected throughput (applying Gamma Approximation) 

    for i=1:length(N_th)       
        disp(strcat('N_th = ',num2str(N_th(i))));             
           
       %% Determining the performance meterics --> Exp_R
       %% Expected values           
       
       %% Gamma Approximation to the non-central chi-squared distribution
       
       %% Determine the controlled power
       if 0 %% Method 1  -- Not Working
           temp_preg = 10^(I_T/10) * 10^(P_tran/10) ./...
               (10^(noise_power/10)/N_th(i) * (2 + 4 * snr)./(1 + snr) .*...
                gammaincinv((1 - epsilon), N_th(i) * (1 + snr).^2./(2 + 4 * snr)) - 10^(noise_power/10))
           func_preg = @(t) 10^(I_T/10) * 10^(P_tran/10) ./...
               (10^(noise_power/10)/N_th(i) * (2 + 4 * t * snr)./(1 + t * snr) .*...
                gammaincinv((1 - epsilon), N_th(i) * (1 + t * snr).^2./(2 + 4 * t * snr))...
                - 10^(noise_power/10)) .* exp(-t);       
           P_reg_th(i) = integral(func_preg, 0, 100);
       end
       
       if 1 %% Method 2  -- Working But takes to long       
           M = 1e3;
           g_p = ones(1, M) .* random('gam',m_p, 1, 1, M);
           P_rcvd_sim = zeros(1, M);
           parfor j=1:M
                % Estimate the received power P_reg with additive noise amplitude (zero mean and noise variance) 
                samples = sqrt(g_p(j) * 10^((P_tran)/10) * 10^(-alpha_true_p/10)) * ones(1, N_th(i)) + random('norm', 0, sqrt(10^(noise_power/10)), 1, N_th(i));
                % Checkpoint 
                % Here the samples must be Gaussian distributed with mean = 10^(P_reg/10), 
                % examine by ploting the histogram, hist(samples, 100)
                P_rcvd_sim(j) = mean(samples.^2);                      

            end 
            test_points = linspace(10^(noise_power/10) * (1 + m_p * snr), max(P_rcvd_sim), num_tp);
            for k = num_tp:-1:1                      
                func_prcvd = @(t) gammainc(test_points(k)./(10^(noise_power/10)/N_th(i)...
                    * (2 + 4 * t * snr)./(1 + t * snr)),...
                    N_th(i) * (1 + t * snr).^2./(2 + 4 * t * snr)) .* exp(-t);
                result = integral(func_prcvd, 0, 100);                
                if result <= (1 - epsilon)
                    test_points(i);
                    P_reg_th(i) = 10^(I_T/10) * 10^(P_tran/10)/...
                        (test_points(k) - 10^(noise_power/10));                    
                    break;
                end
            end
           %% Expected Rate
           func_r_em = @(u) (K- N_th(i))/K * log2(1 + u * P_reg_th(i) * 10^(-alpha_true_s/10) /...
                10^(noise_power/10)) .* exp(-u);
            Exp_R_th(i) = integral(func_r_em, 0, 100);
            disp(strcat('Exp_R_th(i) = ',num2str(Exp_R_th(i)))); 
            disp(strcat('P_reg_th(i) = ',num2str(P_reg_th(i))));
       end
       
       if 1  %% Method 3 -- Approximating the first two moments of the distribution for Prcvd with 
             %% the moments of Gamma distribution -- Not so Accurate for low SNR and for N > 1000

            %% The mean and variance are determined from the moment generating function 
            mean_nc = 10^(noise_power/10) * (1 + m_p * snr);
            var_nc =  (10^(noise_power/10))^2/N_th(i) * (2 * + 4 * m_p *...
                snr + N_th(i) * m_s * (snr)^2);
            % Parameters for Gamma distribution 
            b = var_nc / mean_nc;
            a = mean_nc / b;
            P_reg_ga_th(i) = 10^(I_T/10) * 10^(P_tran/10) / (b * gammaincinv((1 - epsilon),a) - 10^(noise_power/10));       
           
            %% Expected Rate
            func_r_em = @(u) (K- N_th(i))/K * log2(1 + u * P_reg_ga_th(i) * 10^(-alpha_true_s/10) /...
                10^(noise_power/10)) .* exp(-u);
            Exp_R_ga_th(i) = integral(func_r_em, 0, 100);
            disp(strcat('Exp_R_ga_th(i) = ',num2str(Exp_R_ga_th(i)))); 
            disp(strcat('Exp_R_ga_th(i) = ',num2str(P_reg_ga_th(i))));
       end 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Finding optimum estimation time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimum throughput for the ideal model
    [Exp_R_opt index] = max(Exp_R_th);
    N_opt = N_th(index);

    save('results_thr_est_time_tradeoff_fading_wo_nu_snr_m00_log_e01_th.mat');
    quit;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Epsilon = 0.01
load('results_thr_est_time_tradeoff_fading_wo_nu_snr_m00_log_e01_th.mat');
%pbaspect([1 0.3 1]);
[value indexOfInterest] = find(N_th);
Fontsize = 9;

h2 = semilogx(N_th(indexOfInterest) * 1e-3, Exp_R_th(indexOfInterest), 'k', 'LineWidth', 1);
Ymin = min(Exp_R_th);

%% Ideal Model
if 0 %% Average Constraint 
    func_r_id = @(t,u) log2(1 + t * 10^(I_T/10) * 10^(P_tran/10) * 10^(-alpha_true_s/10)...
        ./ (u  * 10^(-alpha_true_p/10) *  10^(noise_power/10))) .* exp(-t) .* exp(-u);
    Exp_R_id = integral2(func_r_id, 0, 100, 0,100);
end

if 1 %% Ouatge Constarint
    P_reg_id = 10^(I_T/10) / 10^(-alpha_true_p/10) / log(1/(epsilon));
   
    func_r_id = @(t) log2(1 + P_reg_id * t *  10^(-alpha_true_s/10)...
        ./ (10^(noise_power/10))) .* exp(-t);
    Exp_R_id = integral(func_r_id, 0, 100);
end
hold on,
h1 = plot(N_th(indexOfInterest) * 1e-3, Exp_R_id * ones(1, length(N_th(indexOfInterest))), 'c', 'LineWidth', 2);

%% Gamma Approximation
hold on,
h5 = plot(N_th(indexOfInterest) * 1e-3, Exp_R_ga_th(indexOfInterest), 'k--', 'LineWidth', 1);
hold on,
[Exp_R_ga_opt index] = max(Exp_R_ga_th);
N_ga_opt = N_th(index);
plot(N_ga_opt * 1e-3, Exp_R_ga_opt, 'rs', 'LineWidth', 1,'handlevisibility','off');

%% Optimum Value
hold on,
h3 = semilogx(N_opt * 1e-3, Exp_R_opt, 'rs', 'LineWidth', 1);


%% Epsilon = 0.10
load('results_thr_est_time_tradeoff_fading_wo_nu_snr_m00_log_e10_th.mat');
[value indexOfInterest] = find(N_th);
Fontsize = 9;
semilogx(N_th(indexOfInterest) * 1e-3, Exp_R_th(indexOfInterest), 'k', 'LineWidth', 1);

%% Ideal Model
if 0 %% Average Constraint 
    func_r_id = @(t,u) log2(1 + t * 10^(I_T/10) * 10^(P_tran/10) * 10^(-alpha_true_s/10)...
        ./ (u  * 10^(-alpha_true_p/10) *  10^(noise_power/10))) .* exp(-t) .* exp(-u);
    Exp_R_id = integral2(func_r_id, 0, 100, 0,100);
end

if 1 %% Ouatge Constarint
    P_reg_id = 10^(I_T/10) / 10^(-alpha_true_p/10) / log(1/(epsilon));
   
    func_r_id = @(t) log2(1 + P_reg_id * t *  10^(-alpha_true_s/10)...
        ./ (10^(noise_power/10))) .* exp(-t);
    Exp_R_id = integral(func_r_id, 0, 100);
end
hold on,
h1 = plot(N_th(indexOfInterest) * 1e-3, Exp_R_id * ones(1, length(N_th(indexOfInterest))), 'c', 'LineWidth', 2);

%% Gamma Approximation

hold on,
h5 = plot(N_th(indexOfInterest) * 1e-3, Exp_R_ga_th(indexOfInterest), 'k--', 'LineWidth', 1);
hold on,
[Exp_R_ga_opt index] = max(Exp_R_ga_th);
N_ga_opt = N_th(index);
plot(N_ga_opt * 1e-3, Exp_R_ga_opt, 'rs', 'LineWidth', 1,'handlevisibility','off');

Ymax = max(Exp_R_id);

hold on,
semilogx(N_opt * 1e-3, Exp_R_opt, 'rs', 'LineWidth', 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_thr_est_time_tradeoff_fading_wo_nu_snr_m00_log_e01_sim.mat');
hold on,
h4 = plot(N_sim * 1e-3, Exp_R_sim, 'ko', 'LineWidth', 1);

load('results_thr_est_time_tradeoff_fading_wo_nu_snr_m00_log_e10_sim.mat');
hold on,
plot(N_sim * 1e-3, Exp_R_sim, 'ko', 'LineWidth', 1);


grid on;
axis([0 max(N_th) * 1e-3 Ymin  Ymax * 1.03]);
ylabel('$\rs(\tau)$ [bits/sec/Hz]','FontSize',Fontsize);
xlabel('$\tau$ [ms]','FontSize',Fontsize);
hl = legend([h1 h2 h5 h3 h4],'IM', 'EM','Lemma 2', '$\trs(\ttau)$', 'Simulated');
%set(hl, 'position',[0.725 0.12 0.15 0.31]);
set(hl, 'position',[0.57 0.12 0.29 0.28]);
set(gca,'FontSize',Fontsize);

if 0
    %%%% Zommed Image %%%%%
    % create a new pair of axes inside current figure
    ax = axes('Position',[.30 .47 .28 .27]);
    box on; % put box around new pair of axes
    child_gca = get(gca,'Children');

    [value indexOfInterest] = find(N_th <= 600 & N_th > 30);
    plot(N_th(indexOfInterest) * 1e-3, Exp_R_th(indexOfInterest), 'k', 'LineWidth', 1); % plot on new axes
    hold on,
    plot(N_th(indexOfInterest) * 1e-3, Exp_R_ga_th(indexOfInterest), 'k--', 'LineWidth', 1); % plot on new axes
    hold on,
    plot(N_opt * 1e-3, Exp_R_opt, 'rs', 'LineWidth', 1);
    hold on,
    plot(N_ga_opt * 1e-3, Exp_R_ga_opt, 'rs', 'LineWidth', 1);
    ylabel('$\rs(\tau)$','FontSize',Fontsize);
    xlabel('$\tau$ [ms]','FontSize',Fontsize);
    set(ax, 'XTick', [0.2 0.4 0.6]);
    set(ax, 'YTick', [1.72 1.75 1.78]);
    axis(ax, [min(N_th(indexOfInterest))* 1e-3 max(N_th(indexOfInterest))* 1e-3...
        min(Exp_R_th(indexOfInterest))  max(Exp_R_ga_th(indexOfInterest)) * 1.001]);
end
laprint(1, '../figures/fig_thr_est_time_tradeoff_fading', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');
