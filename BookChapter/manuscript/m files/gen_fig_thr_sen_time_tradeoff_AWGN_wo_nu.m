%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: This m-file considers the Throughput-Sesning time 
%        tradeoff with SNR estimation for the AWGN channel. In this m-script  
%        the short term analysis is investigted that is, considering invidual 
%        frame. Hence, the SNR is considered constant for the time interval 
%        (tau_est + tau_sen). 
%        This assumption however will become weak for lower SNR, because low
%        SNR regimes presents high tau_sen, leading to a low probability
%        that the channel remains constant.
%
%        Clearly the recieved power is distributed according to non central
%        chi-2 distribution, which can be approximated as Gaussian 
%        distribution for large sample count (This again needs an analysis, 
%        on how good is the approximation is?)
%        
%        The simulation is performed to show the following:
%        1) Analyse the throughput vs sensing for different cases 
%           - Ideal case, when SNR is known at the ST (i)
%           - Average Constarint
%           - Outage Constraint   
%        2) Confirm the theoretical analysis vs simulation results
%
%        Created on: 05.03.15
%        Last modified: 05.03.15
%        Revision History: 05.03.15 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
%warning off;


sim = 00;                                   % Enable( = 1) to perform simulation, theoretical analysis
th = 00;                                    % disable ( = 0) to plot curves, data is read from a file
                                            % Theoretical anaylsis is also included of simulation as
                                            % numerical integration is involved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_s = 10^(-10/10);                          % Power transmitted by ST, the SNR received at SR 
P_p = 10^(-10/10);                          % Power transmitted by PT, the SNR received at SR 
                                            % varried using this parameter
noise_power = 10^(-100/10);                 % noise power -100 dBm
f_s = 1e6;                                  % 1 MHz one band
K = 0.1 * f_s;                              % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha_p_1 = 10^(-100/10);                   % True Path loss between ST and PR   
alpha_p_2 = 10^(-100/10);                   % True Path loss between PR and SR   
alpha_s = 10^(-080/10);                     % True Path loss between ST and SR 
P_H0 = 0.8;                                 % Probability of Hypothesis 0
P_d_d = 0.90;                               % Constraint on Probability of detection P_d
mu = 0.05;                                  % Outage Probability on probability of detection

num_test_points = 10000;                    % Test points for evaluating threshold 

[C_0, C_1] = calc_capacities(alpha_p_2 * P_p / noise_power,... 
    alpha_s * P_s / noise_power,....
    noise_power);                           % Capacities at SR, with and without intereference from PT  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sim
     M = 1e4;                                           % Number of realizations 
     %% Simulations parameters
     N_est_sim = 5 * 1000;                              % N_est_th = number of samples used for estimation = tau * f_s 

     N_sen_sim = [1000:1500:12000];                     % N = Total of samples used for sensing = tau * f_s 

     epsilon_id_sim = zeros(1,length(N_sen_sim));       % Ideal threshold 
     epsilon_ac_sim = zeros(1,length(N_sen_sim));       % Accurate threshold 
     epsilon_oc_sim = zeros(1,length(N_sen_sim));       % Inaccurate threshold due to incorrect estimation of received power

     P_f_id_sim = zeros(1,length(N_sen_sim));           % Probability of false alarm at ST for the ideal case
     P_f_ac_sim = zeros(1,length(N_sen_sim));           % Probability of false alarm at ST due to accurate estimation of threshold
     P_f_oc_sim = zeros(1,length(N_sen_sim));           % Probability of false alarm at ST due to inaccurate estimation of threshold

     P_d_id_sim = zeros(1,length(N_sen_sim));           % Probability of detection at ST for the ideal case
     P_d_ac_sim = zeros(1,length(N_sen_sim));           % Probability of detection at ST due to accurate estimation of threshold
     P_d_oc_sim = zeros(1,length(N_sen_sim));           % Probability of detection at ST due to inaccurate estimation of threshold
    
     R_id_sim = zeros(1,length(N_sen_sim));             % Throughput at SR ideal case
     R_ac_sim = zeros(1,length(N_sen_sim));             % Throughput at SR with accurate estimation
     R_oc_sim = zeros(1,length(N_sen_sim));             % Throughput at SR with upper limit of the energy wall

     P_rcvd_H1_sen_sim = zeros(1,M);                    % Rcvd power samples under Hypothesis 1 (sensing) 
     P_rcvd_H0_sen_sim = zeros(1,M);                    % Rcvd power samples under Hypothesis 0 (sensing)
     P_rcvd_H1_est_sim = zeros(1,M);                    % Rcvd power samples under Hypothesis 1 (estimation) 
     %P_rcvd_H0_est_sim = zeros(1,M);                    % Rcvd power samples under Hypothesis 0 (estimation)
     
          
     
     %% Computation of the estimation 
     % Accurarte energy  
     %Acc_energy = (noise_power + alpha_p_1 * P_p);
     
     for i=1:length(N_sen_sim)
        disp(strcat('N = ',num2str(N_sen_sim(i)))); 
        
        %% Ideal Case -- Remember (no estimation)    
        parfor j=1:M
            P_rcvd_H1_sen_sim(j) = mean(random('norm',...
                random('norm',0, sqrt(P_p * alpha_p_1), 1, N_sen_sim(i)),...
                sqrt(noise_power), 1, N_sen_sim(i)).^2);
        end 
        parfor j=1:M
            P_rcvd_H0_sen_sim(j) = mean(random('norm', 0, sqrt(noise_power), 1, N_sen_sim(i)).^2);
        end 


        %% Determining the performance meteric R for all the different cases
        %% Threshold     
        test_points = linspace(min(P_rcvd_H1_sen_sim), max(P_rcvd_H1_sen_sim),...
            num_test_points);
        for j = length(test_points):-1:1
            if (length(find(P_rcvd_H1_sen_sim > test_points(j)))/M) >= P_d_d   
                epsilon_id_sim(i) = test_points(j);
                break;
            end
        end           

        %% Probability of false alarm and probability of detection
        P_f_id_sim(i) = length(find(P_rcvd_H0_sen_sim > epsilon_id_sim(i)))/M;
        P_d_id_sim(i) = length(find(P_rcvd_H1_sen_sim > epsilon_id_sim(i)))/M;

        %% Expected Rate
        C_0_sim = log2(1 + alpha_s * P_s / noise_power); 
        C_1_sim = log2(1 + alpha_s * P_s / (P_p * alpha_p_2 +...
            noise_power) );
        
        R_id_sim(i) = (K - N_sen_sim(i))/K * (P_H0 * (1 -  P_f_id_sim(i)) *...
            C_0_sim +  (1 - P_H0) * (1 - P_d_d) * C_1_sim);
        
        %% non ideal case -- includes estimation followed by sesning
        if N_sen_sim(i) > N_est_sim
            N_sen_ni_sim = N_sen_sim(i);% - N_est_sim;        % Sensing time for non ideal case

            %% Estimation Phase 
            parfor j=1:M
                P_rcvd_H1_est_sim(j) = mean(random('norm',...
                random('norm',0, sqrt(P_p * alpha_p_1), 1, N_est_sim),...
                sqrt(noise_power), 1, N_est_sim).^2);                               
            end 
            
            %% Sensing Phase
            parfor j=1:M
                P_rcvd_H1_sen_sim(j) = mean(random('norm',...
                    random('norm',0, sqrt(P_p * alpha_p_1), 1, N_sen_ni_sim),...
                    sqrt(noise_power), 1, N_sen_ni_sim).^2);
            end 
            parfor j=1:M
                P_rcvd_H0_sen_sim(j) = mean(random('norm', 0, sqrt(noise_power), 1, N_sen_ni_sim).^2);
            end 
            
            %% Generate test points to determine the threshold
            test_points = linspace(2 * epsilon_id_sim(i), 0.1 * epsilon_id_sim(i), num_test_points);

            %% Average constraint
            if 1
               %% Threshold  
                for j = 1:length(test_points)
                    P_d_test_points = gammainc( N_sen_ni_sim * test_points(j)./...
                        (2 * P_rcvd_H1_est_sim), N_sen_ni_sim/2, 'upper');
                    if  mean(P_d_test_points) >= P_d_d   
                        epsilon_ac_sim(i) = test_points(j);
                        %% Expected probability of detection is calculated
                        %% as P_d is random variable;
                        P_d_ac_sim(i) = mean(P_d_test_points);
                        break;
                    end
                end
                                  
                %% Probability of false alarm and probability of detection
                P_f_ac_sim(i) = length(find(P_rcvd_H0_sen_sim > epsilon_ac_sim(i)))/M;

                R_ac_sim(i) = (K - N_sen_sim(i))/K * (P_H0 * (1 -  P_f_ac_sim(i)) *...
                    C_0 +  (1 - P_H0) * (1 - P_d_ac_sim(i)) * C_1);   
            end
            
            %% Outage constraint            
            if 1
                %% Threshold                  
                for j = 1:length(test_points)
                    P_d_test_points = gammainc(N_sen_ni_sim * test_points(j)./...
                        (2 * P_rcvd_H1_est_sim), N_sen_ni_sim/2, 'upper');
                    if (length(find(P_d_test_points <= P_d_d))/M <= mu)    
                        epsilon_oc_sim(i) = test_points(j);
                        %% Expected probability of detection is calculated
                        %% as P_d is random variable;
                        P_d_oc_sim(i) = mean(P_d_test_points);
                        break;
                    end
                end
                   
                %% Probability of false alarm and probability of detection
                P_f_oc_sim(i) = length(find(P_rcvd_H0_sen_sim > epsilon_oc_sim(i)))/M;

                R_oc_sim(i) = (K - N_sen_sim(i))/K * (P_H0 * (1 -  P_f_oc_sim(i)) *...
                    C_0 +  (1 - P_H0) * (1 - P_d_oc_sim(i)) * C_1);   
            end
        end
            
        disp(strcat('epsilon_id_sim(i) = ',num2str(epsilon_id_sim(i))));   
        disp(strcat('epsilon_ac_sim(i) = ',num2str(epsilon_ac_sim(i))));   
        disp(strcat('epsilon_oc_sim(i) = ',num2str(epsilon_oc_sim(i))));   
        disp(strcat('P_f_id_sim(i) = ',num2str(P_f_id_sim(i)))); 
        disp(strcat('P_f_ac_sim(i) = ',num2str(P_f_ac_sim(i)))); 
        disp(strcat('P_f_oc_sim(i) = ',num2str(P_f_oc_sim(i))));
        disp(strcat('P_d_id_sim(i) = ',num2str(P_d_id_sim(i)))); 
        disp(strcat('P_d_ac_sim(i) = ',num2str(P_d_ac_sim(i)))); 
        disp(strcat('P_d_oc_sim(i) = ',num2str(P_d_oc_sim(i)))); 
        disp(strcat('R_id_sim(i) = ',num2str(R_id_sim(i)))); 
        disp(strcat('R_ac_sim(i) = ',num2str(R_ac_sim(i)))); 
        disp(strcat('R_oc_sim(i) = ',num2str(R_oc_sim(i)))); 
    end
    save('results_thr_sen_time_tradeoff_AWGN_wo_nu_snr_m10_sim2.mat');
    if ~th      % If theoretical analysis is already simulated the computation ends here.
        quit;
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th
   %% Theoretical parameters
    N_est_th = 5 * 1000;                              % N_est_th = number of samples used for estimation = tau * f_s 
    N_sen_th = [100:100:4900,4900:1:5500,5500:100:20000]; % N = includes the time in samples, for the ideal case their is no
                                                      % estimation time hence sensing time = N_th and the throughput will 
                                                      % attain a nonzero value. For the accurate, the energy walls sensing 
                                                      % time begin after the estimation of the rcvd energy is performed.
                                                      % N_sen = N_th - N_est
%    epsilon_id_th = zeros(1,length(N_sen_th));        % Ideal threshold 
%    epsilon_ac_th = zeros(1,length(N_sen_th));        % AC threshold 
%    epsilon_oc_th = zeros(1,length(N_sen_th));        % OC threshold

    P_f_id_th = zeros(1,length(N_sen_th));            % Probability of false alarm at ST for the ideal case
    P_f_ac_th = zeros(1,length(N_sen_th));            % Probability of false alarm at ST due to accurate estimation of threshold
    P_f_oc_th = zeros(1,length(N_sen_th));            % Probability of false alarm at ST due to inaccurate estimation of threshold

    P_d_id_th = zeros(1,length(N_sen_th));            % Probability of detection at ST for the ideal case
    P_d_ac_th = zeros(1,length(N_sen_th));            % Probability of detection at ST due to accurate estimation of threshold
    P_d_oc_th = zeros(1,length(N_sen_th));            % Probability of detection at ST due to inaccurate estimation of threshold
    
    
    C_0_th = 0;                                       % Capacity under scenario 0 at SR
    C_1_th = 0;                                       % Capacity under scenario 1 at SR
    R_id_th = zeros(1,length(N_sen_th));              % Throughput at SR ideal case
    R_ac_th = zeros(1,length(N_sen_th));              % Throughput at SR with accurate estimation
    R_oc_th = zeros(1,length(N_sen_th));              % Throughput at SR with upper limit of the energy wall
    
    
    % Accurarte energy  
    Acc_energy = (noise_power + alpha_p_1 * P_p);      
 
    parfor i=1:length(N_sen_th)
        disp(strcat('N = ',num2str(N_sen_th(i))));     
        epsilon_id_th = 0;
        epsilon_ac_th = 0;
        epsilon_oc_th = 0;

        %% Determining the performance meteric R for all the different cases
        %% Ideal         
        %% Threshold  
        epsilon_id_th = 2 * Acc_energy/N_sen_th(i) *...
            gammaincinv(P_d_d, N_sen_th(i)/2, 'upper');                  
       
        %% Probability of false alarm and probability of dectection
        P_f_id_th(i) = gammainc(N_sen_th(i) * epsilon_id_th/(2 * noise_power), N_sen_th(i)/2, 'upper');
        P_d_id_th(i) =  P_d_d;

       
        %% Expected Rate
        C_0_th = log2(1 + alpha_s * P_s / noise_power); 
        C_1_th = log2(1 + alpha_s * P_s / (P_p * alpha_p_2 +...
            noise_power) );
        
        R_id_th(i) = (K - N_sen_th(i))/K * (P_H0 * (1 -  P_f_id_th(i)) * C_0_th +...
            (1 - P_H0) * (1 - P_d_id_th(i)) * C_1_th);  
        

        if N_sen_th(i) > N_est_th
           N_sen_ni_th = N_sen_th(i);% - N_est_th;
           %% Adjacent case
           if 1
               %% Threshold  
               test_points = linspace(epsilon_id_th, 0.1 * epsilon_id_th, num_test_points);
               Exp_P_d_ac = zeros(1,num_test_points);
               for j=1:num_test_points
                    func_exp_P_d = @(t) gammainc( N_sen_ni_th/2 * test_points(j) ./ t, N_sen_ni_th/2, 'upper') *...
                        N_est_th/Acc_energy .* exp(-N_est_th/2 * log(2) - gammaln(N_est_th/2)  +...
                        (N_est_th/2 - 1) * log(N_est_th * t/Acc_energy) +...
                        (-t * N_est_th/(2 * Acc_energy)));                    
                    Exp_P_d_ac(j) =  integral(func_exp_P_d, 0, test_points(j) * 100);
                    if Exp_P_d_ac(j) >= P_d_d
                       epsilon_ac_th = test_points(j);
                       %% Expected P_d 
                       P_d_ac_th(i) = Exp_P_d_ac(j);
                       break;
                    else
                       epsilon_ac_th = test_points(j); 
                    end
               end

               %% Probability of false alarm and probability of dectection
               P_f_ac_th(i) = gammainc( N_sen_ni_th * epsilon_ac_th/(2 * noise_power), N_sen_ni_th/2, 'upper');  

               R_ac_th(i) = (K - N_sen_th(i))/K * (P_H0 * (1 -  P_f_ac_th(i)) * C_0 +...
                   (1 - P_H0) * (1 - P_d_ac_th(i)) * C_1); 
           end

 
            %% Outage case
            if 1
                %% Threshold  
                epsilon_oc_th = 4 * Acc_energy * gammaincinv(1 - mu, N_est_th/2, 'upper') * ...
                    gammaincinv(P_d_d, N_sen_ni_th/2, 'upper') / (N_est_th * N_sen_ni_th);            

                %% Probability of false alarm and expected probability of dectection
                P_f_oc_th(i) = gammainc(N_sen_ni_th * epsilon_oc_th/(2 * noise_power), N_sen_ni_th/2, 'upper');
                
                func_exp_P_d = @(t) gammainc(N_sen_ni_th/2, N_sen_ni_th/2 * epsilon_oc_th ./ t) *...
                        N_est_th/Acc_energy .* exp(-N_est_th/2 * log(2) - gammaln(N_est_th/2)  +...
                        (N_est_th/2 - 1) * log(N_est_th * t/Acc_energy) +...
                        (-t * N_est_th/(2 * Acc_energy)));                    
                P_d_oc_th(i) =  integral(func_exp_P_d, 0, epsilon_oc_th * 100);
                
                R_oc_th(i) = (K - N_sen_th(i))/K * (P_H0 * (1 -  P_f_oc_th(i)) * C_0 +...
                (1 - P_H0) * (1 - P_d_oc_th(i)) * C_1);
            end
        end        
        
        disp(strcat('epsilon_id_th = ',num2str(epsilon_id_th)));
        disp(strcat('epsilon_ac_th = ',num2str(epsilon_ac_th)));
        disp(strcat('epsilon_oc_th = ',num2str(epsilon_oc_th)));
        disp(strcat('P_f_id_th(i) = ',num2str(P_f_id_th(i)))); 
        disp(strcat('P_f_ac_th(i) = ',num2str(P_f_ac_th(i)))); 
        disp(strcat('P_f_oc_th(i) = ',num2str(P_f_oc_th(i))));
        disp(strcat('P_d_id_th(i) = ',num2str(P_d_id_th(i)))); 
        disp(strcat('P_d_ac_th(i) = ',num2str(P_d_ac_th(i)))); 
        disp(strcat('P_d_oc_th(i) = ',num2str(P_d_oc_th(i)))); 
        disp(strcat('R_id_th(i) = ',num2str(R_id_th(i)))); 
        disp(strcat('R_ac_th(i) = ',num2str(R_ac_th(i)))); 
        disp(strcat('R_oc_th(i) = ',num2str(R_oc_th(i)))); 


    end
    % Optimum throughput for the ideal model
    [R_id_opt index] = max(R_id_th);
    N_id_opt = N_sen_th(index);
    
    % Optimum throughput for the estimation model -- average constraint case
    [R_ac_opt index] = max(R_ac_th);
    N_ac_opt = N_sen_th(index);
    
    % Optimum throughput for the estimation model -- outage constraint case
    [R_oc_opt index] = max(R_oc_th);
    N_oc_opt = N_sen_th(index);
    
    save('results_thr_sen_time_tradeoff_AWGN_wo_nu_snr_m10_th2.mat');
    quit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Sensing throughput tradeoff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_thr_sen_time_tradeoff_AWGN_wo_nu_snr_m10_th2.mat');
Fontsize = 8;
[values index] = find(N_sen_th <= 10000);
grey=[0.7,0.7,0.7];
h1 = plot(N_sen_th(index) * 1e-3, R_id_th(index),'Color', grey, 'LineWidth',2);
hold on,
h2 = plot(N_sen_th(index) * 1e-3, R_ac_th(index), 'k--', 'LineWidth', 1);
hold on,
h3 = plot(N_sen_th(index) * 1e-3, R_oc_th(index), 'k-.', 'LineWidth', 1);

if 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves simulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('results_thr_sen_time_tradeoff_AWGN_wo_nu_snr_m10_sim2.mat');
    if 1
        C_0_sim = log2(1 + alpha_s * P_s / noise_power); 
        C_1_sim = log2(1 + alpha_s * P_s / (P_p * alpha_p_2 +...
            noise_power) );
        
        R_id_sim = (K - N_sen_sim)/K .* (P_H0 * (1 -  P_f_id_sim) *...
            C_0_sim +  (1 - P_H0) * (1 - P_d_d) * C_1_sim);
    end
    hold on,
    h4 = plot(N_sen_sim * 1e-3, R_id_sim, 'rx', 'LineWidth',1);
    hold on,
    plot(N_sen_sim * 1e-3, R_ac_sim, 'rx', 'LineWidth',1, 'HandleVisibility','off');
    hold on,
    plot(N_sen_sim * 1e-3, R_oc_sim, 'rx', 'LineWidth',1, 'HandleVisibility','off');
end 
hold on,
h5 = plot(N_id_opt * 1e-3, R_id_opt, 'rs', 'LineWidth',1);
hold on,
plot(N_ac_opt * 1e-3, R_ac_opt, 'rs', 'LineWidth',1, 'HandleVisibility','off');
hold on,
plot(N_oc_opt * 1e-3, R_oc_opt, 'rs', 'LineWidth',1, 'HandleVisibility','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves cosmetic makeover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid on;
axis([0.1 max(N_sen_th(index)) * 1e-3 0  max(R_id_opt) * 1.03]);
ylabel('$\rs(\test = \SI{5}{ms}, \tsen)$ [bits/sec/Hz]','FontSize',Fontsize);
xlabel('$\tsen$ [ms]','FontSize',Fontsize);
%hl = legend([h1 h2 h3 h5 h4],'$\rs$', '$\rsac$', '$\rsoc$', '$\ttsen$', 'Sim');
hl = legend([h1 h2 h3 h5 h4],'IM', 'EM-AC, Problem 1', 'EM-OC, Problem 2', '$\trs(\test,\ttsen)$', 'Simulated');
%set(hl, 'position',[0.725 0.12 0.15 0.31]);
set(hl, 'position',[0.16 0.2 0.32 0.31]);
set(gca,'FontSize',Fontsize);


if 1
    %%%% Zommed Image %%%%%
    % create a new pair of axes inside current figure
    ax = axes('Position',[.595 .26 .31 .31]);
    box on; % put box around new pair of axes
    child_gca = get(gca,'Children');

    [value indexOfInterest] = find(N_sen_th <= 5600 & N_sen_th >= 4990);
    plot(N_sen_th(indexOfInterest) * 1e-3, R_id_th(indexOfInterest), 'Color', grey, 'LineWidth',2);
    hold on,
    plot(N_sen_th(indexOfInterest) * 1e-3, R_ac_th(indexOfInterest), 'k--', 'LineWidth', 1);
    hold on,
    plot(N_sen_th(indexOfInterest) * 1e-3, R_oc_th(indexOfInterest), 'k-.', 'LineWidth', 1);
    
    plot(N_id_opt * 1e-3, R_id_opt, 'rs', 'LineWidth',1);  
    hold on,
    plot(N_ac_opt * 1e-3, R_ac_opt, 'rs', 'LineWidth',1, 'HandleVisibility','off');
    hold on,
    plot(N_oc_opt * 1e-3, R_oc_opt, 'rs', 'LineWidth',1, 'HandleVisibility','off');
    
    index = 4;
    hold on,
    h4 = plot(N_sen_sim(index) * 1e-3, R_id_sim(index), 'rx', 'LineWidth',1);
    hold on,
    plot(N_sen_sim(index) * 1e-3, R_ac_sim(index), 'rx', 'LineWidth',1, 'HandleVisibility','off');
    hold on,
    plot(N_sen_sim(index) * 1e-3, R_oc_sim(index), 'rx', 'LineWidth',1, 'HandleVisibility','off');
    
    % plot on new axes
    %ylabel(gca, '$\rs$','FontSize',Fontsize);
    %xlabel(gca, '$\test + \tsen$','FontSize',Fontsize);
    %set(ax, 'XTick', [-1.3 -0.9 -0.5]);
    %set(ax, 'YTick', [3.2 3.4 3.7]);
    axis(ax, [min(N_sen_th(indexOfInterest)) * 1e-3 max(N_sen_th(indexOfInterest)) * 1e-3...
        2.54 max(R_id_th(indexOfInterest)) * 1]);
    set(gca,'FontSize',Fontsize);
    grid(ax, 'on');
    handle = title(ax, 'Zoom');
end

laprint(1, '../figures/fig_thr_sen_time_tradeoff_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');
