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


sim = 00;                                  % Enable( = 1) to perform simulation, theoretical analysis
th = 01;                                  % disable ( = 0) to plot curves, data is read from a file
                                          % Theoretical anaylsis is also included of simulation as
                                          % numerical integration is involved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_s = 10^(-10/10);                        % Power transmitted by PR, the SNR received at ST 
P_p = 10^(-10/10);                        % Power transmitted by SR, the SNR received at SR 
                                          % varried using this parameter
noise_power = 10^(-100/10);               % noise power -100 dBm
f_s = 1e6;                                % 1 MHz one band
K = 0.1 * f_s;                            % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha_p_1 = 10^(-100/10);                 % True Path loss between ST and PR   
alpha_p_2 = 10^(-100/10);                 % True Path loss between PR and SR   
alpha_s = 10^(-080/10);                   % True Path loss between ST and SR 
P_H0 = 0.8;                               % Probability of Hypothesis 0
P_d_d = 0.90;                             % Constraint on Probability of detection P_d
mu = 0.05;                                % Outage Probability on probability of detection

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sim
     M = 1e4;                                       % Number of realizations 
     %% Simulations parameters
     N_est_sim = 5 * 1000;                           % N_est_th = number of samples used for estimation = tau * f_s 

     N_sim = [1:2.5:20] * 1000;                       % N = Total of samples used for estimation = tau * f_s 

     epsilon_id_sim = zeros(1,length(N_sim));       % Ideal threshold 
     epsilon_ac_sim = zeros(1,length(N_sim));       % Accurate threshold 
     epsilon_oc_sim = zeros(1,length(N_sim));       % Inaccurate threshold due to incorrect estimation of received power

     P_f_id_sim = zeros(1,length(N_sim));           % Probability of false alarm at ST for the ideal case
     P_f_ac_sim = zeros(1,length(N_sim));           % Probability of false alarm at ST due to accurate estimation of threshold
     P_f_oc_sim = zeros(1,length(N_sim));           % Probability of false alarm at ST due to inaccurate estimation of threshold

     P_d_id_sim = zeros(1,length(N_sim));           % Probability of detection at ST for the ideal case
     P_d_ac_sim = zeros(1,length(N_sim));           % Probability of detection at ST due to accurate estimation of threshold
     P_d_oc_sim = zeros(1,length(N_sim));           % Probability of detection at ST due to inaccurate estimation of threshold
    
     C_0_sim = 0;                                   % Capacity under scenario 0 at SR
     C_1_sim = 0;                                   % Capacity under scenario 1 at SR
     R_id_sim = zeros(1,length(N_sim));             % Throughput at SR ideal case
     R_ac_sim = zeros(1,length(N_sim));             % Throughput at SR with accurate estimation
     R_oc_sim = zeros(1,length(N_sim));             % Throughput at SR with upper limit of the energy wall

     P_rcvd_H1_sen_sim = zeros(1,M);                % Rcvd power samples under Hypothesis 1 (sensing) 
     P_rcvd_H0_sen_sim = zeros(1,M);                % Rcvd power samples under Hypothesis 0 (sensing)
     P_rcvd_H1_est_sim = zeros(1,M);                % Rcvd power samples under Hypothesis 1 (estimation) 
     P_rcvd_H0_est_sim = zeros(1,M);                % Rcvd power samples under Hypothesis 0 (estimation)
     
     num_test_points = 100000;                      % Test point for evaluating P_d, P_f
     
     
     
     %% Computation of the estimation 
     % Accurarte energy  
     %Acc_energy = (noise_power + alpha_p_1 * P_p);
     
     for i=1:length(N_sim)
        disp(strcat('N = ',num2str(N_sim(i)))); 
        
        %% Ideal Case -- (no estimation)    
        parfor j=1:M
            % Estimate the received power P_reg with additive noise amplitude (zero mean and noise variance) 
            samples = sqrt(P_p * alpha_p_1) * ones(1, N_sim(i)) + random('norm', 0, sqrt(noise_power), 1, N_sim(i));
            P_rcvd_H1_sen_sim(j) = mean(samples.^2);                                  
        end 
        parfor j=1:M
            % Estimate the received power P_reg without additive noise amplitude (zero mean and noise variance) 
            samples = random('norm', 0, sqrt(noise_power), 1, N_sim(i));
            P_rcvd_H0_sen_sim(j) = mean(samples.^2);                                  
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

        R_id_sim(i) = (K - N_sim(i))/K * (P_H0 * (1 -  P_f_id_sim(i)) *...
            C_0_sim +  (1 - P_H0) * (1 - P_d_d) * C_1_sim);
        
        %% non ideal case -- includes estimation followed by sesning
        if N_sim(i) > N_est_sim
            N_sen_sim = N_sim(i) - N_est_sim;

            %% Estimation
            parfor j=1:M
                % Estimate the received power P_reg with additive noise amplitude (zero mean and noise variance) 
                samples = sqrt(P_p * alpha_p_1) * ones(1, N_est_sim) + random('norm', 0, sqrt(noise_power), 1, N_est_sim);
                P_rcvd_H1_est_sim(j) = mean(samples.^2);                                  
            end 
            
            %% Sensing
            if 0 %--- not used for the simulation (P_d is evaluated using qfunc)
                for j=1:M
                    samples = sqrt(P_p * alpha_p_1) * ones(1, N_sen_sim) + random('norm', 0, sqrt(noise_power), 1, N_est_sim);
                    P_rcvd_H1_sen_sim(j) = mean(samples.^2);                                  
                end
            end
            parfor j=1:M
                samples = random('norm', 0, sqrt(noise_power), 1, N_sen_sim);
                P_rcvd_H0_sen_sim(j) = mean(samples.^2);                                  
            end
            
            
            Acc_energy = (noise_power + alpha_p_1 * P_p);
            test_points = linspace(Acc_energy* (1 + sqrt(2/N_sen_sim)* qfuncinv(0.999999)),...
                epsilon_id_sim(i), num_test_points);
            %test_points = linspace(min(P_rcvd_H1_sim), max(P_rcvd_H1_sim),...
                %num_test_points);
            %% Average constraint
            if 1
               %% Threshold  
               %% TO DO --- can be replaced with the non-central chi sqaured distribution (MarcumQ)
                %test_points = linspace(min(P_rcvd_H1_sim), max(P_rcvd_H1_sim),...
                %num_test_points);
                for j = length(test_points):-1:1
                    P_d_test_points = qfunc((test_points(j) - P_rcvd_H1_est_sim)./(sqrt(2/N_sen_sim) *...
                        P_rcvd_H1_est_sim));
                    if  mean(P_d_test_points) >= P_d_d   
                        epsilon_ac_sim(i) = test_points(j);
                        break;
                    end
                end
                
                %% test
                if 0
                    bins = 100;
                    ts = sqrt(2/N_sen_sim);
                    te = sqrt(2/N_est_sim);
                    Acc_energy = (noise_power + alpha_p_1 * P_p);   
                    [f x_pts] = hist(P_d_test_points, bins); 
                    pdf_pd_ana = epsilon_ac_sim(i) .*...
                        exp( (qfuncinv(x_pts)/sqrt(2)).^2 - (epsilon_ac_sim(i)./(1 + ts * qfuncinv(x_pts)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
                        * ts ./ (Acc_energy * te * ( 1 + ts * qfuncinv(x_pts)).^2);
                    
                    % Plotting curves
                    bar(x_pts,f/sum(f));    
                    hold on,
                    plot(x_pts, pdf_pd_ana * ((x_pts(3) - x_pts(2))), 'r', 'Linewidth',2.5);                   
                    
                end
                   
                %% Probability of false alarm and probability of detection
                P_f_ac_sim(i) = length(find(P_rcvd_H0_sen_sim > epsilon_ac_sim(i)))/M;
                P_d_ac_sim(i) = mean(qfunc((epsilon_ac_sim(i) - P_rcvd_H1_est_sim)./(sqrt(2/N_sen_sim) *...
                        P_rcvd_H1_est_sim)));

                %% Expected Rate
                C_0_sim = log2(1 + alpha_s * P_s / noise_power); 
                C_1_sim = log2(1 + alpha_s * P_s / (P_p * alpha_p_2 +...
                    noise_power) );

                R_ac_sim(i) = (K - N_sim(i))/K * (P_H0 * (1 -  P_f_ac_sim(i)) *...
                    C_0_sim +  (1 - P_H0) * (1 - P_d_ac_sim(i)) * C_1_sim);   
            end
            
            %% Outage constraint            
            if 1
                %% Threshold  
                %% TO DO --- can be replaced with the non-central chi sqaured distribution (MarcumQ)
                
                %linspace(min(P_rcvd_H1_sim), max(P_rcvd_H1_sim),...
                %num_test_points);
                for j = length(test_points):-1:1
                    P_d_test_points = qfunc((test_points(j) - P_rcvd_H1_est_sim)./(sqrt(2/N_sen_sim) *...
                        P_rcvd_H1_est_sim));
                    if (length(find(P_d_test_points <= P_d_d))/M <= mu)    
                        epsilon_oc_sim(i) = test_points(j);
                        break;
                    end
                end
                   
                %% Probability of false alarm and probability of detection
                P_f_oc_sim(i) = length(find(P_rcvd_H0_sen_sim > epsilon_oc_sim(i)))/M;
                P_d_oc_sim(i) = mean(qfunc((epsilon_oc_sim(i) - P_rcvd_H1_est_sim)./(sqrt(2/N_sen_sim) *...
                        P_rcvd_H1_est_sim)));

                %% Expected Rate
                C_0_sim = log2(1 + alpha_s * P_s / noise_power); 
                C_1_sim = log2(1 + alpha_s * P_s / (P_p * alpha_p_2 +...
                    noise_power) );

                R_oc_sim(i) = (K - N_sim(i))/K * (P_H0 * (1 -  P_f_oc_sim(i)) *...
                    C_0_sim +  (1 - P_H0) * (1 - P_d_oc_sim(i)) * C_1_sim);   
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
    save('results_thr_sen_time_tradeoff_AWGN_wo_nu_snr_m10_sim.mat');
    if ~th      % If theoretical analysis is already simulated the computation ends here.
        %quit;
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th
   %% Theoretical parameters
    N_est_th = 1700;%5 * 1000;                               % N_est_th = number of samples used for estimation = tau * f_s 
    N_th = [0.5:0.5:19] * 1000;%[0.1:0.1:35] * 1000;                        % N = includes the time in samples, for the ideal case their is no
                                                       % estimation time hence sensing time = N_th and the throughput will 
                                                       % attain a nonzero value. For the accurate, the energy walls sensing 
                                                       % time begin after the estimation of the rcvd energy is performed.
                                                       % N_sen = N_th - N_est
    epsilon_id_th = zeros(1,length(N_th));             % Ideal threshold 
    epsilon_ac_th = zeros(1,length(N_th));             % AC threshold 
    epsilon_oc_th = zeros(1,length(N_th));             % OC threshold

    P_f_id_th = zeros(1,length(N_th));                % Probability of false alarm at ST for the ideal case
    P_f_ac_th = zeros(1,length(N_th));                % Probability of false alarm at ST due to accurate estimation of threshold
    P_f_oc_th = zeros(1,length(N_th));                % Probability of false alarm at ST due to inaccurate estimation of threshold

    P_d_id_th = zeros(1,length(N_th));                % Probability of detection at ST for the ideal case
    P_d_ac_th = zeros(1,length(N_th));                % Probability of detection at ST due to accurate estimation of threshold
    P_d_oc_th = zeros(1,length(N_th));                % Probability of detection at ST due to inaccurate estimation of threshold
    
    
    C_0_th = 0;                                       % Capacity under scenario 0 at SR
    C_1_th = 0;                                       % Capacity under scenario 1 at SR
    R_id_th = zeros(1,length(N_th));                  % Throughput at SR ideal case
    R_ac_th = zeros(1,length(N_th));                  % Throughput at SR with accurate estimation
    R_oc_th = zeros(1,length(N_th));                  % Throughput at SR with upper limit of the energy wall
    
    
    % Accurarte energy  
    Acc_energy = (noise_power + alpha_p_1 * P_p);      
 
    for i=1:length(N_th)
        disp(strcat('N = ',num2str(N_th(i))));     
        
       %% Determining the performance meteric R for all the different cases
       %% Ideal         
       %% Threshold  
        epsilon_id_th(i) = qfuncinv(P_d_d) * sqrt(2/N_th(i)) * Acc_energy +...
            Acc_energy;                  
       
       %% Probability of false alarm and probability of dectection
        P_f_id_th(i) = qfunc((epsilon_id_th(i) - noise_power)/...
            (sqrt(2/N_th(i)) * noise_power));        
       
       %% Expected Rate
        C_0_th = log2(1 + alpha_s * P_s / noise_power); 
        C_1_th = log2(1 + alpha_s * P_s / (P_p * alpha_p_2 +...
            noise_power) );
        P_d_id_th(i) =  P_d_d;
        
        R_id_th(i) = (K - N_th(i))/K * (P_H0 * (1 -  P_f_id_th(i)) * C_0_th +...
            (1 - P_H0) * (1 - P_d_id_th(i)) * C_1_th);  
        

        if N_th(i) > N_est_th
           N_sen_th = N_th(i) - N_est_th;
           %% Adjacent case
           if 1
               %% Threshold  
               num_test_points = 1000;
               test_points = linspace(Acc_energy* (1 + sqrt(2/N_sen_th)* qfuncinv(0.85)), ...
                   Acc_energy* (1 + sqrt(2/N_sen_th)* sqrt(2) * (-erfcinv(10^-320))), num_test_points);
               Exp_P_d_ac = zeros(1,num_test_points);
               for j=1:num_test_points
                    ts = sqrt(2/N_sen_th);
                    te = sqrt(2/N_est_th);
                    func_exp_P_d = @(t)  t .*  test_points(j) .*...
                        exp( (qfuncinv(t)/sqrt(2)).^2 - (test_points(j)./(1 + ts * qfuncinv(t)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
                        * ts ./ (Acc_energy * te * ( 1 + ts * qfuncinv(t)).^2);                    
                    Exp_P_d_ac(j) =  integral(func_exp_P_d, 0, 1);
                    if Exp_P_d_ac(j) >= P_d_d
                       epsilon_ac_th(i) = test_points(j); 
                       break;
                    else
                       epsilon_ac_th(i) = test_points(j); 
                    end
               end


               %% Probability of false alarm and probability of dectection
               P_f_ac_th(i) = qfunc((epsilon_ac_th(i) - noise_power)/...
                   (sqrt(2/N_sen_th) * noise_power));  

               func_exp_P_d = @(t)  t .*  epsilon_ac_th(i) .*...
                    exp( (erfcinv(2 * t)).^2 - (epsilon_ac_th(i)./(1 + sqrt(2) * ts * erfcinv(2 * t)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
                    * ts ./ (Acc_energy * te * ( 1 + sqrt(2) * ts * erfcinv(2 * t)).^2);
               P_d_ac_th(i) = integral(func_exp_P_d, 0, 1);
               if  P_d_ac_th(i) < 0.88
                   ts = sqrt(2/N_sen_th);
                   te = sqrt(2/N_est_th);
                   my_int =  my_integral_pd(Acc_energy, epsilon_ac_th(i), ts, te, eps*10^-8, eps,10^7);
                   P_d_ac_th(i) = integral(func_exp_P_d, 0, 1 -eps) + my_int;
               end


               %% Expected Rate
               C_0_th = log2(1 + alpha_s * P_s / noise_power); 
               C_1_th = log2(1 + alpha_s * P_s / (P_p * alpha_p_2 +...
                   noise_power) );

               R_ac_th(i) = (K - N_th(i))/K * (P_H0 * (1 -  P_f_ac_th(i)) * C_0_th +...
                   (1 - P_H0) * (1 - P_d_ac_th(i)) * C_1_th); 
           end

 
            %% Outage case
            if 1
                %% Threshold  
                epsilon_oc_th(i) = Acc_energy * (1 + qfuncinv(1 - mu) * (sqrt(2/N_est_th)))...
                     * (1 + qfuncinv(P_d_d) * (sqrt(2/N_sen_th)));            

                %% Probability of false alarm and  expected probability of dectection
                P_f_oc_th(i) = qfunc((epsilon_oc_th(i) - noise_power)/...
                    (sqrt(2/N_sen_th) * noise_power));  
                
                %% Test
                if 1 %> 2500
                    if 1
                        ts = sqrt(2/N_sen_th);
                        te = sqrt(2/N_est_th);
                        Acc_energy = (noise_power + alpha_p_1 * P_p);                       
                        
%                         scale = 1;
%                         func_pd = @(t) scale * epsilon_oc_th(i) .* exp( (erfinv(1 - 2 * scale * t)).^2 - (epsilon_oc_th(i)./(1 + ts * sqrt(2) * erfcinv(1 - 2 * scale * t)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
%                         .* ts ./ (Acc_energy * te * ( 1 + sqrt(2) * ts * erfinv(1 - 2 * scale * t)).^2); 
%                         integral(func_pd, 0, 1/scale -eps)  
%                         integral(func_pd, 0, 1/scale, 'Waypoints',[0 1/scale - 10^-6 1/scale])  

                        %% Test Maclurian series approximation
                        if 0
                           x = double(0.7);
                           t = 0.5 * sqrt(pi) * x;
                           erfinv(x)
                           series = double((t + 1/3 * t^3 + 7/30 * t^5 + 127/630 * t^7 + 4369/22680 * t^9))
                            
                        end
                        
                        scale = 0.1;
                        func_pd = @(t) scale * epsilon_oc_th(i) .* exp( (erfcinv(2 * scale * t)).^2 - (epsilon_oc_th(i)./(1 + ts * sqrt(2) * erfcinv(2 * scale * t)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
                        .* ts ./ (Acc_energy * te * ( 1 + sqrt(2) * ts * erfcinv(2 * scale * t)).^2); 
                        integral(func_pd, 0, 1/scale);  
                        integral(func_pd, 0, 1/scale, 'Waypoints',[0 1/scale - 10^-6 1/scale]); 
                        %integral(func_pd, 0, 1/scale - 10^-6) + integral(func_pd, 1/scale - 10^-6, 1/scale) 

                        
                        %scale = 1;                        
                        %func_pd = @(t) scale * epsilon_oc_th(i) .* exp( (qfuncinv(scale * t)/sqrt(2)).^2 - (epsilon_oc_th(i)./(1 + ts * qfuncinv(scale * t)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
                        %.* ts ./ (Acc_energy * te * ( 1 + ts * qfuncinv(scale * t)).^2); 
                        %integral(func_pd, 0, 1/scale - eps)  
                        %integral(func_pd, 0, 1/scale, 'Waypoints',[0 1/scale - 10^-3 1/scale])  

                        %integral(func_pd, 0, 1, 'AbsTol',1e-18)
                        %integral(func_pd, 0, 1/scale - 1) + integral(func_pd, 1/scale - 1, 1/scale -eps)
                        %integral(func_pd, 0, 1-eps) + integral(func_pd, 1-eps, 1)
                        if 0
                            close all;
                            sample_pts = 1000; 
                            x = 1/sample_pts:1/sample_pts:1-eps;
                            test_value = linspace(1- 10^-6,1, 1000);
                            func_value = epsilon_oc_th(i) .* exp( (erfcinv(2 * test_value)).^2 - (epsilon_oc_th(i)./(1 + ts * erfcinv(2 * test_value)) - Acc_energy).^2 /(Acc_energy * te)^2) * ts...
                                ./ (Acc_energy * te * ( 1 + ts * erfcinv(2 * test_value)).^2)
                            figure(1);
                            plot(test_value, ts * erfcinv(2 * test_value));
                            figure(2);
                            plot(test_value,( 1 + sqrt(2) * ts * erfcinv(2 * test_value)).^2,'r');
                            hold on,
                            plot(test_value,(sqrt(2) * ts * erfcinv(2 * test_value)).^2,'b');
                            figure(3);
                            plot(test_value,(erfcinv(2 * test_value)).^2 - (epsilon_oc_th(i)./(1 + ts * erfcinv(2 * test_value)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2,'r');
                            hold on,
                            plot(test_value,(erfcinv(2 * test_value)).^2,'b');
                            desi_int = sum(1/sample_pts * epsilon_oc_th(i) .* exp( (qfuncinv(x)/sqrt(2)).^2 - (epsilon_oc_th(i)./(1 + ts * qfuncinv(x)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
                                .* ts ./ (Acc_energy * te * ( 1 + ts * qfuncinv(x).^2)));
                            plot(epsilon_oc_th(i) .* exp( (qfuncinv(x)/sqrt(2)).^2 - (epsilon_oc_th(i)./(1 + ts * qfuncinv(x)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
                                .* ts ./ (Acc_energy * te * ( 1 + ts * qfuncinv(x).^2)));
                        end
                        
                    end
                end

                ts = sqrt(2/N_sen_th);
                te = sqrt(2/N_est_th);
                func_exp_P_d = @(t) epsilon_oc_th(i) .*...
                    exp( (qfuncinv(t)/sqrt(2)).^2 - (epsilon_oc_th(i)./(1 + ts * qfuncinv(t)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
                    * ts ./ (Acc_energy * te * ( 1 + ts * qfuncinv(t)).^2);
                int_11 = integral(func_exp_P_d, 0, 1);  
                    
                %% Integration to test the covergenece of the probability density  
                func_P_d_1 = @(t)  epsilon_oc_th(i) .*...
                    exp( (erfcinv(2 * t)).^2 - (epsilon_oc_th(i)./(1 + sqrt(2) * ts * erfcinv(2 * t)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
                    * ts ./ (Acc_energy * te * ( 1 + sqrt(2) * ts * erfcinv(2 * t)).^2);
                int_12 = integral(func_P_d_1, 0, 1)   
                
                func_P_d_2 = @(t)  epsilon_oc_th(i) .*...
                    exp( ((-erfcinv(2 - 2 * t))).^2 - (epsilon_oc_th(i)./(1 + sqrt(2) * ts * (-erfcinv(2 - 2 * t))) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
                    * ts ./ (Acc_energy * te * ( 1 + sqrt(2) * ts * (-erfcinv(2 - 2 * t))).^2);
                my_int = ...%my_integral(Acc_energy, epsilon_oc_th(i), ts, te, eps^6, eps^5,10^6) +...
                    ...%my_integral(Acc_energy, epsilon_oc_th(i), ts, te, eps^5, eps^4,10^6) +...
                    ...%my_integral(Acc_energy, epsilon_oc_th(i), ts, te, eps^4, eps^3,10^6) +...
                    ...%my_integral(Acc_energy, epsilon_oc_th(i), ts, te, eps^3, eps^2,10^6 +...
                    my_integral(Acc_energy, epsilon_oc_th(i), ts, te, eps*10^-8, eps,10^7)
                    
                    
                if my_int < 1
                    int_13 = integral(func_P_d_1, 0, 1 -eps) +  my_int
                else
                    int_13 = integral(func_P_d_1, 0, 1 -eps)
                end
                %my_integral(Acc_energy, epsilon_oc_th(i), ts, te, eps^5, eps^2,10^6)  
                %% Integration to test the covergenece of the expected value  
                func_exp_P_d = @(t)  t .*  epsilon_oc_th(i) .*...
                    exp( (erfcinv(2 * t)).^2 - (epsilon_oc_th(i)./(1 + sqrt(2) * ts * erfcinv(2 * t)) - Acc_energy).^2 /(sqrt(2) * Acc_energy * te)^2)...
                    * ts ./ (Acc_energy * te * ( 1 + sqrt(2) * ts * erfcinv(2 * t)).^2);
                %my_integral_pd(Acc_energy, epsilon_oc_th(i), ts, te, eps*10^-8, eps,10^7)
                my_int =  my_integral_pd(Acc_energy, epsilon_oc_th(i), ts, te, eps*10^-8, eps,10^7)
                P_d_oc_th(i) =  integral(func_exp_P_d, 0, 1-eps) + my_int;
                integral(func_exp_P_d, 0, 1)
                if 0
                    a = 1 - 10^-12;
                    int1 = integral(func_P_d, 0, 0.99);
                    int2 = 1 - int1;
                    integral(func_exp_P_d, 0, 0.99);
                    P_d_oc_th(i) = integral(func_exp_P_d, 0, a) + 0.113425 *0.5 * epsilon_oc_th(i)./(Acc_energy *  ts * te); 
                end
               
                
                %integral(func_exp_P_d, 0, 1)
                %integral(func_exp_P_d, 0, 1, 'Waypoints',[0 1/scale - 10^-6 1/scale])

                %% Expected Rate
                C_0_th = log2(1 + alpha_s * P_s / noise_power); 
                C_1_th = log2(1 + alpha_s * P_s / (P_p * alpha_p_2 +...
                    noise_power) );

                R_oc_th(i) = (K - N_th(i))/K * (P_H0 * (1 -  P_f_oc_th(i)) * C_0_th +...
                (1 - P_H0) * (1 - P_d_oc_th(i)) * C_1_th);
            end
        end        
        
        disp(strcat('epsilon_id_th(i) = ',num2str(epsilon_id_th(i))));
        disp(strcat('epsilon_ac_th(i) = ',num2str(epsilon_ac_th(i))));
        disp(strcat('epsilon_oc_th(i) = ',num2str(epsilon_oc_th(i))));
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
    N_id_opt = N_th(index);
    
    % Optimum throughput for the estimation model -- average constraint case
    [R_ac_opt index] = max(R_ac_th);
    N_ac_opt = N_th(index);
    
    % Optimum throughput for the estimation model -- outage constraint case
    [R_oc_opt index] = max(R_oc_th);
    N_oc_opt = N_th(index);
    
    %save('results_thr_sen_time_tradeoff_AWGN_wo_nu_snr_m10_th.mat');
    %quit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Sensing throughput tradeoff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load('results_thr_sen_time_tradeoff_AWGN_wo_nu_snr_m10_th.mat');
Fontsize = 9;
plot(N_th * 1e-3, R_id_th, 'b', 'LineWidth',1.5);
hold on,
plot(N_th * 1e-3, R_ac_th,'LineStyle', '--', 'LineWidth', 1.5);
hold on,
plot(N_th * 1e-3, R_oc_th, 'LineStyle', '-.', 'LineWidth', 1.5);
hold on,
plot(N_id_opt * 1e-3, R_id_opt, 'ks', 'LineWidth',1);
hold on,
plot(N_ac_opt * 1e-3, R_ac_opt, 'ks', 'LineWidth',1, 'HandleVisibility','off');
hold on,
plot(N_oc_opt * 1e-3, R_oc_opt, 'ks', 'LineWidth',1, 'HandleVisibility','off');
if 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves simulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('results_thr_sen_time_tradeoff_AWGN_wo_nu_snr_m10_sim.mat');
    hold on,
    plot(N_sim * 1e-3, R_id_sim, 'ro', 'LineWidth',1);
    hold on,
    plot(N_sim * 1e-3, R_ac_sim, 'ro', 'LineWidth',1, 'HandleVisibility','off');
    hold on,
    plot(N_sim * 1e-3, R_oc_sim, 'ro', 'LineWidth',1, 'HandleVisibility','off');
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves cosmetic makeover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid on;
axis([0.1 max(N_th) * 1e-3 0  max(R_id_opt) * 1.02]);
ylabel('$\rs =$ [bits/sec/Hz]','FontSize',Fontsize);
xlabel('$\test + \tsen$ = [ms]','FontSize',Fontsize);
hl = legend('$\rs$', '$\rsac$', '$\rsoc$', 'opt', 'sim');
set(hl, 'Location', 'SouthEast', 'FontSize', Fontsize);
set(gca,'FontSize',Fontsize);
%laprint(1, '../figures/fig_thr_sen_time_tradeoff_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
%    'on', 'factor',0.5, 'keepfontprops', 'on');