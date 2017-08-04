%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: This m-file considers the Throughput-Sesning time 
%        tradeoff with SNR estimation for the fading channel. In this m-script  
%        the long term analysis is investigted that is, considering invidual 
%        frame. Hence, the SNR is considered constant for the time interval 
%        (tau_est + tau_sen). 
%
%        Clearly the recieved power is distributed according to non central
%        chi-2 distribution, which can be approximated as Gaussian 
%        distribution for large sample count (This again needs an analysis, 
%        on how good is the approximation is?)
%        
%        The simulation is performed to show the following:
%        1) Analyse the throughput vs sensing for different cases 
%           - Ideal case, when SNR is known at the ST (i)
%           - Outage Constraint   
%        2) Confirm the theoretical analysis vs simulation results
%
%        Created on: 17.01.16
%        Last modified: 17.01.16
%        Revision History: 17.01.16 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
%warning off;


sim = 00;                                                                       % Enable( = 1) to perform simulation, theoretical analysis
th = 00;                                                                        % disable ( = 0) to plot curves, data is read from a file
                                                                                % Theoretical anaylsis is also included of simulation as
                                                                                % numerical integration is involved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_s = 10^(-10/10);                                                              % Power transmitted by ST, the SNR received at SR 
P_p = 10^(-00/10);                                                              % Power transmitted by PT, the SNR received at SR 
                                                                                % varried using this parameter
noise_power = 10^(-100/10);                                                     % noise power -100 dBm
f_s = 1e6;                                                                      % 1 MHz one band
K = 0.1 * f_s;                                                                  % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
g_p1_true = 10^(-100/10);                                                       % True Path loss between ST and PR   
g_p2_true = 10^(-100/10);                                                       % True Path loss between PR and SR   
g_s_true = 10^(-080/10);                                                        % True Path loss between ST and SR 
P_H0 = 0.8;                                                                     % Probability of Hypothesis 0
P_d_d = 0.90;                                                                   % Constraint on Probability of detection P_d
mu = 0.05;                                                                       % Outage Probability on probability of detection
m_s = 1;%[0.5 1 2 5 100];                                                          % Nakagam-m parameter
m_p2 = 1;%[0.5 1 2 5 100];                                                         % Nakagam-m parameter
N_s = 1;                                                                       % Number of Pilot Symobols 
est_er = noise_power/N_s;                                                       % Variance of the Estimation error for h_s 

num_test_points = 10000;                    % Test points for evaluating threshold 

%[C_0, C_1] = calc_capacities(alpha_p_2 * P_p / noise_power,... 
%    alpha_s * P_s / noise_power,....
%    noise_power);                           % Capacities at SR, with and without intereference from PT  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sim
     M = 1e5;                                                                   % Number of realizations 
     %% Simulations parameters
     N_est_sim = 1000;                                                          % N_est_th = number of samples used for estimation = tau * f_s 
     N_sen_sim = [50:500:15000];                                              % N = Total of samples used for sensing = tau * f_s 

     epsilon_id_sim = zeros(1,length(N_sen_sim));                               % Ideal threshold 
     epsilon_oc_sim = zeros(1,length(N_sen_sim));                               % Inaccurate threshold due to incorrect estimation of received power

     P_f_id_sim = zeros(1,length(N_sen_sim));                                   % Probability of false alarm at ST for the ideal case
     P_f_oc_sim = zeros(1,length(N_sen_sim));                                   % Probability of false alarm at ST due to inaccurate estimation of threshold

     P_d_id_sim = zeros(1,length(N_sen_sim));                                   % Probability of detection at ST for the ideal case
     P_d_oc_sim = zeros(1,length(N_sen_sim));                                   % Probability of detection at ST due to inaccurate estimation of threshold
    
     R_id_sim = zeros(1,length(N_sen_sim));                                     % Throughput at SR ideal case
     R_oc_sim = zeros(1,length(N_sen_sim));                                     % Throughput at SR with upper limit of the energy wall

 
          
     
    %% Computation of the estimation 
     % Accurarte energy  
     %Acc_energy = (noise_power + alpha_p_1 * P_p);
     
     for i=1:length(N_sen_sim)
        disp(strcat('N = ',num2str(N_sen_sim(i)))); 
        P_rcvd_ST_H1_sen_sim = zeros(1,M);                                      % Rcvd power samples under Hypothesis 1 (sensing) at ST 
        P_rcvd_ST_H0_sen_sim = zeros(1,M);                                      % Rcvd power samples under Hypothesis 0 (sensing) at ST
        P_rcvd_ST_H1_est_sim = zeros(1,M);                                      % Rcvd power samples under Hypothesis 1 (estimation) at ST
        %P_rcvd_H0_est_sim = zeros(1,M);                                        % Rcvd power samples under Hypothesis 0 (estimation)
        P_rcvd_SR_H1_est_sim = zeros(1,M);                                      % Rcvd power samples under Hypothesis 1 (estimation) at SR 
        g_s_est = zeros(1,M);                                                   % Pilot based estimation at SR 

        
        g_p1_sim = random('gam',m_p2, g_p1_true/m_p2, 1, M);                    % Power gain for the channel g_p2 (including path loss)
        g_p2_sim = random('gam',m_p2, g_p2_true/m_p2, 1, M);                    % Power gain for the channel g_p2 (including path loss)
        g_s_sim = random('gam', m_s, g_s_true/m_s, 1, M);                       % Power gain for the channel g_s (including path loss)   

       %% Ideal Case -- Remember (no estimation)    
        parfor j=1:M
            P_rcvd_ST_H1_sen_sim(j) = mean(random('norm',...
                random('norm',0, sqrt(P_p * g_p1_sim(j)), 1, N_sen_sim(i)),...
                sqrt(noise_power), 1, N_sen_sim(i)).^2);
        end 
        parfor j=1:M  
            P_rcvd_ST_H0_sen_sim(j) = mean(random('norm', 0, sqrt(noise_power), 1, N_sen_sim(i)).^2);
        end 


       %% Determining the performance meteric R for all the different cases
       %% Threshold     
        test_points = linspace(max(P_rcvd_ST_H1_sen_sim), min(P_rcvd_ST_H1_sen_sim),...
            num_test_points);
        
        if 1   % 1) Determining threshold by considering the outage constraint on the P_d
            
            for k= 1:num_test_points
                P_d_test_points = gammainc(N_sen_sim(i) * test_points(k)./...
                 (2 * (P_p * g_p1_sim + noise_power)), N_sen_sim(i)/2, 'upper');
                if (length(find(P_d_test_points <= P_d_d))/M <= mu)    
                    epsilon_id_sim(i) = test_points(k);
                    P_d_id_sim(i) = mean(P_d_test_points);
                    break;
                end
            end
            
        else   % 2) Conventional approach: Determining threshold directly from (P_rcvd_H1 > mu) = P_d_d      
            for k = 1:num_test_points
                if (length(find(P_rcvd_ST_H1_sen_sim > test_points(k)))/M) <= P_d_d   
                    epsilon_id_sim(i) = test_points(k);
                    P_d_id_sim(i) = length(find(P_rcvd_ST_H1_sen_sim > test_points(k)))/M;
                    break;
                end
            end   
        end

       %% Probability of false alarm and probability of detection
        P_f_id_sim(i) = length(find(P_rcvd_ST_H0_sen_sim > epsilon_id_sim(i)))/M;

       %% Expected Rate
        C_0_id_sim = mean(log2(1 + g_s_sim * P_s / noise_power)); 
        C_1_id_sim = mean(log2(1 + g_s_sim * P_s ./ (P_p * g_p2_sim +...
            noise_power)));
        
        R_id_sim(i) = (K - N_sen_sim(i))/K * (P_H0 * (1 -  P_f_id_sim(i)) *...
            C_0_id_sim +  (1 - P_H0) * (1 - P_d_id_sim(i)) * C_1_id_sim);
        
       %% non ideal case -- includes estimation followed by sesning
        if N_sen_sim(i) > N_est_sim
            N_sen_ni_sim = N_sen_sim(i);% - N_est_sim;        % Sensing time for non ideal case

           %% Estimation Phase
            parfor j=1:M            
                P_rcvd_ST_H1_est_sim(j) = mean(random('norm',...
                    random('norm',0, sqrt(P_p * g_p1_sim(j)), 1, N_est_sim),...
                    sqrt(noise_power), 1, N_est_sim).^2);   
                P_rcvd_SR_H1_est_sim(j) = mean(random('norm',...
                    random('norm',0, sqrt(P_p * g_p2_sim(j)), 1, N_est_sim),...
                    sqrt(noise_power), 1, N_est_sim).^2);
                g_s_est(j) = mean(normrnd(sqrt(P_s * g_s_sim(j)), sqrt(noise_power), 1, N_s).^2);

            end 
            
           %% Generate test points to determine the threshold
            test_points = linspace(epsilon_id_sim(i), 0.1 * epsilon_id_sim(i), num_test_points);

          
           %% Outage constraint            
            if 1
              %% Threshold                  
                for k = 1:length(test_points)
                    P_d_test_points = gammainc(N_sen_ni_sim * test_points(k)./...
                        (2 * P_rcvd_ST_H1_est_sim), N_sen_ni_sim/2, 'upper');
                    if (length(find(P_d_test_points <= P_d_d))/M <= mu)    
                        epsilon_oc_sim(i) = test_points(k);
                      %% Expected probability of detection is calculated
                      %% as P_d is random variable;
                        P_d_oc_sim(i) = mean(P_d_test_points);
                        break;
                    end
                end
                   
                C_0_oc_sim = mean(log2(1 + g_s_est/noise_power));
                C_1_oc_sim = mean(log2(1 + g_s_est./P_rcvd_SR_H1_est_sim));
                
                
              %% Probability of false alarm and probability of detection
                P_f_oc_sim(i) = length(find(P_rcvd_ST_H0_sen_sim > epsilon_oc_sim(i)))/M;

                R_oc_sim(i) = (K - N_sen_sim(i))/K * (P_H0 * (1 -  P_f_oc_sim(i)) *...
                    C_0_oc_sim +  (1 - P_H0) * (1 - P_d_oc_sim(i)) * C_1_oc_sim);   
            end
        end
            
        disp(strcat('epsilon_id_sim(i) = ',num2str(epsilon_id_sim(i))));   
        disp(strcat('epsilon_oc_sim(i) = ',num2str(epsilon_oc_sim(i))));   
        disp(strcat('P_f_id_sim(i) = ',num2str(P_f_id_sim(i)))); 
        disp(strcat('P_f_oc_sim(i) = ',num2str(P_f_oc_sim(i))));
        disp(strcat('P_d_id_sim(i) = ',num2str(P_d_id_sim(i)))); 
        disp(strcat('P_d_oc_sim(i) = ',num2str(P_d_oc_sim(i)))); 
        disp(strcat('R_id_sim(i) = ',num2str(R_id_sim(i)))); 
        disp(strcat('R_oc_sim(i) = ',num2str(R_oc_sim(i)))); 
    	save('results_thr_sen_time_tradeoff_fading_m1_sim.mat');
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
    N_est_th = 1000;                                  % N_est_th = number of samples used for estimation = tau * f_s 
    N_sen_th = [10:15:1000,1001,1100:50:15000];    % N = includes the time in samples, for the ideal case their is no
                                                      % estimation time hence sensing time = N_th and the throughput will 
                                                      % attain a nonzero value. For the accurate, the energy walls sensing 
                                                      % time begin after the estimation of the rcvd energy is performed.
                                                      % N_sen = N_th - N_est
%    epsilon_id_th = zeros(1,length(N_sen_th));        % Ideal threshold 
%    epsilon_ac_th = zeros(1,length(N_sen_th));        % AC threshold 
%    epsilon_oc_th = zeros(1,length(N_sen_th));        % OC threshold

    P_f_id_th = zeros(1,length(N_sen_th));            % Probability of false alarmC_0 at ST for the ideal case
    P_f_oc_th = zeros(1,length(N_sen_th));            % Probability of false alarm at ST due to inaccurate estimation of threshold

    P_d_id_th = zeros(1,length(N_sen_th));            % Probability of detection at ST for the ideal case
    P_d_oc_th = zeros(1,length(N_sen_th));            % Probability of detection at ST due to inaccurate estimation of threshold
    
    
    R_id_th = zeros(1,length(N_sen_th));              % Throughput at SR ideal case
    R_ac_th = zeros(1,length(N_sen_th));              % Throughput at SR with accurate estimation
    R_oc_th = zeros(1,length(N_sen_th));              % Throughput at SR with upper limit of the energy wall
    
   %% Expected Rate
    func_C0 = @(t) log2(1 + t * g_s_true * P_s / noise_power) .*...
        1/gamma(m_s) * m_s^(m_s) .* t.^(m_s - 1) .* exp(-m_s * t);
    func_C1 = @(t, u) log2(1 + t * g_s_true * P_s ./ (P_p * u * g_p2_true + noise_power)) .*...
        1/gamma(m_s) * m_s^(m_s) .* t.^(m_s - 1) .* exp(-m_s * t) .*...
        1/gamma(m_p2) * m_p2^(m_p2) .* u.^(m_p2 - 1) .* exp(-m_p2 * u);
    C_0_id_th = integral(func_C0, 0, 100);
    C_1_id_th = integral2(func_C1, 0, 100, 0, 100);
                    
    [C_0_oc_th, C_1_oc_th] = calc_capacities_fad(P_p, P_s, g_p2_true, g_s_true,...
        noise_power, noise_power/N_s, N_est_th, N_s, m_p2, m_s);
 
    parfor i=1:length(N_sen_th)
        disp(strcat('N = ',num2str(N_sen_th(i))));     
        M = 1e2;                                                                   % Number of realizations 
        epsilon_id_th = 0;
        epsilon_oc_th = 0;
        P_rcvd_sim = zeros(1, M);

        %% Determining the performance meteric R for 0all the different cases
        %% Ideal         
        %% Threshold  
         if 1   % 1) Determining threshold by considering the outage constraint on the P_d
            epsilon_id_th = (2/N_sen_th(i)*(P_p * g_p1_true * gammaincinv(1 - mu, m_p2, 'upper')/m_p2 +...
            noise_power) * gammaincinv(P_d_d, N_sen_th(i)/2, 'upper')); 
            
            %% Probability of detection
            func_exp_P_d = @(t) gammainc(N_sen_th(i) * epsilon_id_th./...
                (2 *( t * P_p * g_p1_true + noise_power)), N_sen_th(i)/2, 'upper') .*...
                1/gamma(m_p2) * m_p2^(m_p2) .* t.^(m_p2 - 1) .* exp(-m_p2 * t);
            P_d_id_th(i) = integral(func_exp_P_d, 0, 100);
         
         else  % 2) Conventional approach: Determining threshold directly from (P_rcvd_H1 > mu) = P_d_d     
            g_p1_sim = random('gam',m_p2, g_p1_true/m_p2, 1, M);                       % Power gain for the channel g_p2 (including path loss)
             for k=1:M
                P_rcvd_sim(k) = mean(random('norm',...
                    random('norm',0, sqrt(P_p * g_p1_sim(k)), 1, N_sen_th(i)),...
                    sqrt(noise_power), 1, N_sen_th(i)).^2);
             end
             test_points = linspace(min(P_rcvd_sim), max(P_rcvd_sim), num_test_points);
             P_d_id = zeros(1,num_test_points);
             for k=1:num_test_points
                func_P_rcvd = @(t) gammainc(N_sen_th(i) * test_points(k)./...
                    (2 * (t * P_p * g_p2_true +  noise_power)), N_sen_th(i)/2, 'upper') .*...
                    1/gamma(m_p2) * m_p2^(m_p2) .* t.^(m_p2 - 1) .* exp(- m_p2 * t);
                P_d_id(k) = integral(func_P_rcvd, 0, 100);
                if P_d_id(k) < P_d_d
                    epsilon_id_th = test_points(k);
                    P_d_id_th(i) = P_d_id(k);
                    break;
                end
             end
        end 
       
        %% Probability of false alarm 
        P_f_id_th(i) = gammainc(N_sen_th(i) * epsilon_id_th/(2 * noise_power), N_sen_th(i)/2, 'upper');
      

        
        R_id_th(i) = (K - N_sen_th(i))/K * (P_H0 * (1 -  P_f_id_th(i)) * C_0_id_th +...
            (1 - P_H0) * (1 - P_d_id_th(i)) * C_1_id_th);  
        

        if N_sen_th(i) >= N_est_th
           N_sen_ni_th = N_sen_th(i);% - N_est_th;

           %% Outage case
            if 1
              %% Threshold  
                test_points = linspace(0.1 * epsilon_id_th, 2 * epsilon_id_th, num_test_points);
                P_d_oc = zeros(1,num_test_points);
                for j=1:num_test_points
                   func_P_d = @(t) gammainc(N_est_th * N_sen_ni_th * test_points(j)./...
                       (4 * (P_p * g_p1_true * t + noise_power) * gammaincinv(P_d_d, N_sen_ni_th/2, 'upper')),  N_est_th/2, 'upper') .*...
                       1/gamma(m_p2) * m_p2^(m_p2) .* t.^(m_p2 - 1) .* exp(- m_p2 * t); 
                   P_d_oc(j) = integral(func_P_d, 0, 100);
                   %P_d_oc(j) = gammainc(N_est_th * N_sen_ni_th * test_points(j)./...
                   %    (4 * (P_p * g_p1_true + noise_power) * gammaincinv(P_d_d, N_sen_ni_th/2, 'upper')),  N_est_th/2, 'upper') .*...
                   %    1; 
                   
                    if P_d_oc(j) <= 1 - mu 
                        epsilon_oc_th = test_points(j);
                        %epsilon_oc_or = 4 * (P_p * g_p1_true + noise_power) * gammaincinv(P_d_d, N_sen_ni_th/2, 'upper')...
                        %   * gammaincinv(1 - mu, N_est_th/2, 'upper') / (N_sen_ni_th * N_est_th);
                        if 0
                            x_pts = linspace(1e-6,1-1e-6,10000);
                            PDF_pd_ana = zeros(1, length(x_pts));
                            for k=1:length(x_pts)
                                func_pd = @(t) exp(gammaln(N_sen_ni_th/2) - gammaln(N_est_th/2)) *...  
                                    exp(-epsilon_oc_th * N_est_th * N_sen_ni_th./ (4 * (t * P_p * g_p1_true + noise_power) *...
                                   gammaincinv(x_pts(k),N_sen_ni_th/2,'upper')) + gammaincinv(x_pts(k),N_sen_ni_th/2,'upper')) .*...
                                   (epsilon_oc_th * N_sen_ni_th * N_est_th./ (4 * (t * P_p * g_p1_true + noise_power) *...
                                   gammaincinv(x_pts(k),N_sen_ni_th/2,'upper'))).^(N_est_th/2) .*...
                                    (1/gammaincinv(x_pts(k),N_sen_ni_th/2,'upper'))^(N_sen_ni_th/2) .* ... 
                                   1/gamma(m_p2) * m_p2^(m_p2) .* t.^(m_p2 - 1) .* exp(-m_p2 * t);
                                PDF_pd_ana(k) = integral(func_pd, 0 , 100);    
                                
                            end
                            plot(x_pts, PDF_pd_ana * ((x_pts(3) - x_pts(2))), 'r', 'Linewidth',2.5);
                        end
                        
                      %% Evaluate the intergral using the distribution function
                        if 1
                            x_pts = linspace(0,1,10000);
                            CDF_pts = zeros(1,length(x_pts));
                            for k=1:length(x_pts)
                                func_P_d = @(t) 1 - gammainc(N_est_th * N_sen_ni_th * epsilon_oc_th./...
                                    (4 * (P_p * g_p1_true * t + noise_power) * gammaincinv(x_pts(k), N_sen_ni_th/2, 'upper')),  N_est_th/2, 'upper') .*...
                                    1/gamma(m_p2) * m_p2^(m_p2) .* t.^(m_p2 - 1) .* exp(- m_p2 * t); 
                                CDF_pts(k) = integral(func_P_d, 0, 100);
                            end
                            PDF_pts = diff(CDF_pts);
                            P_d_oc_th(i) = sum(PDF_pts .* x_pts(2:end));
                        else
                        
                            func_P_d = @(u, t) u .* exp(gammaln(N_sen_ni_th/2) - gammaln(N_est_th/2) -... 
                                epsilon_oc_th * N_est_th * N_sen_ni_th./ (4 * (t * P_p * g_p1_true + noise_power) .*...
                                gammaincinv(u ,N_sen_ni_th/2,'upper')) + gammaincinv(u, N_sen_ni_th/2,'upper') +...
                                N_est_th/2 * log(epsilon_oc_th * N_sen_ni_th * N_est_th./ (4 * (t * P_p * g_p1_true + noise_power) .*...
                                gammaincinv(u, N_sen_ni_th/2,'upper'))) +...
                                (N_sen_ni_th/2) * log(1./gammaincinv(u, N_sen_ni_th/2,'upper'))) .*...
                                1/gamma(m_p2) * m_p2^(m_p2) .* t.^(m_p2 - 1) .* exp(-m_p2 * t);
                             P_d_oc_th(i) = integral2(func_P_d, 0, 1, 0, 100);
                        end
                        break;
                    end
                end

              %% Probability of false alarm and expected probability of dectection
                P_f_oc_th(i) = gammainc(N_sen_ni_th * epsilon_oc_th/(2 * noise_power), N_sen_ni_th/2, 'upper');
                
                R_oc_th(i) = (K - N_sen_th(i))/K * (P_H0 * (1 -  P_f_oc_th(i)) * C_0_oc_th +...
                (1 - P_H0) * (1 - P_d_oc_th(i)) * C_1_oc_th);
            end
        end        
        
        disp(strcat('epsilon_id_th = ',num2str(epsilon_id_th)));
        disp(strcat('epsilon_oc_th = ',num2str(epsilon_oc_th)));
        disp(strcat('P_f_id_th(i) = ',num2str(P_f_id_th(i)))); 
        disp(strcat('P_f_oc_th(i) = ',num2str(P_f_oc_th(i))));
        disp(strcat('P_d_id_th(i) = ',num2str(P_d_id_th(i)))); 
        disp(strcat('P_d_oc_th(i) = ',num2str(P_d_oc_th(i)))); 
        disp(strcat('R_id_th(i) = ',num2str(R_id_th(i)))); 
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
    save('results_thr_sen_time_tradeoff_fading_m1_5_th.mat');
    quit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Sensing throughput tradeoff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_thr_sen_time_tradeoff_fading_m1_th.mat');
Fontsize = 8;
[values index] = find(N_sen_th <= 20000);
h1 = plot(N_sen_th(index) * 1e-3, R_id_th(index), 'c', 'LineWidth',2);
hold on,
h2 = plot(N_sen_th(index) * 1e-3, R_oc_th(index), 'LineStyle', '-', 'LineWidth', 1);

if 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves simulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('results_thr_sen_time_tradeoff_fading_m1_sim.mat');
    diff = 3;
    index = 1:diff:length(N_sen_sim);
    hold on,
    h3 = plot(N_sen_sim(index) * 1e-3, R_id_sim(index), 'kx', 'LineWidth',1);
    hold on,
    plot(N_sen_sim(index) * 1e-3, R_oc_sim(index), 'kx', 'LineWidth',1, 'HandleVisibility','off');
end 
hold on,
h4 = plot(N_id_opt * 1e-3, R_id_opt, 'rs', 'LineWidth',1);
hold on,
plot(N_oc_opt * 1e-3, R_oc_opt, 'rs', 'LineWidth',1, 'HandleVisibility','off');

%% m  = 1.5
load('results_thr_sen_time_tradeoff_fading_m1_5_th.mat');
Fontsize = 8;
[values index] = find(N_sen_th <= 20000);
h1 = plot(N_sen_th(index) * 1e-3, R_id_th(index), 'c', 'LineWidth',2);
hold on,
h2 = plot(N_sen_th(index) * 1e-3, R_oc_th(index), 'LineStyle', '-', 'LineWidth', 1);

if 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves simulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('results_thr_sen_time_tradeoff_fading_m1_5_sim.mat');
    diff = 3;
    index = 1:diff:length(N_sen_sim);
    hold on,
    h3 = plot(N_sen_sim(index) * 1e-3, R_id_sim(index), 'kx', 'LineWidth',1);
    hold on,
    plot(N_sen_sim(index) * 1e-3, R_oc_sim(index), 'kx', 'LineWidth',1, 'HandleVisibility','off');
end 
hold on,
h4 = plot(N_id_opt * 1e-3, R_id_opt, 'rs', 'LineWidth',1);
hold on,
plot(N_oc_opt * 1e-3, R_oc_opt, 'rs', 'LineWidth',1, 'HandleVisibility','off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves cosmetic makeover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid on;
axis([0 max(N_sen_th) * 1e-3 0  max(R_id_opt) * 1.03]);
ylabel('$\rs(\test = \SI{1}{ms}, \tsen)$ [bits/sec/Hz]','FontSize',Fontsize);
xlabel('$\tsen$ [ms]','FontSize',Fontsize);
%hl = legend([h1 h2 h3 h5 h4],'$\rs$', '$\rsac$', '$\rsoc$', '$\ttsen$', 'Sim');
hl = legend([h1 h2 h4 h3],'IM, Problem 3', 'EM, Problem 4', '$\trs(\test,\ttsen)$', 'Simulated');
%set(hl, 'position',[0.725 0.12 0.15 0.31]);
set(hl, 'position',[0.58 0.12 0.32 0.28]);
set(gca,'FontSize',Fontsize);

laprint(1, '../figures/fig_thr_sen_time_tradeoff_fading', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');
