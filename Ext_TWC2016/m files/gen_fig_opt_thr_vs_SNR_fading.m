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
%        Created on: 22.01.16
%        Last modified: 22.01.16
%        Revision History: 22.01.16 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
%warning off;


                                                                                % Enable( = 1) to perform simulation, theoretical analysis
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
g_p1_true = 10.^(-[115:-0.5:90]/10);                                                % True Path loss between ST and PR   
g_p2_true = 10^(-100/10);                                                       % True Path loss between PR and SR   
g_s_true = 10^(-080/10);                                                        % True Path loss between ST and SR 
P_H0 = 0.8;                                                                     % Probability of Hypothesis 0
P_d_d = 0.90;                                                                   % Constraint on Probability of detection P_d
mu = 0.05;                                                                      % Outage Probability on probability of detection
m_s = [1 1.5];                                                                  % Nakagam-m parameter
m_p2 = [1 1.5];                                                                 % Nakagam-m parameter
N_s = 1;                                                                        % Number of Pilot Symobols 
est_er = noise_power/N_s;                                                       % Variance of the Estimation error for h_s 

num_test_points = 10000;                                                        % Test points for evaluating threshold 

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th
   %% Theoretical parameters
    % m  = [0.5 0.8]
    N_est_th = [1000];                                                           % N_est_th = number of samples used for estimation = tau * f_s 
    N_sen_th = [100:10:500,550:50:5000,5000:100:10000,10500:500:35000];                                                    % N = includes the time in samples, for the ideal3case their is no
    
    % m  = [1 1.5]
    %N_est_th = [100:50:500,600:200:12000];                                      % N_est_th = number of samples used for estimation = tau * f_s 
    %N_sen_th = [100:100:25000];                                                 % N = includes the time in samples, for the ideal3case their is no
    
    %m = [2 3]
    %N_est_th = [10:10:50,60:30:2500];                                           % N_est_th = number of samples used for estimation = tau * f_s 
    %N_sen_th = [10:10:2600];                                                    % N = includes the time in samples, for the ideal case their is no
    
    %m = [4 5 10 50 100]
    %N_est_th = [10:10:50,60:20:600];                                            % N_est_th = number of samples used for estimation = tau * f_s 
    %N_sen_th = [10:5:600];                                                      % N = includes the time in samples, for the ideal case their is no
    
    % estimation time hence sensing time = N_th and the throughput will 
                                                                                % attain a nonzero value. For the accurate, the energy walls sensing 
                                                                                % time begin after the estimation of the rcvd energy is performed.
                                                                                % N_sen = N_th - N_est
    P_f_id_th = zeros(length(m_p2),length(g_p1_true),length(N_sen_th));         % Probability of false alarm ideal case
    P_f_oc_th = zeros(length(m_p2),length(g_p1_true),length(N_sen_th));         % Probability of false alarm outage case

    P_d_id_th = zeros(length(m_p2),length(g_p1_true),length(N_sen_th));         % Probability of detection at ST ideal case
    P_d_oc_th = zeros(length(m_p2),length(g_p1_true),length(N_sen_th));         % Probability of detection at ST outage case
    
    
    R_id_th = zeros(length(m_p2),length(g_p1_true),length(N_sen_th));           % Throughput at SR ideal case
    R_oc_th = zeros(length(m_p2),length(g_p1_true),length(N_sen_th));           % Throughput at SR outage case

    P_f_id_opt_th = zeros(length(m_p2),length(g_p1_true));                      % Probability of false alarm ideal case @ Opt sen time
    P_f_oc_opt_th = zeros(length(m_p2),length(g_p1_true));                      % Probability of false alarm outage case @ Opt sen time

    P_d_id_opt_th = zeros(length(m_p2),length(g_p1_true));                      % Probability of detection at ST ideal case @ Opt sen time
    P_d_oc_opt_th = zeros(length(m_p2),length(g_p1_true));                      % Probability of detection at ST outage case @ Opt sen time
    
    
    R_id_opt_th = zeros(length(m_p2), length(g_p1_true));                       % Throughput at SR ideal case @ Opt sen time
    R_oc_opt_th = zeros(length(m_p2), length(g_p1_true));                       % Throughput at SR outage case @ Opt sen time

    N_id_opt_th = zeros(length(m_p2), length(g_p1_true));                       % Sensing time at SR ideal case @ Opt sen time
    N_oc_opt_th = zeros(length(m_p2), length(g_p1_true));                       % Sensing time at SR outage case @ Opt sen time
    
    for i=1:length(m_p2) 
        disp(strcat('mp2 = ',num2str(m_p2(i)))); 
        

        
        for j=1:length(g_p1_true)
            disp(strcat('SNR_p1 = ',num2str(10 * log10( P_p * g_p1_true(j)/noise_power)))); 
            
           %% Expected Rate
            func_C0 = @(t) log2(1 + t * g_s_true * P_s / noise_power) .*...
                1/gamma(m_s(i)) * m_s(i)^(m_s(i)) .* t.^(m_s(i) - 1) .* exp(-m_s(i) * t);
            func_C1 = @(t, u) log2(1 + t * g_s_true * P_s ./ (P_p * u * g_p2_true + noise_power)) .*...
                1/gamma(m_s(i)) * m_s(i)^(m_s(i)) .* t.^(m_s(i) - 1) .* exp(-m_s(i) * t) .*...
                1/gamma(m_p2(i)) * m_p2(i)^(m_p2(i)) .* u.^(m_p2(i) - 1) .* exp(-m_p2(i) * u);
            C_0_id_th = integral(func_C0, 0, 100);
            C_1_id_th = integral2(func_C1, 0, 100, 0, 100);
            [C_0_oc_th, C_1_oc_th] = calc_capacities_fad(P_p, P_s, g_p2_true, g_s_true,...
            noise_power, noise_power/N_s, N_est_th, N_s, m_p2(i), m_s(i));
            parfor k=1:length(N_sen_th)
                %disp(strcat('N_sen_th = ',num2str(N_sen_th(k))));     
                M = 1e2;                                                                % Number of realizations 
                epsilon_id_th = 0;
                epsilon_oc_th = 0;    
                P_rcvd_sim = zeros(1, M);

                if 1 %% This part is executed only once since ideal case is insensitive to 
                    %% Determining the performance meteric R for all the different cases
                    %% Ideal         
                    %% Threshold  
                     if 1       % 1) Determining threshold by considering the outage constraint on the P_d
                        epsilon_id_th = (2/(N_sen_th(k))*(P_p * g_p1_true(j) * gammaincinv(1 - mu, m_p2(i), 'upper')/m_p2(i)  +...
                        noise_power) * gammaincinv(P_d_d, N_sen_th(k)/2, 'upper')); 

                     %% Probability of detection
                        func_exp_P_d = @(t) gammainc(N_sen_th(k) * epsilon_id_th./...
                            (2 *( t *g_p1_true(j)  + noise_power)), N_sen_th(k)/2, 'upper') .*...
                            1/gamma(m_p2(i)) * m_p2(i)^(m_p2(i)) .* t.^(m_p2(i) - 1) .* exp(-m_p2(i) * t);
                        P_d_id_th(i,j,k) = integral(func_exp_P_d, 0, 100);

                     else       % 2) Conventional way: Determining threshold for a certain value of P_d = P_d_d      
                         g_p1_sim = random('gam',m_p2(i), g_p1_true(j)/m_p2(i), 1, M);                       % Power gain for the channel g_p2 (including path loss)
                         for l=1:M
                            P_rcvd_sim(l) = mean(random('norm',...
                                random('norm',0, sqrt(P_p * g_p1_sim(l)), 1, N_sen_th(k)),...
                                sqrt(noise_power), 1, N_sen_th(k)).^2);
                         end
                         test_points = linspace(min(P_rcvd_sim), max(P_rcvd_sim), num_test_points);
                         P_d_id = zeros(1,num_test_points);
                         for l=1:num_test_points
                            func_P_rcvd = @(t) gammainc(N_sen_th(k) * test_points(l)./...
                                (2 * (t * P_p * g_p2_true +  noise_power)), N_sen_th(k)/2, 'upper') .*...
                                1/gamma(m_p2(i)) * m_p2(i)^(m_p2(i)) .* t.^(m_p2(i) - 1) .* exp(- m_p2(i) * t);
                            P_d_id(l) = integral(func_P_rcvd, 0, 100);
                            if P_d_id(l) < P_d_d
                                epsilon_id_th = test_points(l);
                                P_d_id_th(i,j,k) = P_d_id(l);
                                break;
                            end
                         end
                    end 

                   %% Probability of false alarm 
                    P_f_id_th(i,j,k) = gammainc(N_sen_th(k) * epsilon_id_th/(2 * noise_power), N_sen_th(k)/2, 'upper');



                    R_id_th(i,j,k) = (K - N_sen_th(k))/K * (P_H0 * (1 -  P_f_id_th(i,j,k)) * C_0_id_th +...
                        (1 - P_H0) * (1 - P_d_id_th(i,j,k)) * C_1_id_th);                     
                end


                if N_sen_th(k) >= N_est_th
                   N_sen_ni_th = N_sen_th(k);% - N_est_th;

                   %% Outage case
                    if 1
                      %% Threshold  
                        test_points = linspace(0.1 * epsilon_id_th, 1.5 * epsilon_id_th, num_test_points);
                        P_d_oc = zeros(1,num_test_points);
                        for l=1:num_test_points
                           func_P_d = @(t) gammainc(N_est_th * N_sen_ni_th * test_points(l)./...
                               (4 * (P_p * g_p1_true(j) * t + noise_power) * gammaincinv(P_d_d, N_sen_ni_th/2, 'upper')),...
                               N_est_th/2, 'upper') .*...
                               1/gamma(m_p2(i)) * m_p2(i)^(m_p2(i)) .* t.^(m_p2(i) - 1) .* exp(- m_p2(i) * t); 
                           P_d_oc(l) = integral(func_P_d, 0, 100);
                           %P_d_oc(j) = gammainc(N_est_th * N_sen_ni_th * test_points(j)./...
                           %    (4 * (P_p * g_p1_true + noise_power) * gammaincinv(P_d_d, N_sen_ni_th/2, 'upper')),  N_est_th/2, 'upper') .*...
                           %    1; 

                            if P_d_oc(l) <= 1 - mu 
                                epsilon_oc_th = test_points(l);


                              %% Evaluate the intergral using the distribution function
                                if 1
                                    x_pts = linspace(0,1,10000);
                                    CDF_pts = zeros(1,length(x_pts));
                                    for l=1:length(x_pts)
                                        func_P_d = @(t) 1 - gammainc(N_est_th * N_sen_ni_th * epsilon_oc_th./...
                                            (4 * (P_p * g_p1_true(j)  * t + noise_power) * gammaincinv(x_pts(l),...
                                            N_sen_ni_th/2, 'upper')),  N_est_th/2, 'upper') .*...
                                            1/gamma(m_p2(i)) * m_p2(i)^(m_p2(i)) .* t.^(m_p2(i) - 1) .* exp(- m_p2(i) * t); 
                                        CDF_pts(l) = integral(func_P_d, 0, 100);
                                    end
                                    PDF_pts = diff(CDF_pts);
                                    P_d_oc_th(i,j,k) = sum(PDF_pts .* x_pts(2:end));
                                else

                                    func_P_d = @(u, t) u .* exp(gammaln(N_sen_ni_th/2) - gammaln(N_est_th/2) -... 
                                        epsilon_oc_th * N_est_th * N_sen_ni_th./ (4 * (t * P_p * g_p1_true(j)  + noise_power) .*...
                                        gammaincinv(u ,N_sen_ni_th/2,'upper')) + gammaincinv(u, N_sen_ni_th/2,'upper') +...
                                        N_est_th/2 * log(epsilon_oc_th * N_sen_ni_th * N_est_th./ (4 * (t * P_p * g_p1_true(j) + noise_power) .*...
                                        gammaincinv(u, N_sen_ni_th/2,'upper'))) +...
                                        (N_sen_ni_th/2) * log(1./gammaincinv(u, N_sen_ni_th/2,'upper'))) .*...
                                        1/gamma(m_p2(i)) * m_p2(i)^(m_p2(i)) .* t.^(m_p2(i) - 1) .* exp(-m_p2(i) * t);
                                     P_d_oc_th(i,j,k) = integral2(func_P_d, 0, 1, 0, 100);
                                end
                                break;
                            end
                        end


                      %% Probability of false alarm and expected probability of dectection
                        P_f_oc_th(i,j,k) = gammainc(N_sen_ni_th * epsilon_oc_th/(2 * noise_power), N_sen_ni_th/2, 'upper');

                        R_oc_th(i,j,k) = (K - N_sen_th(k))/K * (P_H0 * (1 -  P_f_oc_th(i,j,k)) * C_0_oc_th +...
                        (1 - P_H0) * (1 - P_d_oc_th(i,j,k)) * C_1_oc_th);

                    end   

                end
            end    
            % Optimum throughput for the ideal model
            [R_id_opt_th(i,j) index_id] = max(R_id_th(i,j,:));
            N_id_opt_th(i,j) = N_sen_th(index_id);

            % Optimum throughput for the estimation model -- outage constraint case
            [R_oc_opt_th(i,j) index_oc] = max(R_oc_th(i,j,:));
            N_oc_opt_th(i,j) = N_sen_th(index_oc);
            
            disp(strcat('R_id_opt_th(i,j)  = ',num2str(R_id_opt_th(i,j)))); 
            disp(strcat('N_id_opt_th(i,j)  = ',num2str(N_id_opt_th(i,j)))); 
            disp(strcat('R_oc_opt_th(i,j)  = ',num2str(R_oc_opt_th(i,j)))); 
            disp(strcat('N_oc_opt_th(i,j)  = ',num2str(N_oc_opt_th(i,j))));
            save('results_opt_thr_vs_SNR_fading.mat');
        end
    end
    quit;
end

if 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves theoretical analysis --  Optimum throughput vs SNR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('results_opt_thr_vs_SNR_fading.mat'); 
    Fontsize = 8;
    SNR = (10*log10(P_p * g_p1_true/noise_power));
    index = [1:1:length(SNR)];% [1,5:3:14,15:4:length(N_est_th)];
    figure(1);
        
    
    if 0
        plot(N_sen_th * 1e-3, reshape(R_id_th(1,1,:),1, length(N_sen_th)), 'LineWidth',1); 
    end
    
    %% m = 1
    k = 1;
    h1  = plot(SNR(index), R_id_opt_th(k,index), 'c', 'LineWidth',2);    
    hold on,
    h2 = plot(SNR(index), R_oc_opt_th(k,index), 'LineStyle', '-', 'LineWidth', 1);

    
    %% m = 1
    k = 2;
    hold on,
    plot(SNR(index), R_id_opt_th(k,index), 'c', 'LineWidth',2);    
    hold on,
    plot(SNR(index), R_oc_opt_th(k,index), 'LineStyle', '-', 'LineWidth', 1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves cosmetic makeover
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    grid on;
    axis([min(SNR(index)) max(SNR(index)) 0  max(R_id_opt_th(2, index)) * 1.05]);
    ylabel('$\rs(\test,\ttsen)$','FontSize',Fontsize);
    xlabel('$\snrrcvd$ [dB]','FontSize',Fontsize);
    %hl = legend([h1 h2 h3 h5 h4],'$\rs$', '$\rsac$', '$\rsoc$', '$\ttsen$', 'Sim');
    hl = legend([h1 h2],'IM, Problem 3', 'EM, Problem 4');
    %set(hl, 'position',[0.725 0.12 0.15 0.31]);
    set(hl, 'position',[0.58 0.12 0.32 0.22]);
    set(gca,'FontSize',Fontsize);

    laprint(1, '../figures/fig_opt_thr_vs_SNR_fading', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
end



