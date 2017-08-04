%                                
%        Description: This m-file considers the variation of optimum throughput
%        determined from the Estimation-Throughput tradeoff for the fading channel.  
%        (When an outage constraint is established over the power received at PR)      
%
%        Clearly the regulated power is distributed according non central
%        chi-2 distribution. However to simplify the analysis we
%        approximate the non central distribution against the Gamma
%        distribution.
%        
%        The simulation is performed to show the following:
%        1) Analyse the optimum throughput vs received SNR
%
%        Created on: 04.08.15
%        Revision History: 04.08.15 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

                                              % Enable( = 1) to perform simulation, theoretical analysis
th = 00;                                      % disable ( = 0) to plot curves, data is read from a file
alt = 00;                                     % Theoretical anaylsis is also included of simulation as
                                              % numerical integration is involved
                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_tran = +00;                                        % Power transmitted by PR, the SNR received at ST can be 
                                                     % varried using this parameter
noise_power = -100;                                  % noise power -100 dBm
I_T = -110;                                          % Interference temperature -80 dBm
f_s = 1e6;                                           % 1 MHz one band
K = 0.1 * f_s;                                       % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
%alpha_true_p1 = [120:-0.5:112.5,112.4:-0.01:90];    % (alt, m = 5) True Path loss between ST and PR  to achieve P_p = IT  
alpha_true_p1 = [120:-0.5:114.5,114.4:-0.01:90];    % (alt, m = 1) True Path loss between ST and PR  to achieve P_p = IT  
%alpha_true_p1 = [120:-0.5:90];                       % True Path loss between ST and PR  to achieve P_p = IT  
alpha_true_p2 = 090;
alpha_true_s = 080;                                  % True Path loss between ST and SR  to achieve P_p = IT  
epsilon = [0.10];                                    % Outage Probability constraint at PR


m_p = 1;                                             % Number of realizations of channel for the link ST-PR
m_s = 1;                                             % Number of realizations of channel for the link ST-SR 
P_reg_max = 10.^(-[00 10]/10);                        % Maximum Transmit Power Constraint
P_reg_min = 10^(-30/10);                             % Minimum Transmit Power Constraint

N_s = 10;                                            % Number of Pilot Symobols 

Opt_Exp_R_th = zeros(length(P_reg_max),...
    length(alpha_true_p1));                           % Optimum Expected throughput  
Opt_N_th = zeros(length(P_reg_max),...
    length(alpha_true_p1));                           % Optimum Estimation time  
Opt_P_reg_th = zeros(length(P_reg_max),...
    length(alpha_true_p1));                           % Optimum Controlled Power 
Opt_Exp_R_ga_th = zeros(length(P_reg_max),...
    length(alpha_true_p1));                           % Optimum Expected throughput  
Opt_N_ga_th = zeros(length(P_reg_max),...
    length(alpha_true_p1));                           % Optimum Estimation time  
Opt_P_reg_ga_th = zeros(length(P_reg_max),...
    length(alpha_true_p1));                           % Optimum Controlled Power 

snr = 10^(P_tran/10) * 10.^(-alpha_true_p1/10)...     % snr received at the PR
    / 10^(noise_power/10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th  
    %% Theoretical parameters
    N_th = [10:1:39,40:5:100,110:10:500,520:30:3000];%,3000:100:10000];  % N = Total of samples used for estimation = tau * f_s (ms)
    num_tp = 5000;                             % Number of test points
    P_reg_th = zeros(1, length(N_th));          % Expected power regulated by ST
    Exp_R_th = zeros(1, length(N_th));          % Expected Throughput
    P_reg_ga_th = zeros(1,length(N_th));        % Expected power regulated by ST (applying Gamma Approximation)
    Exp_R_ga_th = zeros(1,length(N_th));        % Expected throughput (applying Gamma Approximation) 
    
    for i = 1:length(P_reg_max)    
       disp(strcat('P_reg_max = ',num2str(P_reg_max(i))));             

       for j = 1:length(snr)
           disp(strcat('snr = ',num2str(10*log10(snr(j)))));      
           parfor k=1:length(N_th)       
               %disp(strcat('N_th = ',num2str(N_th(k))));             

               P_reg_tmp = 0;
               P_reg_ga_tmp = 0;
               if 1 %% Method 2  -- Working But takes to long       
                    P_reg_tp = linspace(P_reg_max(i), P_reg_min, num_tp);
                    for l = 1:1:num_tp                      
                        func_prcvd = @(t) gammainc((10^(I_T/10) * 10^(P_tran/10)/P_reg_tp(l)...
                            + 10^(noise_power/10))./(10^(noise_power/10)/N_th(k)...
                            * (2 + 4 * t * snr(j))./(1 + t * snr(j))),...
                            N_th(k) * (1 + t * snr(j)).^2./(2 + 4 * t * snr(j)))*...
                            1/gamma(m_p) .* t.^(m_p - 1) * m_p^(m_p) .* exp(-m_p * t);
                        result = integral(func_prcvd, 0, 100);                
                        if result >= (1 - epsilon)
                            P_reg_tmp = P_reg_tp(l);                    
                            break;
                        end
                    end
                    P_reg_tmp =  min(P_reg_max(i), P_reg_tmp);
                               

                  %% Expected Rate
                    if 0
                        func_r_em = @(u) (K- N_th(k))/K * log2(1 + u * P_reg_tmp * 10^(-alpha_true_s/10) /...
                            10^(noise_power/10)) *...
                        1/gamma(m_s) * m_s^(m_s) .* u.^(m_s - 1) .* exp(-m_s * u);
                        Exp_R_th(k) = integral(func_r_em, 0, 100);
                    end
                    C_1 = calc_capacities_fad(10^(-alpha_true_p2/10)/10^(noise_power/10),...
                       10^(-alpha_true_s/10), P_reg_tmp, 10^(noise_power/10), N_th(k), N_s, m_p, m_s);
                    Exp_R_th(k) = (K- N_th(k))/K * C_1;
    
                    P_reg_th(k) = P_reg_tmp;
                        
                    %disp(strcat('Exp_R_th(k) = ',num2str(Exp_R_th(k)))); 
                    %disp(strcat('P_reg_th(k) = ',num2str(P_reg_th(k))));
               end

               if 0  %% Method 3 -- Approximating the first two moments of the distribution for Prcvd with 
                     %% the moments of Gamma distribution -- Not so Accurate for low SNR and for N > 1000

                    %% The mean and variance are determined from the moment generating function 
                    mean_nc = 10^(noise_power/10) * (1 + m_p * snr(j));
                    var_nc =  (10^(noise_power/10))^2/N_th(k) * (2 * + 4 * m_p *...
                        snr(j) + N_th(k) * m_s * (snr(j))^2);
                    % Parameters for Gamma distribution 
                    b = var_nc / mean_nc;
                    a = mean_nc / b;
                    P_reg_ga_tmp = 10^(I_T/10) * 10^(P_tran/10) / (b * gammaincinv((1 - epsilon),a) - 10^(noise_power/10));       
                %    P_reg_ga_tmp =  min(P_reg_max(i), P_reg_ga_tmp);
                    P_reg_ga_tmp =  max(P_reg_max(i), P_reg_ga_tmp);

                    
                  %% Expected Rate
                    if 0
                        func_r_em = @(u) (K- N_th(k))/K * log2(1 + u * P_reg_ga_tmp * 10^(-alpha_true_s/10) /...
                            10^(noise_power/10)) *...
                            1/gamma(m_s) * m_s^(m_s) .* u.^(m_s - 1) .* exp(-m_s * u);
                        Exp_R_ga_th(k) = integral(func_r_em, 0, 100);
                    end
                    C_1 = calc_capacities_fad(10^(-alpha_true_p2/10)/10^(noise_power/10),...
                        10^(-alpha_true_s/10), P_reg_ga_tmp, 10^(noise_power/10), N_th(k), N_s, m_p, m_s);
                    Exp_R_ga_th(k) = (K- N_th(k))/K * C_1;
                    P_reg_ga_th(k) = P_reg_ga_tmp;
                    %disp(strcat('Exp_R_ga_th(k) = ',num2str(Exp_R_ga_th(k)))); 
                    %disp(strcat('P_reg_ga_th(k) = ',num2str(P_reg_ga_th(k))));
               end 
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Finding optimum estimation time
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Optimum throughput for the ideal model
            [Opt_Exp_R_th(i,j) index] = max(Exp_R_th);
            Opt_P_reg_th(i,j) = P_reg_th(index);
            Opt_N_th(i,j) = N_th(index);            
            [Opt_Exp_R_ga_th(i,j) index] = max(Exp_R_ga_th);
            Opt_P_reg_ga_th(i,j) = P_reg_ga_th(index);
            Opt_N_ga_th(i,j) = N_th(index);
            
            disp(strcat('Opt_Exp_R_th = ',num2str(Opt_Exp_R_th(i,j)))); 
            disp(strcat('Opt_N_th = ',num2str(Opt_N_th(i,j))));
            disp(strcat('Opt_P_reg_th = ',num2str(Opt_P_reg_th(i,j))));
            disp(strcat('Opt_Exp_R_ga_th = ',num2str(Opt_Exp_R_ga_th(i,j)))); 
            disp(strcat('Opt_N_ga_th = ',num2str(Opt_N_ga_th(i,j))));                     
            disp(strcat('Opt_P_reg_ga_th = ',num2str(Opt_P_reg_ga_th(i,j))));
            save('results_opt_thr_vs_SNR_fading_m1_e10_SI_10_th.mat');
       end
   end
   quit;
end
if alt
    N_alt = 10:10:100000;   
    Exp_opt_R_alt = zeros(length(snr), length(P_reg_max));
    N_alt_opt = zeros(length(snr), length(P_reg_max));
    temp = zeros(1, length(N_alt));
    for j = 1:length(P_reg_max)
        for k = 1:length(snr)
           disp(strcat('snr = ',num2str(10*log10(snr(k)))));      
            parfor l=1:length(N_alt)
               warning off;     
               func_prcvd = @(t) gammainc((10^(I_T/10) * 10^(P_tran/10)/P_reg_max(j)...
                    + 10^(noise_power/10))./(10^(noise_power/10)/N_alt(l)...
                    * (2 + 4 * t * snr(k))./(1 + t * snr(k))),...
                    N_alt(l) * (1 + t * snr(k)).^2./(2 + 4 * t * snr(k)), 'upper')*...
                    1/gamma(m_p) .* t.^(m_p - 1) * m_p^(m_p) .* exp(-m_p * t);

              %% Determine the SNR that achieves the outage probability
                temp(l) = integral(func_prcvd, 0, 100);
            end                          

           [value index] = min(abs(temp - epsilon));  
           index;
           N_alt_opt(k,j) = N_alt(index); 
           if index ~= length(N_alt) & index ~= 1
               C_1  = calc_capacities_fad(10^(-alpha_true_p2/10)/10^(noise_power/10),...
                        10^(-alpha_true_s/10), P_reg_max(j), 10^(noise_power/10), N_alt(index), N_s, m_p, m_s);
               Exp_opt_R_alt(k,j) = (K - N_alt(index))/K * C_1;
           else
               break;
           end
        disp(strcat('N_alt_opt(k,j) = ',num2str(N_alt_opt(k,j))));
        disp(strcat('Exp_opt_R_alt(k,j) = ',num2str(Exp_opt_R_alt(k,j))));
        save('results_opt_thr_vs_SNR_fading_m1_e10_SI_10_alt.mat');
        end
    end
    quit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Secondary Interference (SI) = 0 dBm 
ver_paper = 1;       % = 1, Corresponds to the changes according to TCCN paper          
                     % = 0, Correspoonds to the changes as per the diss
if ver_paper == 1
    ideal_linewidth = 1;
else
    ideal_linewidth = 2;
end 
                    
if 1
    %% m = 5
    load('results_opt_thr_vs_SNR_fading_m5_e10_SI_00_th.mat');
    marker_index = 1:5:length(snr);
    SNR = 10*log10(snr);
    

    %pbaspect([1 0.3 1]);
    Fontsize = 8.5;
    %% Ideal Model
    %% Outage constraint on the contolled power
    Exp_R_id = zeros(length(epsilon), length(SNR));
    P_reg_id =  zeros(length(epsilon), length(SNR));
    for i = 1:length(epsilon) 
        for j = 1:length(SNR)

               %P_reg_id = 10^(I_T/10) / 10^(-alpha_true_p(j)/10) / log(1/(epsilon(i)));                
               P_reg_id(j) = 10^(I_T/10) * m_p / 10^(-alpha_true_p1(j)/10) / gammaincinv(epsilon(i), m_p, 'upper');
               %P_reg_id =  min(P_reg_max, P_reg_id);
               P_reg_id(j) =  max(P_reg_min, P_reg_id(j));
               func_r_id = @(u, t) log2(1 + P_reg_id(j) * u *  10^(-alpha_true_s/10)...
                    ./ (10^(noise_power/10)  + t * 10^(-alpha_true_p2/10))) .*...
                    1/gamma(m_s) * m_s^(m_s) .* u.^(m_s - 1) .* exp(-m_s * u) .*...
                    1/gamma(m_p) * m_p^(m_p) .* t.^(m_p - 1) .* exp(-m_p * t);
               Exp_R_id(i,j) = integral2(func_r_id, 0, 100, 0, 100); 
        end
    end
    h1 = plot(SNR, Exp_R_id(1,:), 'c', 'LineWidth', ideal_linewidth);
    hold on,
    h2 = plot(SNR, Opt_Exp_R_th(1,:), 'k', 'LineWidth', 1);
    if ver_paper
        hold on,
        gap = 12;
        downsampled = 1:ceil(length(SNR)/gap):length(SNR);
        plot(SNR(downsampled), Exp_R_id(1, downsampled), 'cd', 'LineWidth', 1);
        lp1 = plot(0,0,'-c','Marker','d', 'LineWidth', 1, 'visible','off');
    end

    if ver_paper            
        hold on,
        plot(SNR(downsampled), Opt_Exp_R_th(1,downsampled), 'ko', 'LineWidth', 1);
        lp2 = plot(0,0,'-k','Marker','o', 'LineWidth', 1, 'visible','off');        
    end

    %hold on,
    %plot(SNR, Opt_Exp_R_th(2,:), 'k', 'LineWidth', 1);
    if 0
        hold on,
        plot(SNR(marker_index), Exp_R_id(1,marker_index), 's', 'LineWidth', 1);
        hold on,
        plot(SNR(marker_index), Opt_Exp_R_th(1,marker_index), 's', 'LineWidth', 1);
        hold on,
        plot(SNR(marker_index), Opt_Exp_R_th(2,marker_index), 's', 'LineWidth', 1);
    end
    %hold on,
    %h3 = plot(SNR, Opt_Exp_R_ga_th(1,:), 'k--', 'LineWidth', 1);
    Ymax = max(Exp_R_id(1,:));
    
    if 01
        %% Alternative Strategy --- With no power control
        load('results_opt_thr_vs_SNR_fading_m5_e10_SI_00_alt.mat');
        hold on,
        SNR_alt = (10*log10(snr))';
        h3  = plot(SNR_alt, Exp_opt_R_alt(:,1), 'r--', 'LineWidth', 1);   
        if ver_paper
            downsampled_alt = [];
            for i = downsampled
                downsampled_alt(end + 1) = find(SNR_alt == SNR(i));
            end
            hold on,
            plot(SNR_alt(downsampled_alt), Exp_opt_R_alt(downsampled_alt,1), 'rs', 'LineWidth', 1);
            lp3 = plot(0,0,'--r','Marker','s', 'LineWidth', 1, 'visible','off');
        end
        %hold on,
        %plot(SNR_alt, Exp_opt_R_alt(:,2), 'r--', 'LineWidth', 1);
    end
end

if 1

    %% m = 1
    load('results_opt_thr_vs_SNR_fading_m1_e10_SI_00_th.mat');
    SNR = 10*log10(snr);
    %% Ideal Model
    %% Outage constraint on the contolled power
    Exp_R_id = zeros(length(epsilon), length(SNR));
    P_reg_id =  zeros(length(epsilon), length(SNR));   
    for i = 1:length(epsilon) 
        for j = 1:length(SNR)

               %P_reg_id = 10^(I_T/10) / 10^(-alpha_true_p(j)/10) / log(1/(epsilon(i)));                
               P_reg_id(j) = 10^(I_T/10) * m_p / 10^(-alpha_true_p1(j)/10) / gammaincinv(epsilon(i), m_p, 'upper');
               %P_reg_id =  min(P_reg_max, P_reg_id);
               P_reg_id(j) =  max(P_reg_min, P_reg_id(j));
               func_r_id = @(u, t) log2(1 + P_reg_id(j) * u *  10^(-alpha_true_s/10)...
                    ./ (10^(noise_power/10)  + t * 10^(-alpha_true_p2/10))) .*...
                    1/gamma(m_s) * m_s^(m_s) .* u.^(m_s - 1) .* exp(-m_s * u) .*...
                    1/gamma(m_p) * m_p^(m_p) .* t.^(m_p - 1) .* exp(-m_p * t);
               Exp_R_id(i,j) = integral2(func_r_id, 0, 100, 0, 100); 
        end
    end
    hold on,
    plot(SNR, Exp_R_id(1,:), 'c', 'LineWidth', ideal_linewidth);
    hold on,
    plot(SNR, Opt_Exp_R_th(1,:), 'k', 'LineWidth', 1);
    if ver_paper
        hold on,
        plot(SNR(downsampled), Exp_R_id(1, downsampled), 'cd', 'LineWidth', 1);
        hold on,
        plot(SNR(downsampled), Opt_Exp_R_th(1,downsampled), 'ko', 'LineWidth', 1);
    end
    %hold on,
    %plot(SNR, Opt_Exp_R_th(2,:), 'k', 'LineWidth', 1);
    if 0
        hold on,
        plot(SNR(marker_index), Exp_R_id(1,marker_index), 's', 'LineWidth', 1);
        hold on,
        plot(SNR(marker_index), Opt_Exp_R_th(1,marker_index), 's', 'LineWidth', 1);
        hold on,
        plot(SNR(marker_index), Opt_Exp_R_th(2,marker_index), 's', 'LineWidth', 1);
    end
    %hold on,
    %plot(SNR, Opt_Exp_R_ga_th(1,:), 'k--', 'LineWidth', 1);

end

if 01
    %% Alternative Strategy --- With no power control
    load('results_opt_thr_vs_SNR_fading_m1_e10_SI_00_alt.mat');
    hold on,
    SNR_alt = (10*log10(snr))';
    plot(SNR_alt, Exp_opt_R_alt(:,1), 'r--', 'LineWidth', 1);
    if ver_paper
        downsampled_alt = [];
        for i = downsampled
            downsampled_alt(end + 1) = find(SNR_alt == SNR(i));
        end
        hold on,
        plot(SNR_alt(downsampled_alt), Exp_opt_R_alt(downsampled_alt,1), 'rs', 'LineWidth', 1);    
    end
    %hold on,
    %plot(SNR_alt, Exp_opt_R_alt(:,2), 'r--', 'LineWidth', 1);
end


grid on;
axis([min(SNR) max(SNR) 0 8.5 * 1.01]);
ylabel('$\rs(\ttau)$ [bits/sec/Hz]','FontSize',Fontsize);
xlabel('$\gamma$ [dB]','FontSize',Fontsize);
if ver_paper
    hl = legend([lp1 lp2 lp3],'IM', 'EM', 'Coro 3');
else
    hl = legend([h1 h2 h3],'IM', 'EM', 'Coro 3');
end
%set(hl, 'position',[0.725 0.12 0.15 0.31]);
set(hl, 'position',[0.66 0.72 0.24 0.2]);
set(gca,'FontSize',Fontsize);
 
if 0
    %%%% Zommed Image %%%%%
    % create a new pair of axes inside current figure
    ax = axes('Position',[.52 .52 .32 .32]);
    box on; % put box around new pair of axes
    child_gca = get(gca,'Children');

    [value indexOfInterest] = find(SNR <= 1 & SNR >= 0);
    index = 2;
    plot(SNR(indexOfInterest), Exp_opt_R_th(index,indexOfInterest), 'k', 'LineWidth', 1); % plot on new axes
    hold on,
    plot(SNR(indexOfInterest), Exp_opt_R_ga_th(index,indexOfInterest), 'k--', 'LineWidth', 1); 
    hold on,
    plot(SNR(indexOfInterest), Exp_R_id(index,indexOfInterest), 'c', 'LineWidth', 1); 
    ylabel(gca, '$\rs(\ttau)$','FontSize',Fontsize);
    xlabel(gca, '$\gamma$','FontSize',Fontsize);
    set(ax, 'XTick', [0 0.5 1]);
    set(ax, 'YTick', [1.8 1.9 2.0]);
    axis(ax, [min(SNR(indexOfInterest)) max(SNR(indexOfInterest))...
        min(Exp_opt_R_ga_th(index,(indexOfInterest)))  max(Exp_R_id(index, indexOfInterest)) * 1]);
    th = title('Zoom');
    set(ht, 'Position', [0.5 0.2]);
end

if 01
    if ver_paper
        laprint(1, '../figures/fig_opt_thr_vs_SNR_SI_00_fading', 'options', 'factory', 'width', 8, 'scalefonts',...
            'on', 'factor',0.5, 'keepfontprops', 'on');
    else
        laprint(1, '../figures/fig_opt_thr_vs_SNR_SI_00_fading_diss', 'options', 'factory', 'width', 8, 'scalefonts',...
            'on', 'factor',0.5, 'keepfontprops', 'on');
    end
end

%% Secondary Interference (SI) = -10 dBm 
figure(2);
if 1
    %% m = 5
    load('results_opt_thr_vs_SNR_fading_m5_e10_SI_10_th.mat');
    marker_index = 1:5:length(snr);
    SNR = 10*log10(snr);

    %pbaspect([1 0.3 1]);
    Fontsize = 8.5;
    %% Ideal Model
    %% Outage constraint on the contolled power
    Exp_R_id = zeros(length(epsilon), length(SNR));
    P_reg_id =  zeros(length(epsilon), length(SNR));
    for i = 1:length(epsilon) 
        for j = 1:length(SNR)

               %P_reg_id = 10^(I_T/10) / 10^(-alpha_true_p(j)/10) / log(1/(epsilon(i)));                
               P_reg_id(j) = 10^(I_T/10) * m_p / 10^(-alpha_true_p1(j)/10) / gammaincinv(epsilon(i), m_p, 'upper');
               %P_reg_id =  min(P_reg_max, P_reg_id);
               P_reg_id(j) =  max(P_reg_min, P_reg_id(j));
               func_r_id = @(u, t) log2(1 + P_reg_id(j) * u *  10^(-alpha_true_s/10)...
                    ./ (10^(noise_power/10)  + t * 10^(-alpha_true_p2/10))) .*...
                    1/gamma(m_s) * m_s^(m_s) .* u.^(m_s - 1) .* exp(-m_s * u) .*...
                    1/gamma(m_p) * m_p^(m_p) .* t.^(m_p - 1) .* exp(-m_p * t);
               Exp_R_id(i,j) = integral2(func_r_id, 0, 100, 0, 100); 
        end
    end
    h1 = plot(SNR, Exp_R_id(1,:), 'c', 'LineWidth', ideal_linewidth);
    hold on,
    if ver_paper
        plot(SNR(downsampled), Exp_R_id(1, downsampled), 'cd', 'LineWidth', 1);
        lp1 = plot(0,0,'-c','Marker','d', 'LineWidth', 1, 'visible','off');
    end
    hold on,
    h2 = plot(SNR, Opt_Exp_R_th(1,:), 'k', 'LineWidth', 1);
    hold on,
    if ver_paper
        plot(SNR(downsampled), Opt_Exp_R_th(1,downsampled), 'ko', 'LineWidth', 1);
        lp2 = plot(0,0,'-k','Marker','o', 'LineWidth', 1, 'visible','off');
    end

    %hold on,
    %plot(SNR, Opt_Exp_R_th(2,:), 'k', 'LineWidth', 1);
    %hold on,
    %h3 = plot(SNR, Opt_Exp_R_ga_th(1,:), 'k--', 'LineWidth', 1);
    Ymax = max(Exp_R_id(1,:));
    if 0
        hold on,
        plot(SNR(marker_index), Exp_R_id(1,marker_index), 's', 'LineWidth', 1);
        hold on,
        plot(SNR(marker_index), Opt_Exp_R_th(1,marker_index), 's', 'LineWidth', 1);
        hold on,
        plot(SNR(marker_index), Opt_Exp_R_th(2,marker_index), 's', 'LineWidth', 1);
        lp3 = plot(0,0,'--r','Marker','s', 'LineWidth', 1, 'visible','off');
    end
end

if 01
    %% Alternative Strategy --- With no power control
    load('results_opt_thr_vs_SNR_fading_m5_e10_SI_10_alt.mat');
    hold on,
    SNR_alt = (10*log10(snr))';
    h3  = plot(SNR_alt, Exp_opt_R_alt(:,1), 'r--', 'LineWidth', 1);
    if ver_paper
        downsampled_alt = [];
        for i = downsampled
            downsampled_alt(end + 1) = find(SNR_alt == SNR(i));
        end
        hold on,
        plot(SNR_alt(downsampled_alt), Exp_opt_R_alt(downsampled_alt,1), 'rs', 'LineWidth', 1); 
    end
    %hold on,
    %plot(SNR_alt, Exp_opt_R_alt(:,2), 'r--', 'LineWidth', 1);
end

if 1

    %% m = 1
    load('results_opt_thr_vs_SNR_fading_m1_e10_SI_10_th.mat');
    SNR = 10*log10(snr);
    %% Ideal Model
    %% Outage constraint on the contolled power
    Exp_R_id = zeros(length(epsilon), length(SNR));
    P_reg_id =  zeros(length(epsilon), length(SNR));   
    for i = 1:length(epsilon) 
        for j = 1:length(SNR)

               %P_reg_id = 10^(I_T/10) / 10^(-alpha_true_p(j)/10) / log(1/(epsilon(i)));                
               P_reg_id(j) = 10^(I_T/10) * m_p / 10^(-alpha_true_p1(j)/10) / gammaincinv(epsilon(i), m_p, 'upper');
               %P_reg_id =  min(P_reg_max, P_reg_id);
               P_reg_id(j) =  max(P_reg_min, P_reg_id(j));
               func_r_id = @(u, t) log2(1 + P_reg_id(j) * u *  10^(-alpha_true_s/10)...
                    ./ (10^(noise_power/10)  + t * 10^(-alpha_true_p2/10))) .*...
                    1/gamma(m_s) * m_s^(m_s) .* u.^(m_s - 1) .* exp(-m_s * u) .*...
                    1/gamma(m_p) * m_p^(m_p) .* t.^(m_p - 1) .* exp(-m_p * t);
               Exp_R_id(i,j) = integral2(func_r_id, 0, 100, 0, 100); 
        end
    end
    hold on,
    plot(SNR, Exp_R_id(1,:), 'c', 'LineWidth', ideal_linewidth);
    if ver_paper
        hold on,
        plot(SNR, Opt_Exp_R_th(1,:), 'k', 'LineWidth', 1);
        hold on,
        plot(SNR(downsampled), Exp_R_id(1, downsampled), 'cd', 'LineWidth', 1);
        hold on,
        plot(SNR(downsampled), Opt_Exp_R_th(1,downsampled), 'ko', 'LineWidth', 1);
    end
    %hold on,
    %plot(SNR, Opt_Exp_R_th(2,:), 'k', 'LineWidth', 1);
    if 0
        hold on,
        plot(SNR(marker_index), Exp_R_id(1,marker_index), 's', 'LineWidth', 1);
        hold on,
        plot(SNR(marker_index), Opt_Exp_R_th(1,marker_index), 's', 'LineWidth', 1);
        hold on,
        plot(SNR(marker_index), Opt_Exp_R_th(2,marker_index), 's', 'LineWidth', 1);
    end
    %hold on,
    %plot(SNR, Opt_Exp_R_ga_th(1,:), 'k--', 'LineWidth', 1);

end

if 01
    %% Alternative Strategy --- With no power control
    load('results_opt_thr_vs_SNR_fading_m1_e10_SI_10_alt.mat');
    hold on,
    SNR_alt = (10*log10(snr))';
    h3  = plot(SNR_alt, Exp_opt_R_alt(:,1), 'r--', 'LineWidth', 1);  
    if ver_paper
        downsampled_alt = [];
        for i = downsampled
            downsampled_alt(end + 1) = find(SNR_alt == SNR(i));
        end
        hold on,
        plot(SNR_alt(downsampled_alt), Exp_opt_R_alt(downsampled_alt,1), 'rs', 'LineWidth', 1); 
    end
    %hold on,
    %plot(SNR_alt, Exp_opt_R_alt(:,2), 'r--', 'LineWidth', 1);
end

grid on;
axis([min(SNR) max(SNR) 0 6 * 1.01]);
ylabel('$\rs(\ttau)$ [bits/sec/Hz]','FontSize',Fontsize);
xlabel('$\gamma$ [dB]','FontSize',Fontsize);
if ver_paper
    hl = legend([lp1 lp2 lp3],'IM', 'EM', 'Coro 3');
else
    hl = legend([h1 h2 h3],'IM', 'EM', 'Coro 3');
end
%set(hl, 'position',[0.725 0.12 0.15 0.31]);
set(hl, 'position',[0.66 0.72 0.24 0.2]);
set(gca,'FontSize',Fontsize);

if 0
    %%%% Zommed Image %%%%%
    % create a new pair of axes inside current figure
    ax = axes('Position',[.52 .52 .32 .32]);
    box on; % put box around new pair of axes
    child_gca = get(gca,'Children');

    [value indexOfInterest] = find(SNR <= 1 & SNR >= 0);
    index = 2;
    plot(SNR(indexOfInterest), Exp_opt_R_th(index,indexOfInterest), 'k', 'LineWidth', 1); % plot on new axes
    hold on,
    plot(SNR(indexOfInterest), Exp_opt_R_ga_th(index,indexOfInterest), 'k--', 'LineWidth', 1); 
    hold on,
    plot(SNR(indexOfInterest), Exp_R_id(index,indexOfInterest), 'c', 'LineWidth', 1); 
    ylabel(gca, '$\rs(\ttau)$','FontSize',Fontsize);
    xlabel(gca, '$\gamma$','FontSize',Fontsize);
    set(ax, 'XTick', [0 0.5 1]);
    set(ax, 'YTick', [1.8 1.9 2.0]);
    axis(ax, [min(SNR(indexOfInterest)) max(SNR(indexOfInterest))...
        min(Exp_opt_R_ga_th(index,(indexOfInterest)))  max(Exp_R_id(index, indexOfInterest)) * 1]);
    th = title('Zoom');
    set(ht, 'Position', [0.5 0.2]);
end
if 1
    if ver_paper
        laprint(2, '../figures/fig_opt_thr_vs_SNR_SI_10_fading', 'options', 'factory', 'width', 8, 'scalefonts',...
            'on', 'factor',0.5, 'keepfontprops', 'on');
    else
        laprint(2, '../figures/fig_opt_thr_vs_SNR_SI_10_fading_diss', 'options', 'factory', 'width', 8, 'scalefonts',...
            'on', 'factor',0.5, 'keepfontprops', 'on');
    end
end

if 0 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Optium estimatio time vs SNR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2)
    %load('results_opt_thr_vs_SNR_fading_th.mat');
    points = [1:2:length(snr),26];
    SNR = 10*log10(snr);
    h2 = plot(SNR(points), Opt_N_th(1,(points)) * 1e-3, 'k', 'LineWidth', 1);
    hold on,
    h3 = plot(SNR(points), Opt_N_ga_th(1,points) * 1e-3, 'k--', 'LineWidth', 1);
    hold on,
    h4 = plot(SNR(points), Opt_N_th(2,points) * 1e-3 , 'k', 'LineWidth', 1);
    hold on,
    h5 = plot(SNR(points), Opt_N_ga_th(2,points) * 1e-3, 'k--', 'LineWidth', 1);
    %hold on,
    %h5 = plot(SNR, Opt_N_th(3,:) * 1e-3, 'k', 'LineWidth', 1);
    grid on;
    axis([min(SNR) max(SNR) 0  max(Opt_N_th(1,(points))) * 1e-3]);
    ylabel('$\ttau$ [ms]','FontSize',Fontsize);
    xlabel('$\gamma$ [dB]','FontSize',Fontsize);
    hl = legend([h1 h2 h3],'IM', 'EM', 'Lemma');
    %set(hl, 'position',[0.725 0.12 0.15 0.31]);
    set(hl, 'position',[0.68 0.71 0.22 0.2]);
    set(gca,'FontSize',Fontsize);
    laprint(2, '../figures/fig_opt_est_time_vs_SNR_fading', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
end
