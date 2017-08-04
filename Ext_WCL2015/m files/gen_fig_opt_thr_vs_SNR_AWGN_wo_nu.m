%                                
%        Description: This m-file considers the variation of optimum throughput
%        determined from the Estimation-Throughput tradeoff for the AWGN channel.  
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
                                              % Theoretical anaylsis is also included of simulation as
                                              % numerical integration is involved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_tran = -00;                                     % Power transmitted by PR, the SNR received at ST can be 
                                                  % varried using this parameter
noise_power = -100;                               % noise power -100 dBm
I_T = -110;                                       % Interference temperature -80 dBm
f_s = 1e6;                                        % 1 MHz one band
K = 0.1 * f_s;                                    % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
alpha_true_p1 = [120:-0.5:90];                    % True Path loss between ST and PR  to achieve P_p = IT  
alpha_true_p2 = 100;                              % True Path loss between PT and SR  
alpha_true_s = 080;                               % True Path loss between ST and SR  to achieve P_p = IT  
epsilon = [0.01, 0.1];                            % Outage Probability constraint at PR

P_reg_max = 10.^([00 -10]/10);                   % Maximum Trasmit Power Constraint


Exp_opt_R_th = zeros(length(epsilon),...
    length(alpha_true_p1), length(P_reg_max));     % Optimum Expected throughput  

Exp_R_id = zeros(1, length(alpha_true_p1));        % Idle throughput  

Opt_N_th = zeros(length(epsilon),...
    length(alpha_true_p1), length(P_reg_max));     % Optimum Estimation time 
snr = 10^(P_tran/10) * 10.^(-alpha_true_p1/10)...  % snr received at the PR
    / 10^(noise_power/10);

N_s = 10;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th  
    %% Theoretical parameters
    N_th = 10:1:7000;                               % N = Total of samples used for estimation = tau * f_s 
    P_reg_th = zeros(1, length(N_th));              % Expected power regulated by ST
    Exp_R_th = zeros(1, length(N_th));              % Expected Throughput
    for i = 1:length(epsilon)
       disp(strcat('epsilon = ',num2str(epsilon(i))));             

       for j = 1:length(snr)
           disp(strcat('snr = ',num2str(10*log10(snr(j)))));  
           for k = 1:length(P_reg_max)
               parfor l=1:length(N_th)       
                   %disp(strcat('N = ',num2str(N_th(k))));             

                   %% Determining the performance meterics --> Exp_R
                   %% Expected values           

                   %% Gamma Approximation to the non-central chi-squared distribution
                   mean = 10^(noise_power/10) * (1 + snr(j));       
                   var = (10^(noise_power/10))^2/N_th(l) * ( 2 + 4 * snr(j)); 

                   b = var/mean;
                   a = mean/b;

                   %% Determine the controlled power
                   P_reg_th(l) = min(P_reg_max(k), 10^(I_T/10) * 10^(P_tran/10) /...
                       (b * gammaincinv((1 - epsilon(i)),a) - 10^(noise_power/10)));   
                   
                   C_1  = calc_capacities(10^(-alpha_true_p2/10)/10^(noise_power/10),...
                        10^(-alpha_true_s/10), P_reg_th(l), 10^(noise_power/10), N_th(l), N_s);

                   %% Expected Rate
                    Exp_R_th(l) = (K- N_th(l))/K * C_1;
               end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   Finding optimum estimation time
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Optimum throughput for the ideal model
                [Exp_opt_R_th(i,j,k) index] = max(Exp_R_th);
                Opt_N_th(i,j, k) = N_th(index);
                
                % Idle Case
                Exp_R_id(j,k) = log2(1 + 10^(I_T/10) ./(10^(P_tran/10) * 10.^(-alpha_true_p1(j)/10))...
                   * 10^(-alpha_true_s/10)./ (10^(P_tran/10) *  (10^(-alpha_true_p2/10) + 10^(noise_power/10))));
                
                disp(strcat('Exp_R_id = ',num2str(Exp_R_id(j,k)))); 
                disp(strcat('Exp_R_th = ',num2str(Exp_opt_R_th(i,j, k)))); 
                disp(strcat('N_opt_th = ',num2str(Opt_N_th(i,j, k))));
           end
           save('results_opt_thr_vs_SNR_AWGN_e01_SI_00_th.mat');
       end
   end
   quit;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ver_paper = 1;       % = 1, Corresponds to the changes according to TCCN paper          
                     % = 0, Correspoonds to the changes as per the diss
if ver_paper == 1
    ideal_linewidth = 1;
else
    ideal_linewidth = 2;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Optium Throughput vs SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
load('results_opt_thr_vs_SNR_AWGN_e01_SI_00_th.mat');
SNR = 10*log10(snr);

%% Ideal Model
Exp_R_id = log2(1 + 10^(I_T/10) ./(10^(P_tran/10) * 10.^(-alpha_true_p1/10))...
    * 10^(-alpha_true_s/10)./ (10^(P_tran/10) *  (10^(-alpha_true_p2/10) + 10^(noise_power/10))));
h11 = plot(SNR, Exp_R_id, 'c', 'LineWidth', ideal_linewidth);
hold on,
if ver_paper
    gap = 12;
    downsampled = 1:ceil(length(SNR)/gap):length(SNR);
    plot(SNR(downsampled), Exp_R_id(downsampled), 'cd', 'LineWidth', 1);
    lp1 = plot(0,0,'-c','Marker','d', 'LineWidth', 1, 'visible','off');
end

%pbaspect([1 0.3 1]);
Fontsize = 8.5;
hold on,
%h2 = plot(SNR, reshape(Exp_opt_R_th(1,:,1), 1, length(SNR)), 'k-.', 'LineWidth', 1);
%hold on,
plot(SNR, reshape(Exp_opt_R_th(2,:,1), 1, length(SNR)), 'k-', 'LineWidth', 1);
hold on,
%plot(SNR, reshape(Exp_opt_R_th(1,:,2), 1, length(SNR)), 'k-.', 'LineWidth', 1);
%hold on,
h33 = plot(SNR, reshape(Exp_opt_R_th(2,:,2), 1, length(SNR)), 'k-', 'LineWidth', 1);
hold on,
if ver_paper
    hold on,
    plot(SNR(downsampled), reshape(Exp_opt_R_th(2,downsampled,1), 1, length(downsampled)), 'ko', 'LineWidth', 1);

    plot(SNR(downsampled), reshape(Exp_opt_R_th(2,downsampled,2), 1, length(downsampled)), 'ko', 'LineWidth', 1);
    lp2 = plot(0,0,'-k','Marker','o', 'LineWidth', 1, 'visible','off');
end
%hold on,
%h4 = plot(SNR, Exp_opt_R_th(3,:), 'k', 'LineWidth', 1);


alt = 00;
if alt
    %% Alternative Strategy --- With no power control
    N_alt = 10:10:100000;   
    alpha_true_p1 = [120:-0.05:90];                    % True Path loss between ST and PR  to achieve P_p = IT  
    snr = 10^(P_tran/10) * 10.^(-alpha_true_p1/10)...  % snr received at the PR
        / 10^(noise_power/10);
    Exp_opt_R_alt = zeros(length(epsilon), length(snr), length(P_reg_max));
    N_alt_opt = zeros(length(epsilon), length(snr), length(P_reg_max));

    for i = 1:length(epsilon)
        for j = 1:length(P_reg_max)
            for k = 1:length(snr)


                %% Gamma Approximation to the non-central chi-squared distribution
                mean = 10^(noise_power/10) * (1 + snr(k));       
                var = (10^(noise_power/10))^2./N_alt * ( 2 + 4 * snr(k)); 

                b = var/mean;
                a = mean./b;

                %% Determine the SNR that achieves the outage probability
                temp1 = gammainc((10^(I_T/10) * 10^(P_tran/10)...
                   /P_reg_max(j) + 10^(noise_power/10)) * 1./b, a, 'upper');
               [value index] = min(abs(temp1 - epsilon(i)));  
               index
               N_alt_opt(i,k,j) = index; 
               if (index~=length(N_alt)) && value < epsilon(i)
                   C_1  = calc_capacities(10^(-alpha_true_p2/10)/10^(noise_power/10),...
                            10^(-alpha_true_s/10), P_reg_max(j), 10^(noise_power/10), N_alt(index), N_s)
                   Exp_opt_R_alt(i,k,j) = (K - N_alt(index))/K * C_1;
               else
                   break;
               end
            end
        end
    end
    save('results_opt_thr_vs_SNR_AWGN_SI_00_alt.mat');
end
load('results_opt_thr_vs_SNR_AWGN_SI_00_alt.mat');
SNR_alt = 10*log10(snr);

%hold on,
%plot(SNR, reshape(Exp_opt_R_alt(1,:,1), 1, length(SNR)), 'r--', 'LineWidth', 1);
hold on,
h44  = plot(SNR_alt, reshape(Exp_opt_R_alt(2,:,1), 1, length(SNR_alt)), 'r--', 'LineWidth', 1);

hold on,
plot(SNR_alt, reshape(Exp_opt_R_alt(2,:,2), 1, length(SNR_alt)), 'r--', 'LineWidth', 1);


if ver_paper
    downsampled_alt = [];
    for i = downsampled
        downsampled_alt(end + 1) = find(SNR_alt == SNR(i));
    end
    hold on,
    plot(SNR_alt(downsampled_alt), reshape(Exp_opt_R_alt(2,downsampled_alt,1), 1, length(downsampled_alt)), 'rs', 'LineWidth', 1);

    %hold on,
    %plot(SNR, reshape(Exp_opt_R_alt(1,:,2), 1, length(SNR)), 'r-', 'LineWidth', 1);
    hold on,
    plot(SNR_alt(downsampled_alt), reshape(Exp_opt_R_alt(2,downsampled_alt,2), 1, length(downsampled_alt)), 'rs', 'LineWidth', 1);
    lp3 = plot(0,0,'--r','Marker','s', 'LineWidth', 1, 'visible','off');
end

grid on;
axis([min(SNR) max(SNR) 0  max(Exp_R_id) * 1.01]);
ylabel('$\rs(\ttau)$ [bits/sec/Hz]','FontSize',Fontsize);
xlabel('$\gamma$ [dB]','FontSize',Fontsize);
%hl = legend([h1 h3 h4],'IM', 'EM, $\opc = 0.1$', 'EM, $\opc = 0.01$');
if ver_paper
    hl = legend([lp1 lp2 lp3],'IM', 'EM','Coro 2');
else
hl = legend([h11 h33 h44],'IM', 'EM','Coro 2');
end
set(hl, 'position',[0.66 0.74 0.24 0.18]);
set(gca,'FontSize',Fontsize);

if 00
    %%%% Zommed Image %%%%%
    % create a new pair of axes inside current figure
    ax = axes('Position',[.56 .51 .33 .33]);
    box on; % put box around new pair of axes
    child_gca = get(gca,'Children');

    [value indexOfInterest] = find(SNR <= -0.0 & SNR >= -1.0);
    %plot(SNR(indexOfInterest), reshape(Exp_opt_R_th(1,indexOfInterest,1), 1, length(SNR(indexOfInterest))),...
    %'k-.', 'LineWidth', 1); % plot on new axes
    %hold on,
    plot(SNR(indexOfInterest), reshape(Exp_opt_R_th(2,indexOfInterest,1), 1, length(SNR(indexOfInterest))),...
    'k-', 'LineWidth', 1); 
    hold on,
    %plot(SNR(indexOfInterest), reshape(Exp_opt_R_th(1,indexOfInterest,2), 1, length(SNR(indexOfInterest))),...
    %'k-.', 'LineWidth', 1); 
    %hold on,
    plot(SNR(indexOfInterest), reshape(Exp_opt_R_th(2,indexOfInterest,2), 1, length(SNR(indexOfInterest))),...
    'k-', 'LineWidth', 1); 
    plot(SNR(indexOfInterest), Exp_R_id(indexOfInterest), 'c', 'LineWidth', 3);


    % plot on new axes
    %plot(SNR(indexOfInterest), Exp_opt_R_th(3,indexOfInterest), 'k', 'LineWidth', 1); 
    %ylabel(gca, '$\rs(\ttau)$','FontSize',Fontsize);
    %xlabel(gca, '$\gamma$','FontSize',Fontsize);
    %set(ax, 'XTick', [-1.3 -0.9 -0.5]);
    %set(ax, 'YTick', [3.2 3.4 3.7]);
    axis(ax, [min(SNR(indexOfInterest)) max(SNR(indexOfInterest))...
        min(reshape(Exp_opt_R_th(1,indexOfInterest,1), 1, length(SNR(indexOfInterest))))...
        max(Exp_R_id(indexOfInterest)) * 1]);
    set(gca,'FontSize',Fontsize);
    th = title(ax, 'Zoom');
    set(th, 'Position', [-0.85    3.4    1.0001]);
end
if ver_paper    
    laprint(1, '../figures/fig_opt_thr_vs_SNR_AWGN_SI_00', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
else    
    laprint(1, '../figures/fig_opt_thr_vs_SNR_AWGN_SI_00_diss', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
end

%% SI = -10 dBm
figure(2);
load('results_opt_thr_vs_SNR_AWGN_e01_SI_10_th.mat');
SNR = 10*log10(snr);

%% Ideal Model
Exp_R_id = log2(1 + 10^(I_T/10) ./(10^(P_tran/10) * 10.^(-alpha_true_p1/10))...
    * 10^(-alpha_true_s/10)./ (10^(P_tran/10) *  (10^(-alpha_true_p2/10) + 10^(noise_power/10))));
h11 = plot(SNR, Exp_R_id, 'c', 'LineWidth', ideal_linewidth);
if ver_paper
    hold on,
    gap = 12;
    downsampled = 1:ceil(length(SNR)/gap):length(SNR);
    plot(SNR(downsampled), Exp_R_id(downsampled), 'cd', 'LineWidth', 1);
    lp1 = plot(0,0,'-c','Marker','d', 'LineWidth', 1, 'visible','off');
end
hold on,
%h2 = plot(SNR, reshape(Exp_opt_R_th(1,:,1), 1, length(SNR)), 'k-.', 'LineWidth', 1);
%hold on,
h33 = plot(SNR, reshape(Exp_opt_R_th(2,:,1), 1, length(SNR)), 'k-', 'LineWidth', 1);
hold on,
plot(SNR, reshape(Exp_opt_R_th(2,:,2), 1, length(SNR)), 'k-', 'LineWidth', 1);


if ver_paper
    hold on,
    plot(SNR(downsampled), reshape(Exp_opt_R_th(2,downsampled,1), 1, length(downsampled)), 'ko', 'LineWidth', 1);
    %hold on,
    %plot(SNR, reshape(Exp_opt_R_th(1,:,2), 1, length(SNR)), 'k-.', 'LineWidth', 1);
    hold on,
    plot(SNR(downsampled), reshape(Exp_opt_R_th(2,downsampled,2), 1, length(downsampled)), 'ko', 'LineWidth', 1);
    lp2 = plot(0,0,'-k','Marker','o', 'LineWidth', 1, 'visible','off');
end

grid on;
axis([min(SNR) max(SNR) 0  max(Exp_R_id) * 1.01]);
ylabel('$\rs(\ttau)$ [bits/sec/Hz]','FontSize',Fontsize);
xlabel('$\gamma$ [dB]','FontSize',Fontsize);

alt = 00;
if alt
    %% Alternative Strategy --- With no power control
    N_alt = 10:10:100000;   
    alpha_true_p1 = [120:-0.05:90];                    % True Path loss between ST and PR  to achieve P_p = IT  
    snr = 10^(P_tran/10) * 10.^(-alpha_true_p1/10)...  % snr received at the PR
        / 10^(noise_power/10);
    Exp_opt_R_alt = zeros(length(epsilon), length(snr), length(P_reg_max));
    N_alt_opt = zeros(length(epsilon), length(snr), length(P_reg_max));

    for i = 1:length(epsilon)
        for j = 1:length(P_reg_max)
            for k = 1:length(snr)


                %% Gamma Approximation to the non-central chi-squared distribution
                mean = 10^(noise_power/10) * (1 + snr(k));       
                var = (10^(noise_power/10))^2./N_alt * ( 2 + 4 * snr(k)); 

                b = var/mean;
                a = mean./b;

                %% Determine the SNR that achieves the outage probability
                temp1 = gammainc((10^(I_T/10) * 10^(P_tran/10)...
                   /P_reg_max(j) + 10^(noise_power/10)) * 1./b, a, 'upper');
               [value index] = min(abs(temp1 - epsilon(i)));  
               index
               N_alt_opt(i,k,j) = index; 
               if (index~=length(N_alt)) && value < epsilon(i)
                   C_1  = calc_capacities(10^(-alpha_true_p2/10)/10^(noise_power/10),...
                            10^(-alpha_true_s/10), P_reg_max(j), 10^(noise_power/10), N_alt(index), N_s)
                   Exp_opt_R_alt(i,k,j) = (K - N_alt(index))/K * C_1;
               else
                   break;
               end
            end
        end
    end
    save('results_opt_thr_vs_SNR_AWGN_SI_10_alt.mat');
end
load('results_opt_thr_vs_SNR_AWGN_SI_10_alt.mat');
SNR_alt = 10*log10(snr);

%hold on,
%plot(SNR, reshape(Exp_opt_R_alt(1,:,1), 1, length(SNR)), 'r--', 'LineWidth', 1);
hold on,
h44  = plot(SNR_alt, reshape(Exp_opt_R_alt(2,:,1), 1, length(SNR_alt)), 'r--', 'LineWidth', 1);


%hold on,
%plot(SNR, reshape(Exp_opt_R_alt(1,:,2), 1, length(SNR)), 'r-', 'LineWidth', 1);
hold on,
plot(SNR_alt, reshape(Exp_opt_R_alt(2,:,2), 1, length(SNR_alt)), 'r--', 'LineWidth', 1);

if ver_paper
    downsampled_alt = [];
    for i = downsampled
        downsampled_alt(end + 1) = find(SNR_alt == SNR(i));
    end
    hold on,
    plot(SNR_alt(downsampled_alt), reshape(Exp_opt_R_alt(2,downsampled_alt,1), 1, length(downsampled_alt)), 'rs', 'LineWidth', 1);
    hold on,
    plot(SNR_alt(downsampled_alt), reshape(Exp_opt_R_alt(2,downsampled_alt,2), 1, length(downsampled_alt)), 'rs', 'LineWidth', 1);
    lp3 = plot(0,0,'--r','Marker','s', 'LineWidth', 1, 'visible','off');
end

%hl = legend([h1 h3 h2],'IM', 'EM, $\opc = 0.1$', 'EM, $\opc = 0.01$');
if ver_paper
    hl = legend([lp1 lp2 lp3], 'IM', 'EM','Coro 2');
else
    hl = legend([h11 h33 h44], 'IM', 'EM','Coro 2');
end
set(hl, 'position',[0.66 0.74 0.24 0.18]);
set(gca,'FontSize',Fontsize);

if 00
    %%%% Zommed Image %%%%%
    % create a new pair of axes inside current figure
    ax = axes('Position',[.56 .51 .33 .33]);
    box on; % put box around new pair of axes
    child_gca = get(gca,'Children');

    [value indexOfInterest] = find(SNR <= -0.0 & SNR >= -1.0);
    plot(SNR(indexOfInterest), reshape(Exp_opt_R_th(1,indexOfInterest,1), 1, length(SNR(indexOfInterest))),...
    'k-.', 'LineWidth', 1); % plot on new axes
    hold on,
    plot(SNR(indexOfInterest), reshape(Exp_opt_R_th(2,indexOfInterest,1), 1, length(SNR(indexOfInterest))),...
    'k-', 'LineWidth', 1); 
    hold on,
    plot(SNR(indexOfInterest), reshape(Exp_opt_R_th(1,indexOfInterest,2), 1, length(SNR(indexOfInterest))),...
    'k-.', 'LineWidth', 1); 
    hold on,
    plot(SNR(indexOfInterest), reshape(Exp_opt_R_th(2,indexOfInterest,2), 1, length(SNR(indexOfInterest))),...
    'k-', 'LineWidth', 1); 
    plot(SNR(indexOfInterest), Exp_R_id(indexOfInterest), 'c', 'LineWidth', 3);


    % plot on new axes
    %plot(SNR(indexOfInterest), Exp_opt_R_th(3,indexOfInterest), 'k', 'LineWidth', 1); 
    %ylabel(gca, '$\rs(\ttau)$','FontSize',Fontsize);
    %xlabel(gca, '$\gamma$','FontSize',Fontsize);
    %set(ax, 'XTick', [-1.3 -0.9 -0.5]);
    %set(ax, 'YTick', [3.2 3.4 3.7]);
    axis(ax, [min(SNR(indexOfInterest)) max(SNR(indexOfInterest))...
        min(reshape(Exp_opt_R_th(1,indexOfInterest,1), 1, length(SNR(indexOfInterest))))...
        max(Exp_R_id(indexOfInterest)) * 1]);
    set(gca,'FontSize',Fontsize);
    th = title(ax, 'Zoom');
    set(th, 'Position', [-0.85    3.4    1.0001]);
end
if ver_paper
    laprint(2, '../figures/fig_opt_thr_vs_SNR_AWGN_SI_10', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
else
    laprint(2, '../figures/fig_opt_thr_vs_SNR_AWGN_SI_10_diss', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
end
if 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Optium estimatio time vs SNR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2)
    load('results_opt_thr_vs_SNR_AWGN_e01_th.mat');
    SNR = 10*log10(snr);
    h2 = plot(SNR, Opt_N_th(1,:) * 1e-3, 'k', 'LineWidth', 1);
    hold on,
    h3 = plot(SNR, Opt_N_th(2,:) * 1e-3 , 'k', 'LineWidth', 1);
    %hold on,
    %h4 = plot(SNR, Opt_N_th(3,:) * 1e-3, 'k', 'LineWidth', 1);
    grid on;
    axis([min(SNR) max(SNR) min(Opt_N_th(3,:)) * 1e-3  max(Opt_N_th(1,:)) * 1e-3]);
    ylabel('$\ttau$ [ms]','FontSize',Fontsize);
    xlabel('$\gamma$ [dB]','FontSize',Fontsize);
    %hl = legend([h2 h3 h4],'IM','EM');
    %set(hl, 'position',[0.725 0.12 0.15 0.31]);
    set(hl, 'position',[0.68 0.71 0.22 0.2]);
    set(gca,'FontSize',Fontsize);
    laprint(2, '../figures/fig_opt_est_time_vs_SNR_AWGN', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'off');
end