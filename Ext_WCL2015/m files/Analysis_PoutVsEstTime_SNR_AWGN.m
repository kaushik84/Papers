%                                
%        Description: This m-file we analyze two issues:
%                     1) We consider that positive values of control power, 
%                     hence depict the variations of the outage
%                     probability w.r.t to the estimation time and SNR
%                     2) We consider that application of maximum transmit power, 
%                     hence depict the variations of the outage
%                     probability w.r.t to the estimation time and SNR   
%        Clearly the regulated power is distributed according non central
%        chi-2 distribution. However to simplify the analysis we
%        approximate the non central distribution against the Gamma
%        distribution.
%        
%        The simulation is performed to show the following:
%        1) Analyse the optimum Controlled Power vs Outage Probability
%
%        Created on: 04.08.15
%        Revision History: 04.08.15 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%% Analysis 1 -- Variation of P_out vs SNR for different P_reg and certain N_th
if 0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P_tran = -00;                                     % Power transmitted by PR, the SNR received at ST can be 
                                                      % varried using this parameter
    noise_power = -100;                               % noise power -100 dBm
    I_T = -110;                                       % Interference temperature -80 dBm
    f_s = 1e6;                                        % 1 MHz one band
    K = 0.1 * f_s;                                    % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
    alpha_true_p = [120:-0.005:90];                     % True Path loss between ST and PR  to achieve P_p = IT  
    alpha_true_s = 080;                               % True Path loss between ST and SR  to achieve P_p = IT  

    snr = 10^(P_tran/10) * 10.^(-alpha_true_p/10)...  % snr received at the PR
        / 10^(noise_power/10);

    P_reg_max = [-10 00 10];                           % Maximum Controlled Power at the ST

    N_th = [1000 5000 10000];                          % N = Total of samples used for estimation = tau * f_s 
    epsilon_I1 = zeros(length(N_th), length(snr), 1);     % Outage Probability constraint at PR
    epsilon_I2 = zeros(length(N_th), length(snr), length(P_reg_max));     % Outage Probability constraint at PR operating at max Tx power
    epsilon_I3 = zeros(length(N_th), length(snr), length(P_reg_max));     % Outage Probability constraint at PR operating at max Tx power
    

    %% Issue 1
     if 1
         for i=1:length(N_th)
             for j=1:length(snr)
                       %disp(strcat('N = ',num2str(N_th(k))));             
    
                       %% Determining the performance meterics --> Exp_R
                       %% Expected values           

                       %% Gamma Approximation to the non-central chi-squared distribution
                       mean = 10^(noise_power/10) * (1 + snr(j));       
                       var = (10^(noise_power/10))^2/N_th(i) * ( 2 + 4 * snr(j)); 

                       b = var/mean;
                       a = mean/b;

                       %% Determine the outage probability
                       epsilon_I1(i, j, 1) = 1 - gammainc((10^(noise_power/10)/b), a);

                 end 
         end
     end


    %% Issue 2
     if 1
         for i=1:length(N_th)
             disp(strcat('N = ',num2str(N_th(i))));             
             for j=1:length(snr)
                 for k=1:length(P_reg_max)       

                       %% Determining the performance meterics --> Exp_R
                       %% Expected values           

                       %% Gamma Approximation to the non-central chi-squared distribution
                       mean = 10^(noise_power/10) * (1 + snr(j));       
                       var = (10^(noise_power/10))^2/N_th(i) * ( 2 + 4 * snr(j)); 

                       b = var/mean;
                       a = mean/b;


                       %% Determine the outage probability Operating Maximum Transmit Power 
                       epsilon_I2(i, j, k) = 1 - gammainc((10^(I_T/10) * 10^(P_tran/10)...
                           /10^(P_reg_max(k)/10) + 10^(noise_power/10))/b, a);

                       if 0  %% No Success
                           %% Determine the outage probability with regulated Transmit Power and Maxmim Transmit Power Constraint

                           P_reg = min(10^(P_reg_max(k)/10), 10^(I_T/10) * 10^(P_tran/10) ./...
                                (b .* gammaincinv((1 - epsilon_I2(j,k)),a) - 10^(noise_power/10))); 

                            epsilon_I3(i, j, k) = 1 - gammainc((10^(I_T/10) * 10^(P_tran/10)...
                                /P_reg + 10^(noise_power/10))/b, a);
                       end
                 end 
             end
         end
     end
    SNR = 10*log10(snr);
    Fontsize = 9;
    figure(1);
    h1 = plot(SNR, reshape(epsilon_I1(1,:,1), 1, length(SNR)), 'k', 'LineWidth', 1);
    %hold on,
    %plot(SNR, reshape(epsilon_I1(2,:,1), 1, length(SNR)), 'k', 'LineWidth', 1);
    hold on, 
    h2 = plot(SNR, reshape(epsilon_I1(3,:,1), 1, length(SNR)), 'k--', 'LineWidth', 1);
    
    %hold on, 
    %plot(SNR, reshape(epsilon_I2(1,:,1), 1, length(SNR)), 'k--', 'LineWidth', 1);
    %hold on,
    %plot(SNR, reshape(epsilon_I2(2,:,1), 1, length(SNR)), 'k--', 'LineWidth', 1);
    %hold on, 
    %plot(SNR, reshape(epsilon_I2(3,:,1), 1, length(SNR)), 'k--', 'LineWidth', 1);   
    hold on, 
    plot(SNR, reshape(epsilon_I2(1,:,2), 1, length(SNR)), 'k-', 'LineWidth', 1);
    %hold on,
    %plot(SNR, reshape(epsilon_I2(2,:,2), 1, length(SNR)), 'k--', 'LineWidth', 1);
    hold on, 
    plot(SNR, reshape(epsilon_I2(3,:,2), 1, length(SNR)), 'k--', 'LineWidth', 1);    
    hold on, 
    plot(SNR, reshape(epsilon_I2(1,:,3), 1, length(SNR)), 'k-', 'LineWidth', 1);
    %hold on,
    %plot(SNR, reshape(epsilon_I2(2,:,3), 1, length(SNR)), 'k--', 'LineWidth', 1);
    hold on, 
    plot(SNR, reshape(epsilon_I2(3,:,3), 1, length(SNR)), 'k--', 'LineWidth', 1);    
    axis([-20 -6 0 1]);
    grid on;
    ylabel('$\opc$','FontSize',Fontsize);
    xlabel('$\gamma$ [dB]','FontSize',Fontsize);
    hl = legend([h1 h2], '$\tau = 1$ ms', '$\tau = 10$ ms');
    %set(hl, 'position',[0.725 0.12 0.15 0.31]);
    %set(hl, 'position',[0.15 0.71 0.27 0.2]);
    set(gca,'FontSize',Fontsize);
    laprint(1, '../figures/fig_opc_vs_SNR_diff_N_diff_maxContPow_th', 'options', 'factory', 'width', 8, 'scalefonts',...
       'on', 'factor',0.5, 'keepfontprops', 'on'); 
    
end


%% Analysis 2 -- Variation of SNR vs N_th for different P_reg and certain P_out
if 1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P_tran = -00;                                     % Power transmitted by PR, the SNR received at ST can be 
                                                      % varried using this parameter
    noise_power = -100;                               % noise power -100 dBm
    I_T = -110;                                       % Interference temperature -80 dBm
    f_s = 1e6;                                        % 1 MHz one band
    K = 0.1 * f_s;                                    % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
    alpha_true_p = [120:-0.005:90];                    % True Path loss between ST and PR  to achieve P_p = IT  
    alpha_true_s = 080;                               % True Path loss between ST and SR  to achieve P_p = IT  

    snr = 10^(P_tran/10) * 10.^(-alpha_true_p/10)...  % snr received at the PR
        / 10^(noise_power/10);

    P_reg_max = [00];                              % Maximum Controlled Power at the ST

    epsilon = [0.1];                              

    
    N_th = [10:50:20000];                             % N = Total of samples used for estimation = tau * f_s 
    SNR_I1 = zeros(1, length(epsilon),...
        length(N_th));                                % Outage Probability constraint at PR
    SNR_I2 = zeros(length(P_reg_max), length(epsilon),....
        length(N_th));                                % Outage Probability constraint at PR
    
    SNR_I3 = zeros(length(P_reg_max), length(epsilon),....
        length(N_th));                                % With Full Power and transmit power constraint  
    
    SNR_PA = zeros(length(P_reg_max), length(epsilon),....
        length(N_th));                                % With Controlled Power and transmit power constraint  
    
    
    th = 00;
    
    if th
        %% Issue 1
         if 0
                 for j=1:length(epsilon)
                     disp(strcat('epsilon = ',num2str(epsilon(j))));   
                     for k=1:length(N_th)       

                           %% Determining the performance meterics --> Exp_R
                           %% Expected values           

                           %% Gamma Approximation to the non-central chi-squared distribution
                           mean = 10^(noise_power/10) * (1 + snr);       
                           var = (10^(noise_power/10))^2/N_th(k) * ( 2 + 4 * snr); 

                           b = var./mean;
                           a = mean./b;

                         %% Determine the SNR that achieves the outage probability at the maximum Transmit Power
                           temp = 1 - gammainc((10^(noise_power/10)./b), a);
                           [value index] = min(abs(temp - epsilon(j)));
                           SNR_I1(1, j,k) = 10 * log10(snr(index));                   

                     end 
                 end
         end


        %% Issue 2
         if 1    
             for i=1:length(P_reg_max)
                 disp(strcat('P_reg_max = ',num2str(P_reg_max(i))));   
                 for j=1:length(epsilon)
                    disp(strcat('epsilon = ',num2str(epsilon(j))));   
                    for k=1:length(N_th)       

                    %% Gamma Approximation to the non-central chi-squared distribution
                       mean = 10^(noise_power/10) * (1 + snr);       
                       var = (10^(noise_power/10))^2/N_th(k) * ( 2 + 4 * snr); 

                       b = var./mean;
                       a = mean./b;

                    %% Determine the SNR that achieves the outage probability
                       temp1 = gammainc((10^(I_T/10) * 10^(P_tran/10)...
                           /10^(P_reg_max(i)/10) + 10^(noise_power/10)) * 1./b, a, 'upper');
                       [value index1] = min(abs(temp1 - epsilon(j)));
                       10 * log10(snr(index1))
                       SNR_I2(i,j,k) = 10 * log10(snr(index1));


                    %% Determine the SNR that achieves the outage probability at the regulated Power                       
                       if 0
                           P_reg = min(10^(P_reg_max(i)/10), 10^(I_T/10) * 10^(P_tran/10) ./...
                                (b .* gammaincinv((1 - epsilon(j)),a) - 10^(noise_power/10)));  
                           temp2 = gammainc((10^(I_T/10) * 10^(P_tran/10)...
                               ./P_reg + 10^(noise_power/10))./b, a, 'upper'); 
                           [value index2] = find(temp2 >= epsilon(j));
                           10 * log10(snr(min(index2)))
                           SNR_PA(i,j,k) = 10 * log10(snr(min(index2)));
                       end
                       
                    end
                 end
             end
             save('results_N_vs_SNR_diff_pout_diff_maxContPow_AWGN_th');
         end
     end
     load('results_N_vs_SNR_diff_pout_diff_maxContPow_AWGN_th');
     Fontsize = 8.5;
     figure(1);
     %h1 = plot(reshape(SNR_I1(1,1,:),1, length(N_th)), N_th * 1e-3,  'c', 'LineWidth', 2); 
     %hold on,
     %h2 = plot(reshape(SNR_I2(1,1,:),1, length(N_th)), N_th * 1e-3,  'k-', 'LineWidth', 1); 
     %hold on, 
     h3 = plot(reshape(SNR_I2(1,1,:),1, length(N_th)), N_th * 1e-3,  'k-', 'LineWidth', 1); 
     %hold on, 
     %plot(reshape(SNR_PA(1,1,:),1, length(N_th)), N_th * 1e-3,  'k-.', 'LineWidth', 1); 
     %h3 = plot(reshape(SNR_I3(1,1,:),1, length(N_th)), N_th * 1e-3,  'k-', 'LineWidth', 1); 
     %hold on, 
     %h4 = plot(reshape(SNR_I3(1,2,:),1, length(N_th)), N_th * 1e-3,  'k--', 'LineWidth', 1);     
     
     %hold on,
     %plot(reshape(SNR_I1(1,2,:),1, length(N_th)), N_th * 1e-3,  'c', 'LineWidth', 2); 
     %hold on,
     %plot(reshape(SNR_I2(2,1,:),1, length(N_th)), N_th * 1e-3,  'k-', 'LineWidth', 1); 
     %hold on, 
     %plot(reshape(SNR_I2(2,2,:),1, length(N_th)), N_th * 1e-3,  'k--', 'LineWidth', 1); 
     %hold on, 
     %plot(reshape(SNR_I2(2,3,:),1, length(N_th)), N_th * 1e-3,  'k-.', 'LineWidth', 1); 
     
     %hold on,
     %plot(reshape(SNR_I2(3,1,:),1, length(N_th)), N_th * 1e-3,  'k-', 'LineWidth', 1); 
     %hold on, 
     %plot(reshape(SNR_I2(3,2,:),1, length(N_th)), N_th * 1e-3,  'k--', 'LineWidth', 1); 
     %hold on, 
     %plot(reshape(SNR_I2(3,3,:),1, length(N_th)), N_th * 1e-3,  'k-.', 'LineWidth', 1); 
     axis([-20 -5 0 20]);
     %grid on;
     ylabel('$\tau$ [ms]','FontSize',Fontsize);
     xlabel('$\gamma$ [dB]','FontSize',Fontsize);
     %hl = legend([h2 h3], '$\opc = 0.01$', '$\opc = 0.1$');
     %set(hl, 'position',[0.725 0.12 0.15 0.31]);
     %set(hl, 'position',[0.14 0.71 0.265 0.2]);
     set(gca,'FontSize',Fontsize);
     laprint(1, '../figures/fig_N_vs_SNR_diff_pout_diff_maxContPow_th', 'options', 'factory', 'width', 8, 'scalefonts',...
         'on', 'factor',0.5, 'keepfontprops', 'on'); 
end


%% Analysis 3 -- Variation of P_reg along the SNR for a fixed N_th
if 0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P_tran = -00;                                     % Power transmitted by PR, the SNR received at ST can be 
                                                      % varried using this parameter
    noise_power = -100;                               % noise power -100 dBm
    I_T = -110;                                       % Interference temperature -80 dBm
    f_s = 1e6;                                        % 1 MHz one band
    K = 0.1 * f_s;                                    % K = Total number of samples in a frame = T * f_s, T = 100 ms WRAN standard
    alpha_true_p = [120:-0.01:90];                    % True Path loss between ST and PR  to achieve P_p = IT  
    alpha_true_s = 080;                               % True Path loss between ST and SR  to achieve P_p = IT  

    snr = 10^(P_tran/10) * 10.^(-alpha_true_p/10)...  % snr received at the PR
        / 10^(noise_power/10);

    P_reg_max = [-10 00 20];                         % Maximum Controlled Power at the ST

    epsilon = [0.01 0.1 0.2];                              

    P_reg = zeros(1, length(snr));    
    N_th = 10000;                              % N = Total of samples used for estimation = tau * f_s 
    P_reg_max = 10^(-00/10);                       % Maximum Trasmit Power Constraint
  
    
     if 1
             for j=1:length(epsilon)
                 disp(strcat('epsilon = ',num2str(epsilon(j))));   
                 for k=1:length(snr)       

                       %% Determining the performance meterics --> Exp_R
                       %% Expected values           

                       %% Gamma Approximation to the non-central chi-squared distribution
                       mean = 10^(noise_power/10) * (1 + snr(k));       
                       var = (10^(noise_power/10))^2/N_th * ( 2 + 4 * snr(k)); 

                       b = var/mean;
                       a = mean/b;

                       %% Determine the SNR that achieves the outage probability
                       P_reg(j,k) = min(P_reg_max, 10^(I_T/10) * 10^(P_tran/10) /...
                            (b * gammaincinv((1 - epsilon(j)),a) - 10^(noise_power/10)));
            
                 end 
             end
             SNR = 10*log10(snr);
             figure(1);
             plot(SNR, 10 * log10(P_reg(1,:)), 'k-', 'LineWidth', 1);
             hold on,
             plot(SNR, 10 * log10(P_reg(2,:)), 'k--', 'LineWidth', 1);
             hold on,
             plot(SNR, 10 * log10(P_reg(3,:)), 'k-.', 'LineWidth', 1); 
             axis([-20  10 -25 1]);
     end     
end
