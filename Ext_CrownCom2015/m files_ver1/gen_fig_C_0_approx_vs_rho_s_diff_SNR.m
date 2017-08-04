%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: This m-file validates through simulation and Gamma 
%        approximation made determining the analytical expressions for the 
%        expected capapacites C_0.
%
%        C_0:       log2(1 + snr_s)   
%
%        According to the channel estimation model the amplitude gain of 
%        the channel estimated at the SR for the PT-PR and ST-SR is 
%        Gausssian distributed 
%            h_s ~ N(h_s * \pho_s, sqrt(1 - rho_s^2)),
%            h_p ~ N(h_p * \pho_s, sqrt(1 - rho_p^2))  
%        
%        Clearly, the power gains (g_s = |h_s|^2 and g_p = |h_p|^2) are
%        distributed according to non central chi-2 distribution with 1
%        degree, which can be approximated as Gamma distribution to ease
%        the analytical analysis. Based on we perform simulation to
%        validate the analytical (Gamma) approximation 
%
%
%        Created on: 31.03.15
%        Last modified: 31.03.15
%        Revision History: 31.03.15 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;



sim = 00;                                       % Enable( = 1) to perform simulation, theoretical analysis
th = 00;                                        % disable ( = 0) to plot curves, data is read from a file
                                                % Theoretical anaylsis is also included of simulation as
                                                % numerical integration is involved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
snr_s = 10.^([-10,0,10]/10);                    % snr received at the SR, over the link ST-SR 
dof_s = 1;                                      % degree of freedom for the estimated channel, usually set to 1  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sim   
    %% Simulation parameters
    M = 1e6;                                                % Number of Realization
    g_s_sim = zeros(1,M);                                   % power gain for the channel g_s 
    rho_s_sim = 0.1:0.1:0.95;                               % Channel correlation coefficient                                            
    Exp_C_0_sim = zeros(length(snr_s),length(rho_s_sim));   % Expected Capacity

    for i = 1:length(snr_s)
        disp(strcat('SIM:snr_s(i) = ',num2str(snr_s(i))));
        for j = 1:length(rho_s_sim)
            disp(strcat('SIM:rho_s_sim(j) = ',num2str(rho_s_sim(j))));
            parfor k = 1:M
                g_s_sim(k) =  sum( normrnd(rho_s_sim(j), sqrt(snr_s(i))...
                    * sqrt(1 - rho_s_sim(j)^2), 1 ,dof_s).^2);
            end
            Exp_C_0_sim(i,j) = mean(log2(1 + g_s_sim));
        end
    end
    save('results_C_0_approx_vs_rho_s_diff_SNR_sim.mat');
    if ~th      % If theoretical analysis is already simulated the computation ends here.
        quit;
    end
end    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th
   %% Theoretical parameters
    rho_s_th = 0.02:0.02:0.98;                              % Channel correlation coefficient                                            
    Exp_C_0_th = zeros(length(snr_s),length(rho_s_th));   % Expected Capacity

    for i = 1:length(snr_s)           
        disp(strcat('y:snr_s(i) = ',num2str(snr_s(i))));        
        for j = 1:length(rho_s_th)
            disp(strcat('TH:rho_s_th(j) = ',num2str(rho_s_th(j))));

            %% Gamma Parameters
            lambda_s = dof_s *  rho_s_th(j)^2/(snr_s(i) * (1 - rho_s_th(j)^2));   
            mean_ana_s = snr_s(i) * (1 - rho_s_th(j)^2) * (dof_s + lambda_s);    
            var_ana_s = snr_s(i)^2 * (1 - rho_s_th(j)^2)^2 * (2 * dof_s + 4 * lambda_s);
    
            b_gamma_s = var_ana_s/mean_ana_s;
            a_gamma_s = mean_ana_s/b_gamma_s;
            
            func_exp_C_0 = @(t) t .* 2.^(t) * log(2) .* 1/b_gamma_s^(a_gamma_s)...
                .* (2.^t - 1).^(a_gamma_s - 1)/gamma(a_gamma_s)...
                .* exp(- (2.^t - 1)/b_gamma_s);
            Exp_C_0_th(i,j) = integral(func_exp_C_0,0, 10);            
        end
    end
    save('results_C_0_approx_vs_rho_s_diff_SNR_th.mat');
    quit;
end

load('results_C_0_approx_vs_rho_s_diff_SNR_th.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis --  Exp_C_0 vs vs Channel
%   estimatoin correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fontsize = 9;
plot(rho_s_th, Exp_C_0_th(1,:),'-', 'LineWidth', 1.5);
hold on, 
plot(rho_s_th, Exp_C_0_th(2,:),'--', 'LineWidth', 1.5, 'HandleVisibility','off');
hold on,
plot(rho_s_th, Exp_C_0_th(3,:),'-.', 'LineWidth', 1.5, 'HandleVisibility','off');

load('results_C_0_approx_vs_rho_s_diff_SNR_sim.mat');

hold on,
plot(rho_s_sim, Exp_C_0_sim(1,:), 'r*');
hold on, 
plot(rho_s_sim, Exp_C_0_sim(2,:), 'r*');
hold on,
plot(rho_s_sim, Exp_C_0_sim(3,:), 'r*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves cosmetic makeover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid on;
%axis tight;
axis([min(rho_s_th) max(rho_s_th)...
    min(Exp_C_0_th(1,:)) * 0.0 max(Exp_C_0_th(3,:)) * 1.05]);
ylabel('$\C0$','FontSize',Fontsize+2);
xlabel('$\rho_s$ = [ms]','FontSize',Fontsize+2);
hl = legend('Theoretical', 'Simulated');
set(hl, 'Location', 'NorthEast', 'FontSize', Fontsize);
set(gca,'FontSize',Fontsize);
laprint(1, '../figures/fig_C_0_approx_vs_rho_s_diff_SNR', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');

