%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: This m-file validates through simulation and Gamma 
%        approximation made determining the analytical expressions for the 
%        expected capapacites C_0.
%
%        C_1:       log2(1 + (snr_s/(1 + snr_p))   
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
snr_s = 10.^([-10 00 10]/10);                   % snr received at the SR, over the link ST-SR 
snr_p = 10.^([-10 00 10]/10);                   % snr received at the SR, over the link PT-SR 
dof_s = 1;                                      % degree of freedom for the estimated channel, usually set to 1  
dof_p = 1;                                      % degree of freedom for the estimated channel, usually set to 1  

rho_s = 0.5;                                    % Channel correlation coefficient ST-SR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sim   
    %% Simulation parameters
    M = 1e5;                                                % Number of Realization
    rho_p_sim = 0.1:0.1:0.95;                               % Channel correlation coefficient PT-SR    

    Exp_C_1_sim = zeros(length(snr_p),length(snr_s),...     % Expected Capacity
        length(rho_p_sim));   

    for i = 1:length(snr_p)
        disp(strcat('SIM:snr_p(i) = ',num2str(snr_p(i))));
        for j = 1:length(snr_s)
            disp(strcat('SIM:snr_s(j) = ',num2str(snr_s(j))));
            for k = 1:length(rho_p_sim)
                disp(strcat('SIM:rho_p_sim(j) = ',num2str(rho_p_sim(k))));
                g_s_sim = zeros(1,M);                                   % power gain for the channel g_s 
                g_p_sim = zeros(1,M);                                   % power gain for the channel g_p    
                parfor l = 1:M
                    g_s_sim(l) =  sum( normrnd(rho_s, sqrt(snr_s(j))...
                        * sqrt(1 - rho_s^2), 1 ,dof_s).^2);
                    g_p_sim(l) =  sum( normrnd(rho_p_sim(k), sqrt(snr_p(i))...
                        * sqrt(1 - rho_p_sim(k)^2), 1 ,dof_p).^2);
                end
                Exp_C_1_sim(i,j,k) = mean(log2(ones(1,M) + g_s_sim./(ones(1,M) + g_p_sim)));  
            end
        end
    end
    save('results_C_1_approx_vs_rho_p_diff_SNR_sim.mat', 'snr_p','snr_s',...
        'Exp_C_1_sim', 'rho_p_sim', 'rho_s');
    if ~th      % If theoretical analysis is already simulated the computation ends here.
        quit;
    end
end    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if th
   %% Theoretical parameters
    rho_p_th = 0.02:0.02:0.98;                              % Channel correlation coefficient                                            
    Exp_C_1_th = zeros(length(snr_s),length(snr_s),...      % Expected Capacity
    length(rho_p_th));    

    for i = 1:length(snr_p)           
        disp(strcat('TH:snr_p(i) = ',num2str(snr_p(i))));           
        for j = 1:length(snr_s)           
            disp(strcat('TH:snr_s(j) = ',num2str(snr_s(j))));  
            for k = 1:length(rho_p_th)
                disp(strcat('TH:rho_p_th(k) = ',num2str(rho_p_th(k))));

                %% Gamma Parameters
                lambda_s = dof_s *  rho_s^2/(snr_s(j) * (1 - rho_s^2));   
                mean_ana_s = snr_s(j) * (1 - rho_s^2) * (dof_s + lambda_s);    
                var_ana_s = snr_s(j)^2 * (1 - rho_s^2)^2 * (2 * dof_s + 4 * lambda_s);

                b_gamma_s = var_ana_s/mean_ana_s;
                a_gamma_s = mean_ana_s/b_gamma_s;

                lambda_p = dof_p *  rho_p_th(k)^2/(snr_p(i) * (1 - rho_p_th(k)^2));   
                mean_ana_p = snr_p(i) * (1 - rho_p_th(k)^2) * (dof_p + lambda_p);    
                var_ana_p = snr_p(i)^2 * (1 - rho_p_th(k)^2)^2 * (2 * dof_p + 4 * lambda_p);

                b_gamma_p2 = var_ana_p/mean_ana_p;
                a_gamma_p2 = mean_ana_p/b_gamma_p2;

                %% Single Integral
                if 0
                    f_exp_C_1_s = @(z) z .* (2.^z) * log(2) * 1/(b_gamma_p2^a_gamma_p2 * b_gamma_s^a_gamma_s) *...
                        exp(1/b_gamma_p2) .* (2.^z - 1).^(a_gamma_s -1) .* ((1/b_gamma_p2 +...
                        (2.^z - 1)/b_gamma_s).^(-a_gamma_s - a_gamma_p2) * gamma(a_gamma_p2 + a_gamma_s)...
                        .* hypergeom(1 - a_gamma_p2, 1 - a_gamma_p2 - a_gamma_s,...
                        -(b_gamma_s + b_gamma_p2 * (2.^z - 1))/(b_gamma_s * b_gamma_p2))...
                        + gamma(-a_gamma_p2 - a_gamma_s) * 1/gamma(-a_gamma_s) * gamma(a_gamma_p2) *...
                        hypergeom(1 + a_gamma_s, 1 + a_gamma_p2 + a_gamma_s,...
                        -(b_gamma_s + b_gamma_p2 * (2.^z - 1))/(b_gamma_s * b_gamma_p2)))/...
                        (gamma(a_gamma_s) * gamma(a_gamma_p2));


                    Exp_C_1_th(i,j) = abs(integral(f_exp_C_1_s,0, 25));  
                end
                %% Double Integral
                if 1
                    f_exp_C_1_s = @(y,z) z .* (2.^z) .* log(2) * 1/(gamma(a_gamma_s) * b_gamma_s^a_gamma_s)...
                       * 1/(gamma(a_gamma_p2) * b_gamma_p2^a_gamma_p2) .* y .* ((2.^z - 1) .* y).^(a_gamma_s-1)...
                       .* exp(- (2.^z- 1)  .* y/ b_gamma_s) .* (y - 1).^(a_gamma_p2 - 1) .*...
                       exp(-(y-1)/b_gamma_p2);
                    if rho_p_th(k) == 0.98
                        Exp_C_1_th(i,j,k) = integral2(f_exp_C_1_s, 1, 100, 0, 25);
                    else
                       Exp_C_1_th(i,j,k) = integral2(f_exp_C_1_s, 1, Inf, 0, 25);
                    end
                end
            end            
        end
    end
    save('results_C_1_approx_vs_rho_p_diff_SNR_th.mat');
    %quit;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis --  Exp_C_0 vs vs Channel
%   estimatoin correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fontsize = 9;
if 1
    load('results_C_1_approx_vs_rho_p_diff_SNR_th.mat');
    plot(rho_p_th, reshape(Exp_C_1_th(1,1,:), 1, length(rho_p_th)),'-', 'LineWidth', 1.5);
    %hold on, 
    %plot(rho_p_th, reshape(Exp_C_1_th(1,2,:), 1, length(rho_p_th)),'-', 'LineWidth', 1.5, 'HandleVisibility','off');
    hold on, 
    plot(rho_p_th, reshape(Exp_C_1_th(1,3,:), 1, length(rho_p_th)),'-', 'LineWidth', 1.5, 'HandleVisibility','off');
    hold on, 
    plot(rho_p_th, reshape(Exp_C_1_th(2,1,:), 1, length(rho_p_th)),'--', 'LineWidth', 1.5, 'HandleVisibility','off');
    %hold on, 
    %plot(rho_p_th, reshape(Exp_C_1_th(2,2,:), 1, length(rho_p_th)),'--', 'LineWidth', 1.5, 'HandleVisibility','off');
    hold on, 
    plot(rho_p_th, reshape(Exp_C_1_th(2,3,:), 1, length(rho_p_th)),'--', 'LineWidth', 1.5, 'HandleVisibility','off');
    hold on,
    plot(rho_p_th, reshape(Exp_C_1_th(3,1,:), 1, length(rho_p_th)),'-.', 'LineWidth', 1.5, 'HandleVisibility','off');
    %hold on,
    %plot(rho_p_th, reshape(Exp_C_1_th(3,2,:), 1, length(rho_p_th)),'-.', 'LineWidth', 1.5, 'HandleVisibility','off');
    hold on,
    plot(rho_p_th, reshape(Exp_C_1_th(3,3,:), 1, length(rho_p_th)),'-.', 'LineWidth', 1.5, 'HandleVisibility','off');    
end
if 1
    load('results_C_1_approx_vs_rho_p_diff_SNR_sim.mat');
    hold on,
    plot(rho_p_sim, reshape(Exp_C_1_sim(1,1,:), 1, length(rho_p_sim)), 'r*');    
    %hold on,
    %plot(rho_p_sim, reshape(Exp_C_1_sim(1,2,:), 1, length(rho_p_sim)), 'r*', 'HandleVisibility','off');    
    hold on,
    plot(rho_p_sim, reshape(Exp_C_1_sim(1,3,:), 1, length(rho_p_sim)), 'r*', 'HandleVisibility','off');
    hold on, 
    plot(rho_p_sim, reshape(Exp_C_1_sim(2,1,:), 1, length(rho_p_sim)), 'r*', 'HandleVisibility','off');    
    %hold on, 
    %plot(rho_p_sim, reshape(Exp_C_1_sim(2,2,:), 1, length(rho_p_sim)), 'r*', 'HandleVisibility','off');
    hold on, 
    plot(rho_p_sim, reshape(Exp_C_1_sim(2,3,:), 1, length(rho_p_sim)), 'r*', 'HandleVisibility','off');
    hold on,
    plot(rho_p_sim, reshape(Exp_C_1_sim(3,1,:), 1, length(rho_p_sim)), 'r*', 'HandleVisibility','off');
    %hold on,
    %plot(rho_p_sim, reshape(Exp_C_1_sim(3,2,:), 1, length(rho_p_sim)), 'r*', 'HandleVisibility','off');
    hold on,
    plot(rho_p_sim, reshape(Exp_C_1_sim(3,3,:), 1, length(rho_p_sim)), 'r*', 'HandleVisibility','off');    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves cosmetic makeover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid on;
%axis tight;
axis([min(rho_p_th) max(rho_p_th)...
    min(reshape(Exp_C_1_th(3,1,:), 1, length(rho_p_th))) * 0.0....
    max(reshape(Exp_C_1_th(1,3,:), 1, length(rho_p_th))) * 1.05]);
ylabel('$\Co$','FontSize',Fontsize+2);
xlabel('$\rho_p$ = [ms]','FontSize',Fontsize+2);
hl = legend('Theoretical', 'Simulated');
set(hl, 'Location', 'NorthEast', 'FontSize', Fontsize);
set(gca,'FontSize',Fontsize);
laprint(1, '../figures/fig_C_1_approx_vs_rho_p_diff_SNR', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');