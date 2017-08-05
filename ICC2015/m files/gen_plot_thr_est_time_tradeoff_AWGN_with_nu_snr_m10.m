%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: This m-file considers the Throughput-Estimation time 
%        tradeoff for the AWGN channel.        
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
%        Created on: 17.09.14
%        Revision History: 17.09.14 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

Fontsize = 18;
Y1_axis_min = 0;
Y1_axis_max = 0;
Y2_axis_min = 0;
Y2_axis_max = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis wo nu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_thr_est_time_tradeoff_AWGN_wo_nu_snr_m10_th.mat');
PC_opt = PC_P_p_opt;
[ax h1 h2] = plotyy(N_th * 1e-3, Exp_R_th, N_th * 1e-3, PC_P_p_th);
set(h1,'LineWidth',1.5);
set(h2,'LineWidth',1.5);
%set(ax,'NextPlot','add');
%plot(ax(2), [N_optimum,max(N_opt)] * 1e-3, PC_opt * ones(1, 2),'LineStyle','--','LineWidth',1.5);
set(ax,'NextPlot','add');
plot(ax(2), N_optimum * ones(1, 2) * 1e-3, [PC_opt, min(PC_P_p_th)],...
    'LineWidth',1.5,'LineStyle','--', 'HandleVisibility','off');
set(ax,'NextPlot','add');
plot(ax(1), N_optimum * 1e-3, Exp_R_opt, 'o','LineWidth',1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves simulation wo nu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_thr_est_time_tradeoff_AWGN_wo_nu_snr_m10_sim.mat');
set(ax,'NextPlot','add');
plot(ax(1), N_sim * 1e-3, Exp_R_sim, 'r*', 'HandleVisibility','off');
set(ax,'NextPlot','add');
plot(ax(2), N_sim * 1e-3, PC_P_p_sim, 'r*', 'HandleVisibility','off');

Y1_axis_min =  min(Exp_R_th);
Y1_axis_max =  max(Exp_R_th);
Y2_axis_min =  min(PC_P_p_th);
Y2_axis_max =  max(PC_P_p_th);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis with nu upper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_thr_est_time_tradeoff_AWGN_nu_u_snr_m10_th.mat');
PC_opt = PC_P_p_opt;
set(ax,'NextPlot','add');
plot(ax(1), N_th * 1e-3, Exp_R_th, 'LineWidth',1.5, 'LineStyle','-.');
set(ax,'NextPlot','add');
h3 = plot(ax(2), N_th * 1e-3, PC_P_p_th, 'LineWidth',1.5, 'LineStyle','-.');
%set(ax,'NextPlot','add');
%plot(ax(2), [N_optimum,max(N_opt)] * 1e-3, PC_opt * ones(1, 2),'LineStyle','--','LineWidth',1.5);
set(ax,'NextPlot','add');
plot(ax(2), N_optimum * ones(1, 2) * 1e-3, [PC_opt, min(PC_R_th)],...
    'LineWidth',1.5,'LineStyle','--', 'HandleVisibility','off');
set(ax,'NextPlot','add');
plot(ax(1), N_optimum * 1e-3, Exp_R_opt, 'o','LineWidth',1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves simulation with nu upper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_thr_est_time_tradeoff_AWGN_nu_u_snr_m10_sim.mat');
set(ax,'NextPlot','add');
plot(ax(1), N_sim * 1e-3, Exp_R_sim, 'r*', 'HandleVisibility','off');
set(ax,'NextPlot','add');
plot(ax(2), N_sim * 1e-3, PC_P_p_sim, 'r*', 'HandleVisibility','off');

Y1_axis_min =  min(Y1_axis_min, min(Exp_R_th));
Y1_axis_max =  max(Y1_axis_max, max(Exp_R_th));
Y2_axis_min =  min(Y1_axis_min, min(PC_P_p_th));
Y2_axis_max =  max(Y1_axis_max, max(PC_P_p_th));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis with nu lower
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_thr_est_time_tradeoff_AWGN_nu_l_snr_m10_th.mat');
PC_opt = PC_P_p_opt;
set(ax,'NextPlot','add');
plot(ax(1), N_th * 1e-3, Exp_R_th, 'LineWidth',1.5, 'LineStyle','-.');
set(ax,'NextPlot','add');
plot(ax(2), N_th * 1e-3, PC_P_p_th, 'LineWidth',1.5, 'LineStyle','-.');
set(ax,'NextPlot','add');
h4 = plot(ax(2), [N_optimum,max(N_opt)] * 1e-3, PC_opt * ones(1, 2),'LineStyle','--','LineWidth',1.5);
set(ax,'NextPlot','add');
plot(ax(2), N_optimum * ones(1, 2) * 1e-3, [PC_opt, min(PC_P_p_th)],...
    'LineWidth',1.5,'LineStyle','--', 'HandleVisibility','off');
set(ax,'NextPlot','add');
h5 = plot(ax(1), N_optimum * 1e-3, Exp_R_opt, 'o','LineWidth',1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves simulation with nu lower
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_thr_est_time_tradeoff_AWGN_nu_l_snr_m10_sim.mat');
set(ax,'NextPlot','add');
h6  = plot(ax(1), N_sim * 1e-3, Exp_R_sim, 'r*', 'HandleVisibility','off');
set(ax,'NextPlot','add');
plot(ax(2), N_sim * 1e-3, PC_P_p_sim, 'r*', 'HandleVisibility','off');

Y1_axis_min =  min(Y1_axis_min, min(Exp_R_th));
Y1_axis_max =  max(Y1_axis_max, max(Exp_R_th));
Y2_axis_min =  min(Y1_axis_min, min(PC_P_p_th));
Y2_axis_max =  max(Y1_axis_max, max(PC_P_p_th));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves cosmetic makeover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid on;
axis(ax(1), [min(N_th * 1e-3) max(N_th * 1e-3) Y1_axis_min Y1_axis_max]);
axis(ax(2), [min(N_th * 1e-3) max(N_th * 1e-3) Y2_axis_min 1]);
ylabel(ax(1),'R_s = [bits/sec/Hz]','FontSize',Fontsize+2);
ylabel(ax(2),'Probability','FontSize',Fontsize+2);
xlabel(ax(1),'\tau = [ms]','FontSize',Fontsize+2);
xlabel(ax(2), '\tau = [ms]','FontSize',Fontsize+2);
hl = legend([h1;h2;h3;h4;h5;h6],'R', 'PC', 'PC with nu', 'PC = \epsilon', 'max \tau', 'sim');
set(hl, 'Position', [.7, .5, .1, 0.2], 'FontSize', Fontsize + 2);
set(ax(1),'FontSize',Fontsize);
set(ax(2),'FontSize',Fontsize);
Ytick_R =  round(linspace(Y1_axis_min, Y1_axis_max, 10)*100)/100;
Ytick_PC = round(linspace(Y2_axis_min, 1, 10)*100)/100;
set(ax(1), 'YTick', Ytick_R);
set(ax(2), 'YTick', Ytick_PC);
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize',[10 7.5]);
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
print(gcf,'-dpdf','../figures/fig_thr_est_time_tradeoff_AWGN_snr_m10');
print(gcf,'-depsc','../figures/fig_thr_est_time_tradeoff_AWGN_snr_m10');