%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: This m-file considers the Optimum throughput-SNR 
%        for the AWGN channel.        
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
%        Created on: 19.09.14
%        Revision History: 19.09.14 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

Y_axis_min = inf;
Y_axis_max = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_opt_thr_vs_snr_AWGN_wo_nu_th.mat');
h1 = plot(P_tran, Exp_R_opt(1,:), 'LineWidth',1.5);
hold on,
plot(P_tran, Exp_R_opt(2,:), 'LineWidth',1.5, 'HandleVisibility','off');
hold on,
plot(P_tran, Exp_R_opt(3,:), 'LineWidth',1.5, 'HandleVisibility','off');

Y_axis_min =  min(Y_axis_min, min(min(Exp_R_opt)));
Y_axis_max =  max(Y_axis_max, max(max(Exp_R_opt)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis with nu upper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_opt_thr_vs_snr_AWGN_nu_u_th.mat');
hold on,
plot(P_tran, Exp_R_opt, 'LineWidth',1.5, 'LineStyle','-.', 'HandleVisibility','off');

Y_axis_min =  min(Y_axis_min, min(Exp_R_opt));
Y_axis_max =  max(Y_axis_max, max(Exp_R_opt));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis with nu lower
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_opt_thr_vs_snr_AWGN_nu_l_th.mat');
hold on,
h2  = plot(P_tran, Exp_R_opt, 'LineWidth',1.5, 'LineStyle','-.');


Y_axis_min =  min(Y_axis_min, min(Exp_R_opt));
Y_axis_max =  max(Y_axis_max, max(Exp_R_opt));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

presentation = 1; %% Enable this flag to plot figure for the presentation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves cosmetic makeover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~presentation
    Fontsize = 20;
    grid on;
    axis([min(P_tran) max(P_tran) Y_axis_min Y_axis_max]);
    ylabel('$\ers$ = [bits/sec/Hz]','FontSize',Fontsize+2);
    xlabel('$\snrrcvd$ = [dB]','FontSize',Fontsize+2);
    hl = legend([h1;h2],'$\npu$ = 0 dB', '$\npu = \pm$ 3 dB');
    set(hl, 'Location', 'SouthEast', 'FontSize', Fontsize + 2);
    set(gca,'FontSize',Fontsize);
    if 0
        set(gcf, 'PaperUnits','inches');
        set(gcf, 'PaperSize',[10 7.5]);
        set(gcf, 'PaperPositionMode','manual');
        set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
        print(gcf,'-dpdf','../figures/fig_opt_thr_vs_snr_AWGN');
        print(gcf,'-depsc','../figures/fig_opt_thr_vs_snr_AWGN');
    end
    laprint(1, '../figures/fig_opt_thr_vs_snr_AWGN', 'options', 'factory', 'width', 16, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'off');
else
    Fontsize = 7;
    grid on;
    axis([min(P_tran) max(P_tran) Y_axis_min Y_axis_max]);
    ylabel('$\ers$ = [bits/sec/Hz]','FontSize',Fontsize+2);
    xlabel('$\snrrcvd$ = [dB]','FontSize',Fontsize+2);
    hl = legend([h1;h2],'$\npu$ = 0 dB', '$\npu = \pm$ 3 dB');
    set(hl, 'Location', 'SouthEast', 'FontSize', Fontsize + 2);
    set(gca,'FontSize',Fontsize);
    if 0
        set(gcf, 'PaperUnits','inches');
        set(gcf, 'PaperSize',[10 7.5]);
        set(gcf, 'PaperPositionMode','manual');
        set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
        print(gcf,'-dpdf','../figures/fig_opt_thr_vs_snr_AWGN');
        print(gcf,'-depsc','../figures/fig_opt_thr_vs_snr_AWGN');
    end
    laprint(1, '../figures/fig_opt_thr_vs_snr_AWGN_pres', 'options', 'factory', 'width', 6, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'off');
end