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


spacing = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    load('results_opt_thr_vs_snr_AWGN_wo_nu_th_CR.mat');
    h1 = plot(-alpha_true_p - noise_power, Exp_R_opt(1,:),'r-', 'LineWidth',1);
    hold on,
    plot(-alpha_true_p(1:spacing:length(alpha_true_p)) - noise_power,...
        Exp_R_opt(1,1:spacing:length(alpha_true_p)), 'ro');
    hold on,
    h2 = plot(-alpha_true_p - noise_power, Exp_R_opt(2,:), 'LineWidth',1, 'HandleVisibility','off');
    hold on,
    plot(-alpha_true_p(1:spacing:length(alpha_true_p)) - noise_power,...
        Exp_R_opt(2,1:spacing:length(alpha_true_p)), 's');
    hold on,
    h3 = plot(-alpha_true_p - noise_power, Exp_R_opt(3,:), 'm-','LineWidth',1, 'HandleVisibility','off');
    hold on,
    plot(-alpha_true_p(1:spacing:length(alpha_true_p)) - noise_power,...
        Exp_R_opt(3,1:spacing:length(alpha_true_p)), 'mx');
    Y_axis_min =  min(Y_axis_min, min(min(Exp_R_opt)));
    Y_axis_max =  max(Y_axis_max, max(max(Exp_R_opt)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Exp_R_conventional = log2(1 + 10./(10.^((P_tran - alpha_true_p- noise_power) /10)));
hold on,
h4 = plot(-alpha_true_p - noise_power, Exp_R_conventional, 'k-', 'LineWidth', 1);
Y_axis_min =  min(Y_axis_min, min(min(Exp_R_opt)));
Y_axis_max =  max(Y_axis_max, max(max(Exp_R_conventional)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis with nu upper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    load('results_opt_thr_vs_snr_AWGN_nu_upper_th_CR.mat');
    hold on,
    plot( - alpha_true_p - noise_power_p  + nu, Exp_R_opt, 'LineWidth',1, 'LineStyle','-.', 'HandleVisibility','off');

    Y_axis_min =  min(Y_axis_min, min(Exp_R_opt));
    Y_axis_max =  max(Y_axis_max, max(Exp_R_opt));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis with nu lower
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    load('results_opt_thr_vs_snr_AWGN_nu_lower_th_CR.mat');
    hold on,
    h2  = plot(- alpha_true_p - noise_power_p  + nu, Exp_R_opt, 'LineWidth',1, 'LineStyle',':');


    Y_axis_min =  min(Y_axis_min, min(Exp_R_opt));
    Y_axis_max =  max(Y_axis_max, max(Exp_R_opt));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

presentation = 1; %% Enable this flag to plot figure for the presentation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves cosmetic makeover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~presentation
    Fontsize = 9;
    grid on;
    axis([-15 10 0.98 8.3]);
    ylabel('$\ers$ = [bits/sec/Hz]','FontSize',Fontsize);
    xlabel('$\gamma$ = [dB]','FontSize',Fontsize);
    hl = legend([h4 h1 h2 h3],'(6)', '$\pcd = 0.92$', '$\pcd = 0.95$', '$\pcd = 0.97$');
    c = get(hl,'Children');
    set(c(7),'Marker','o');
    set(c(4),'Marker','s'); 
    set(c(1),'Marker','x');
    set(hl, 'Location', 'NorthEast', 'FontSize', Fontsize);
    set(gca,'FontSize',Fontsize);
    if 0
        set(gcf, 'PaperUnits','inches');
        set(gcf, 'PaperSize',[10 7.5]);
        set(gcf, 'PaperPositionMode','manual');
        set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
        print(gcf,'-dpdf','../figures/fig_opt_thr_vs_snr_AWGN');
        print(gcf,'-depsc','../figures/fig_opt_thr_vs_snr_AWGN');
    end
    laprint(1, '../figures/fig_opt_thr_vs_snr_AWGN', 'options', 'factory', 'width', 7, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
else
    Fontsize = 5.5;
    grid on;
    axis([-15 10 0.98 8.3]);
    ylabel('$\ers$ = [bits/sec/Hz]','FontSize',Fontsize);
    xlabel('$\gamma$ = [dB]','FontSize',Fontsize);
    hl = legend([h4 h1 h2 h3],'(6)', '$\pcd = 0.92$', '$\pcd = 0.95$', '$\pcd = 0.97$');
    c = get(hl,'Children');
    set(c(7),'Marker','o');
    set(c(4),'Marker','s'); 
    set(c(1),'Marker','x');
    set(hl, 'Location', 'NorthEast', 'FontSize', Fontsize);
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
        'on', 'factor',0.5, 'keepfontprops', 'on');
end