%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: This m-file considers the Optimum throughput-accuracy 
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
X_axis_min = inf;
X_axis_max = 0;
spacing  = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_opt_thr_vs_acc_AWGN_wo_nu_th.mat');
temp = 6:1:length(Exp_R_opt)/2 + 4;
temp2 = 6:spacing:length(Exp_R_opt)/2 + 4;
h1 = plot(mu_P_th(temp), Exp_R_opt(1,(temp)),'r-', 'LineWidth',1);
hold on,
plot(mu_P_th(temp2), Exp_R_opt(1,(temp2)), 'ro');
hold on,
h2 = plot(mu_P_th(temp), Exp_R_opt(2,(temp)), 'LineWidth',1, 'HandleVisibility','off');
hold on,
plot(mu_P_th(temp2), Exp_R_opt(2,(temp2)), 'bs');
hold on,
h3 = plot(mu_P_th(temp), Exp_R_opt(3,(temp)),'m-', 'LineWidth',1, 'HandleVisibility','off');
hold on,
plot(mu_P_th(temp2), Exp_R_opt(3,(temp2)), 'mx');

Exp_R_conventional = log2(1 + 10) * ones(1, length(mu_P_th(temp)));
hold on,
h4 = plot(mu_P_th(temp), Exp_R_conventional, 'k-', 'LineWidth', 1);

Y_axis_min =  min(Y_axis_min, min(min(Exp_R_opt(:,temp))));
Y_axis_max =  max(Y_axis_max, max(Exp_R_conventional));
X_axis_min =  min(Y_axis_min, min(min(mu_P_th(:,temp))));
X_axis_max =  max(Y_axis_max, max(max(mu_P_th(:,temp))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis with nu upper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('results_opt_thr_vs_acc_AWGN_nu_u_th.mat');
temp = 6:1:length(Exp_R_opt)/2 + 4;
%hold on,
%h2 = plot(mu_P_th(temp), Exp_R_opt(temp), 'LineWidth',1.5, 'LineStyle','-.');

Y_axis_min =  min(Y_axis_min, min(Exp_R_opt(temp)));
Y_axis_max =  max(Y_axis_max, max(Exp_R_opt(temp)));
X_axis_min =  min(Y_axis_min, min(mu_P_th(temp)));
X_axis_max =  max(Y_axis_max, max(mu_P_th(temp)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis with nu lower
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results_opt_thr_vs_acc_AWGN_nu_l_th.mat');
temp = 6:1:length(Exp_R_opt)/2 + 4;
%hold on,
%h3 = plot(mu_P_th(temp), Exp_R_opt(temp), 'LineWidth',1.5, 'LineStyle','--');


Y_axis_min =  min(Y_axis_min, min(Exp_R_opt(temp)));
Y_axis_max =  max(Y_axis_max, max(Exp_R_opt(temp)));
X_axis_min =  min(Y_axis_min, min(mu_P_th(temp)));
X_axis_max =  max(Y_axis_max, max(mu_P_th(temp)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

presentation = 1; %% Enable this flag to plot figure for the presentation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves cosmetic makeover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~presentation
        Fontsize = 9;
        grid on;
        axis([min(mu_P_th(temp)) max(mu_P_th(temp)) Y_axis_min Y_axis_max*1.02]);
        ylabel('$\ers$ = [bits/sec/Hz]','FontSize',Fontsize);
        xlabel('Accuracy ($\mu$)','FontSize',Fontsize);
        hl = legend([h4 h1 h2 h3],'(6)','$\pcd = 0.92$','$\pcd = 0.95$','$\pcd = 0.97$');
        set(hl, 'Location', 'SouthEast', 'FontSize', Fontsize);
        c = get(hl,'Children');
        set(c(7),'Marker','o');
        set(c(4),'Marker','s'); 
        set(c(1),'Marker','x');
        if 0
            set(gca,'FontSize',Fontsize);
            set(gcf, 'PaperUnits','inches');
            set(gcf, 'PaperSize',[10 7.5]);
            set(gcf, 'PaperPositionMode','manual');
            set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
            print(gcf,'-dpdf','../figures/fig_opt_thr_vs_acc_AWGN');
            print(gcf,'-depsc','../figures/fig_opt_thr_vs_acc_AWGN');
        end
        laprint(1, '../figures/fig_opt_thr_vs_acc_AWGN', 'options', 'factory', 'width', 7, 'scalefonts',...
            'on', 'factor',0.5, 'keepfontprops', 'on');
else
    Fontsize = 5.5;
    grid on;
    axis([min(mu_P_th(temp)) max(mu_P_th(temp)) Y_axis_min Y_axis_max*1.02]);
    ylabel('$\ers$ = [bits/sec/Hz]','FontSize',Fontsize);
    xlabel('Accuracy ($\mu$)','FontSize',Fontsize);
    hl = legend([h4 h1 h2 h3],'(6)','$\pcd = 0.92$','$\pcd = 0.95$','$\pcd = 0.97$');
    set(hl, 'Location', 'SouthEast', 'FontSize', Fontsize);
    c = get(hl,'Children');
    set(c(7),'Marker','o');
    set(c(4),'Marker','s'); 
    set(c(1),'Marker','x');
    if 0
        set(gca,'FontSize',Fontsize);
        set(gcf, 'PaperUnits','inches');
        set(gcf, 'PaperSize',[10 7.5]);
        set(gcf, 'PaperPositionMode','manual');
        set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
        print(gcf,'-dpdf','../figures/fig_opt_thr_vs_acc_AWGN');
        print(gcf,'-depsc','../figures/fig_opt_thr_vs_acc_AWGN');
    end
    laprint(1, '../figures/fig_opt_thr_vs_acc_AWGN_pres', 'options', 'factory', 'width', 6, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
end