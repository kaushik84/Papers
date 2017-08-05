%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: This m-file considers the Throughput-Estimation time 
%        tradeoff for the fading channel.        
%
%        
%        The simulation is performed to show the following:
%        1) Analyse the throughput vs estimation curves for diff snr
%
%        Created on: 17.09.14
%        Revision History: 17.09.14 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

Y1_axis_min = 0;
Y1_axis_max = 0;
Y2_axis_min = 0;
Y2_axis_max = 0;


g1 = random('exp', 1 , 1, 1e6);
g2 = random('exp', 1 , 1, 1e6);
Exp_R_conventional = mean(log2(1 + g1 * 10./g2));

temp = 3:1:111;
temp2 = [4,15,50,93,100,104,111];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis wo nu snr m00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    load('results_thr_est_time_tradeoff_fading_wo_nu_snr_m00_sim.mat');
    
    % Some values for the PC_P_p_th resulted Inf due to large N, hence
    % these values have already attained a saturation value, hence a
    % workround to replace these Inf values with the last non Inf value.
    nonInf_index = 0;
    for i = temp
        if  PC_P_p_th(i) == Inf
            if nonInf_index == 0
                nonInf_index = i - 1;
            end            
            PC_P_p_th(i) = PC_P_p_th(nonInf_index);
        end
    end
    
    [ax h1 h2] = plotyy(N_sim(temp) * 1e-3, Exp_R_sim(temp), N_sim(temp) * 1e-3, PC_P_p_th(temp));
    set(h1,'LineWidth',1);
    set(h2,'LineWidth',1);
    set(ax,'XScale','log');
 
    %Adding simulation result to represented as theory --> only for display
    if 1
        set(ax,'NextPlot','add');
        h5 = plot(ax(2), N_sim(temp2) * 1e-3, PC_P_p_sim(temp2), 'r*', 'HandleVisibility','off');
        set(ax,'XScale','log');
    end
    Y1_axis_min =  min(Exp_R_sim(temp));
    Y1_axis_max =  max(Exp_R_sim(temp));
    Y2_axis_min =  min(PC_P_p_sim(temp));
    Y2_axis_max =  max(PC_P_p_sim(temp));
    
    %   Plot curves for the convetional model
    set(ax,'NextPlot','add');
    h6 = plot(ax(1), N_sim(temp) * 1e-3, Exp_R_conventional * ones(1, length(N_sim(temp))),...
    'k-', 'LineWidth',1);
    set(ax,'XScale','log');
    Y1_axis_max =  max(Y1_axis_max, max(Exp_R_conventional+0.02));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis wo nu snr m10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    load('results_thr_est_time_tradeoff_fading_wo_nu_snr_m10_sim.mat');
    % Some values for the PC_P_p_th resulted Inf due to large N, hence
    % these values have already attained a saturation value, hence a
    % workround to replace these Inf values with the last non Inf value.
    nonInf_index = 0;
    for i = temp
        if  PC_P_p_th(i) == Inf
            if nonInf_index == 0
                nonInf_index = i - 1;
            end            
            PC_P_p_th(i) = PC_P_p_th(nonInf_index);
        end
    end
    set(ax,'NextPlot','add');
    h3 = plot(ax(1), N_sim(temp) * 1e-3, Exp_R_sim(temp), 'Linewidth', 1.5, 'LineStyle', '-.');
    set(ax,'XScale','log');
    set(ax,'NextPlot','add');
    h4 = plot(ax(2), N_sim(temp) * 1e-3, PC_P_p_th(temp), 'Linewidth', 1.5, 'LineStyle', '-.');
    set(ax,'XScale','log');
    
    %Adding simulation result to represented as theory --> only for display
    if 1
        set(ax,'NextPlot','add');
        plot(ax(2), N_sim(temp2) * 1e-3, PC_P_p_sim(temp2), 'r*');
        set(ax,'XScale','log');
    end
    
    Y1_axis_min =  min(Y1_axis_min, min(Exp_R_sim(temp)));
    Y1_axis_max =  max(Y1_axis_max, max(Exp_R_sim(temp)));
    Y2_axis_min =  min(Y2_axis_min, min(PC_P_p_th(temp)));
    Y2_axis_max =  max(Y2_axis_max, max(PC_P_p_th(temp)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot curves theoretical analysis wo nu snr p10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    load('results_thr_est_time_tradeoff_fading_wo_nu_snr_p10_sim.mat');
    set(ax,'NextPlot','add');
    plot(ax(1), N_sim * 1e-3, Exp_R_sim, 'Linewidth', 1.5);
    set(ax,'NextPlot','add');
    plot(ax(2), N_sim * 1e-3, PC_P_p_th, 'Linewidth', 1.5);

    Y1_axis_min =  min(Y1_axis_min, min(Exp_R_sim));
    Y1_axis_max =  max(Y1_axis_max, max(Exp_R_sim));
    Y2_axis_min =  min(Y2_axis_min, min(PC_P_p_th));
    Y2_axis_max =  max(Y2_axis_max, max(PC_P_p_th));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1

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
    axis(ax(1), [min(N_sim(temp) * 1e-3) max(N_sim(temp) * 1e-3) Y1_axis_min Y1_axis_max + 0.02]);
    axis(ax(2), [min(N_sim(temp) * 1e-3) max(N_sim(temp) * 1e-3) Y2_axis_min 0.032]);
    ylabel(ax(1),'$\ers$ = [bits/sec/Hz]','FontSize',Fontsize);
    ylabel(ax(2),'$\pc$','FontSize',Fontsize);
    xlabel(ax(1),'$\tau$ = [ms]','FontSize',Fontsize);
    xlabel(ax(2), '$\tau$ = [ms]','FontSize',Fontsize); 
    %hl = legend([h6;h1;h3;h2;h4;h5],'(6)', '\gamma = -10 dB   ', '\gamma = 0 dB   ', '\gamma = -10 dB   ',  '\gamma = 0 dB',...
    %    'theory');
    hl = legend([h1 h2 h6 h5],'$\ers$', '$\pc$', '(6)', 'theory');
    set(hl, 'Position', [.58, .19, 0.2, 0.2], 'FontSize', Fontsize);
    set(ax(1),'FontSize',Fontsize);
    set(ax(2),'FontSize',Fontsize);
    Ytick_R =  round(linspace(Y1_axis_min, Y1_axis_max, 6)*100)/100;
    Ytick_PC = round(linspace(Y2_axis_min, .032, 6)*1000)/1000;
    set(ax(1), 'YTick', Ytick_R);
    set(ax(2), 'YTick', Ytick_PC);
    if 0
        set(gcf, 'PaperUnits','inches');
        set(gcf, 'PaperSize',[11 7.5]);
        set(gcf, 'PaperPositionMode','manual');
        set(gcf, 'PaperPosition',[ 0 0 11 7.5]);
        print(gcf,'-dpdf','../figures/fig_thr_est_time_tradeoff_fading');
        print(gcf,'-depsc','../figures/fig_thr_est_time_tradeoff_fading');
    end
    laprint(1, '../figures/fig_thr_est_time_tradeoff_fading', 'options', 'factory', 'width', 7, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
else
    Fontsize = 5.5;
    grid on;
    axis(ax(1), [min(N_sim(temp) * 1e-3) max(N_sim(temp) * 1e-3) Y1_axis_min Y1_axis_max + 0.02]);
    axis(ax(2), [min(N_sim(temp) * 1e-3) max(N_sim(temp) * 1e-3) Y2_axis_min 0.032]);
    ylabel(ax(1),'$\ers$ = [bits/sec/Hz]','FontSize',Fontsize);
    ylabel(ax(2),'$\pc$','FontSize',Fontsize);
    xlabel(ax(1),'$\tau$ = [ms]','FontSize',Fontsize);
    xlabel(ax(2), '$\tau$ = [ms]','FontSize',Fontsize); 
    %hl = legend([h6;h1;h3;h2;h4;h5],'(6)', '\gamma = -10 dB   ', '\gamma = 0 dB   ', '\gamma = -10 dB   ',  '\gamma = 0 dB',...
    %    'theory');
    hl = legend([h1 h2 h6 h5],'$\ers$', '$\pc$', '(6)', 'theory');
    set(hl, 'Position', [.58, .19, 0.2, 0.2], 'FontSize', Fontsize);
    set(ax(1),'FontSize',Fontsize);
    set(ax(2),'FontSize',Fontsize);
    Ytick_R =  round(linspace(Y1_axis_min, Y1_axis_max, 6)*100)/100;
    Ytick_PC = round(linspace(Y2_axis_min, .032, 6)*1000)/1000;
    set(ax(1), 'YTick', Ytick_R);
    set(ax(2), 'YTick', Ytick_PC);
    if 0
        set(gcf, 'PaperUnits','inches');
        set(gcf, 'PaperSize',[11 7.5]);
        set(gcf, 'PaperPositionMode','manual');
        set(gcf, 'PaperPosition',[ 0 0 11 7.5]);
        print(gcf,'-dpdf','../figures/fig_thr_est_time_tradeoff_fading');
        print(gcf,'-depsc','../figures/fig_thr_est_time_tradeoff_fading');
    end
    laprint(1, '../figures/fig_thr_est_time_tradeoff_fading_pres', 'options', 'factory', 'width', 6, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
end