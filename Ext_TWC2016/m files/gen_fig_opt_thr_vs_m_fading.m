%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      
%
%                                
%        Description: This m-file considers the varaition of optimim
%        secondary throughput (optimized over the sensing and estimation
%        time) against the m, the Nakagami-m parameter
%
%        
%        The simulation is performed to show the following:
%        1) Analyse the throughput vs sensing for different cases 
%           - Ideal case, when SNR is known at the ST (i)
%           - Outage Constraint   
%        2) Confirm the theoretical analysis vs simulation results
%
%        Created on: 28.01.16
%        Last modified: 28.01.16
%        Revision History: 28.01.16 --> File generated   
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

if 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves theoretical analysis --  Optimum throughput vs est time    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m_para = [0.8 1 1.5 2 3 4 5 10 50 100]
    R_id_opt_sen_opt_est_th = zeros(1, length(m_para));
    R_oc_opt_sen_opt_est_th = zeros(1, length(m_para));
    
    Fontsize = 8;
    figure(1);
    l = 0;
    if 1
        l = 1;
        %% m = 0.5
        load('results_opt_thr_vs_est_time_sen_time_fading_m_08.mat'); 
        [R_id_opt_sen_opt_est_th(l) index_opt] = max(R_id_opt_th(l,:));        
        [R_oc_opt_sen_opt_est_th(l) index_opt] = max(R_oc_opt_th(l,:));  
    end
        
   %% m = 1.5 1
   if 1
        l = l + 1; 
        load('results_opt_thr_vs_est_time_sen_time_fading.mat'); 
        [R_id_opt_sen_opt_est_th(l + 1) index_opt] = max(R_id_opt_th(1,:));        
        [R_id_opt_sen_opt_est_th(l) index_opt] = max(R_id_opt_th(2,:));  
        [R_oc_opt_sen_opt_est_th(l + 1) index_opt] = max(R_oc_opt_th(1,:));        
        [R_oc_opt_sen_opt_est_th(l) index_opt] = max(R_oc_opt_th(2,:));        
   end
   
   %% m = 2 3
    if 1
        l = l  + 1;
        load('results_opt_thr_vs_est_time_sen_time_fading_m_2_3.mat'); 
        for k=1:1:2
            [R_id_opt_sen_opt_est_th(k + l) index_opt] = max(R_id_opt_th(k,:));
            [R_oc_opt_sen_opt_est_th(k + l) index_opt] = max(R_oc_opt_th(k,:));        
        end
    end

   %% m = 4 5 10 50 100
    if 1
        l = l  + 2;
        load('results_opt_thr_vs_est_time_sen_time_fading_m_4_5_10_50_100.mat'); 
        for k=1:1:length(m_p2)
            [R_id_opt_sen_opt_est_th(k + l) index_opt] = max(R_id_opt_th(k,:));
            [R_oc_opt_sen_opt_est_th(k + l) index_opt] = max(R_oc_opt_th(k,:));        
        end
    end
    
    
    
    h1 = semilogx(m_para, R_id_opt_sen_opt_est_th, 'c', 'LineWidth',2);
    hold on,
    h2 = semilogx(m_para, R_oc_opt_sen_opt_est_th, 'LineStyle', '-', 'LineWidth', 1);

    hold on,
    h3 = semilogx(m_para(2), R_oc_opt_sen_opt_est_th(2), 'kx', 'LineWidth', 1);

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Plot curves cosmetic makeover
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    grid on;
    axis([min(m_para) max(m_para) R_oc_opt_sen_opt_est_th(1) max(R_oc_opt_sen_opt_est_th)* 1.05]);
    ylabel('$\rs(\ttest, \ttsen)$ [bits/sec/Hz]','FontSize',Fontsize);
    xlabel('Nakagami-$m$ parameter','FontSize',Fontsize);
    %hl = legend([h1 h2 h3 h5 h4],'$\rs$', '$\rsac$', '$\rsoc$', '$\ttsen$', 'Sim');
    hl = legend([h1 h2 h3],'IM, Problem 3', 'EM, Problem 4', 'Rayleigh');
    %set(hl, 'position',[0.725 0.12 0.15 0.31]);
    set(hl, 'position',[0.58 0.12 0.32 0.25]);
    set(gca,'FontSize',Fontsize);
    laprint(1, '../figures/fig_opt_thr_vs_m_fading', 'options', 'factory', 'width', 8, 'scalefonts',...
        'on', 'factor',0.5, 'keepfontprops', 'on');
end