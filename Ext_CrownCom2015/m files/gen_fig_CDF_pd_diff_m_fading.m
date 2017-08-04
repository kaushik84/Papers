%% This matlab script validates analytical expression for distribution function for the
%% pd (probability of detection) is charterized using Gaussian approximation for different 
%% value of estimation and sensing time. Nakagam -m fading added

clear all;
close all;

analysis = 0;                                                        % Enable it for performing a complete new analysis
                                                                     % otherwise a mat file is read 
                                                        
if analysis
        M = 1e4;                                                     % Number of realizations
        test = 100;                                                  % Estimation time  
        tsen = [100];                                                % Sensing time  
        noise_power = 10^(-100/10);                                  % Noise power        
        rcvd_energy = zeros(1, M);                                   % Received Energy     
        pdc = 0.90;                                                  % Constraint on probability of detection  
        sam_pts_th = 10000;                                          % Samples for probability of detection  
        snr_rcvd = 10.^([0]/10);                                     % snr for the channel g_p1          
        m_p = [0.5 1 2 5 100];                                       % Nakagam-m parameter
        kappa = 0.3;                                                 % Outage constraint on the             
        P_tran = 1;                                                  % Transmit power at the PR, PT 
        
        alpha_p2_true = 10^(-100/10);                                % Path Loss for the channel  PT - ST 
        rcvd_energy_bar = (P_tran * alpha_p2_true + noise_power);    % true value of received energy
        CDF_pd_sim = zeros(length(m_p), M);                          % Buffer to store the simulated CDF values 
        CDF_pd_ana = zeros(length(m_p), M);                          % Buffer to store the analytic CDF values 
        
        pd = zeros(1, M);
        pd_pts = zeros(length(m_p),M);                               % Buffer to store the C0 points for the distribution function
        
        
        
        for k = 1:length(m_p)
          %% Simulation
          %% Estimation of rcvd_energy is supposed to be non-central chi squared distributed 
           disp(strcat('m = ',num2str(m_p(k))));
           g_p2 = ones(1, M) .* random('gam',m_p(k), 1/m_p(k), 1, M); 
           for i=1:M
                rcvd_energy(i) = mean(random('norm',...
                    random('norm',0, sqrt(g_p2(i) * P_tran * alpha_p2_true),1, test),...
                    sqrt(noise_power), 1, test).^2);
            end

            % Determine the threshold  --- Actually this is a default pick           
            threshold = 2/tsen * ( noise_power + (rcvd_energy_bar  - noise_power) *....
                gammaincinv(kappa, m_p(k), 'upper')/m_p(k)) * gammaincinv(pdc, tsen/2 ,'upper');
          
            %threshold = gammaincinv(pdc, tsen/2,'upper') * 2 * rcvd_energy_bar / tsen;

            
           %% Realization of the random variable
            pd = gammainc( tsen/2 * threshold ./ rcvd_energy, tsen/2, 'upper');        

           %% CDF of pd    

           %% Emperical CDF 
            [CDF_pd_temp, pd_pts_temp,~,~,eid] = cdfcalc(pd);
            CDF_pd_temp = CDF_pd_temp(1:length(pd_pts_temp));
            pd_pts(k,:) = [pd_pts_temp; (pd_pts_temp(end) * ones(1, M - length(pd_pts_temp)))'];
            CDF_pd_sim(k,:) =  [CDF_pd_temp; (CDF_pd_temp(end) * ones(1, M - length(pd_pts_temp)))'];

           %% Anayltical CDF
           for i=1:M
               func_pd = @(t) gammainc( test * tsen * threshold ./...
                    (4 * (t * P_tran * alpha_p2_true + noise_power) *...
                    gammaincinv(pd_pts(k,i), tsen/2, 'upper')), test/2, 'upper') .*...
                    1/gamma(m_p(k)) * m_p(k)^(m_p(k)) .* t.^(m_p(k) - 1) .* exp(-m_p(k) * t);
               CDF_pd_ana(k,i) = 1 - integral(func_pd, 0 , 200);
           end
            
   
        end
        save('results_CDF_pd_diff_m_fading.mat');
        quit;
end
load('results_CDF_pd_diff_m_fading.mat');
        
figure(1);

% Plotting Curves
diff = 500;
sim = 1:diff:M;

k = 1;
plot(pd_pts(k, :), CDF_pd_ana(k, :), 'c-', 'Linewidth',1);
hold on,
plot(pd_pts(k, sim), CDF_pd_sim(k, sim), 'x', 'Linewidth',1);

hold on,
k = 2;
plot(pd_pts(k, :), CDF_pd_ana(k, :), 'c-', 'Linewidth',1);
hold on,
plot(pd_pts(k, sim), CDF_pd_sim(k, sim), 'x', 'Linewidth',1);

hold on,
k = 3;
plot(pd_pts(k, :), CDF_pd_ana(k, :), 'c-', 'Linewidth',1);
hold on,
plot(pd_pts(k, sim), CDF_pd_sim(k, sim), 'x', 'Linewidth',1);

hold on,
k = 4;
plot(pd_pts(k, :), CDF_pd_ana(k, :), 'c-', 'Linewidth',1);
hold on,
plot(pd_pts(k, sim), CDF_pd_sim(k, sim), 'x', 'Linewidth',1);


hold on,
k = 5;
plot(pd_pts(k, :), CDF_pd_ana(k, :), 'c-', 'Linewidth',1);
hold on,
plot(pd_pts(k, sim), CDF_pd_sim(k, sim), 'x', 'Linewidth',1);

Fontsize = 8;

grid on;
xlabel('$\pd$','FontSize',Fontsize);
ylabel('CDF','FontSize',Fontsize);
hl = legend('Theoretical', 'Simulated');
%set(hl, 'Location', 'NorthWest', 'FontSize', Fontsize);
set(hl, 'position',[0.14 0.77 0.28 0.14]);
set(gca,'FontSize',Fontsize);
laprint(1, '../figures/fig_CDF_pd_m_fading', 'options', 'factory', 'width', 7, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');        
