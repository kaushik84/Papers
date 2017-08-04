%% This matlab script validates analytical expression for density function for the
%% pd (probability of detection) is charterized using Gaussian approximation for different 
%% value of estimation and sensing time. Nakagam -m fading added

clear all;
close all;

analysis = 0;                                                        % Enable it for performing a complete new analysis
                                                                     % otherwise a mat file is read 
                                                        
if analysis
        M = 5e3;                                                     % Number of realizations
        test = 100;                                                  % Estimation time  
        tsen = [150];                                                % Sensing time  
        noise_power = 10^(-100/10);                                  % Noise power        
        rcvd_energy = zeros(1, M);                                   % Received Energy     
        pdc = 0.90;                                                  % Constraint on probability of detection  
        m_p = 10;%[0.5 1 2 5 100];                                       % Nakagam-m parameter
        kappa = 0.3;                                                 % Outage constraint on the             
        P_tran = 1;                                                  % Transmit power at the PR, PT 
        
        alpha_p2_true = 10^(-100/10);                                % Path Loss for the channel  PT - ST 
        rcvd_energy_bar = (P_tran * alpha_p2_true + noise_power);    % True value of received energy 
        bins = 200;                                                   
        PDF_pd_ana = zeros(length(m_p), bins);                       % Probability density function
        
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
            threshold = 2/tsen * (noise_power + (rcvd_energy_bar - noise_power) *...
                gammaincinv(kappa, m_p(k), 'upper')/m_p(k)) * gammaincinv(pdc, tsen/2 ,'upper');
            
           %% Emperical PDF
            pd = gammainc( tsen/2 * threshold ./ rcvd_energy, tsen/2, 'upper');    
            [f x_pts] = hist(pd, bins);   

           %% Anayltical PDF
           for i=1:length(x_pts)
                func_pd = @(t) exp(gammaln(tsen/2) - gammaln(test/2)) *...  
                    exp(-threshold * test * tsen./ (4 * (t * P_tran * alpha_p2_true + noise_power) *...
                   gammaincinv(x_pts(i),tsen/2,'upper')) + gammaincinv(x_pts(i),tsen/2,'upper')) .*...
                   (threshold * tsen * test./ (4 * (t * P_tran * alpha_p2_true + noise_power) *...
                   gammaincinv(x_pts(i),tsen/2,'upper'))).^(test/2) .*...
                    (1/gammaincinv(x_pts(i),tsen/2,'upper'))^(tsen/2) .* ... 
                   1/gamma(m_p(k)) * m_p(k)^(m_p(k)) .* t.^(m_p(k) - 1) .* exp(-m_p(k) * t);
                PDF_pd_ana(k,i) = integral(func_pd, 0 , 100);   
               %               func_pd = @(t) gammainc( test * tsen * threshold ./...
%                    (4 * (t * P_tran * alpha_p2_true + noise_power) *...
%                    gammaincinv(pd_pts(k,i), tsen/2, 'upper')), test/2, 'upper') .*...
%                    1/gamma(m_p(k)) * m_p(k)^(m_p(k)) .* t.^(m_p(k) - 1) .* exp(-m_p(k) * t);
%               CDF_pd_ana(k,i) = integral(func_pd, 0 , 200);
           end
            
           % Plotting curves
            bar(x_pts,f/sum(f));    
            hold on,
            plot(x_pts, PDF_pd_ana(k,:) * ((x_pts(3) - x_pts(2))), 'r', 'Linewidth',2.5);
            grid on,
        end
        save('results_CDF_pd_diff_m_fading.mat');
        %quit;
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
