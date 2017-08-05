%% This matlab script validates analytical expression for distribution function for the
%% pd (probability of detection) is charterized using Gaussian approximation for different 
%% value of estimation and sensing time.

clear all;
close all;

analysis = 1;                                           % Enable it for performing a complete new analysis
                                                        % otherwise a mat file is read 
                                                        
if analysis
        M = 1e4;                                        % Number of realizations
        test = 5000;                                    % Estimation   
        tsen = 5000;                                    % Sensing time  
        noise_power = 10^(-100/10);                     % Noise power        
        rcvd_energy = zeros(1, M);                      % Received Energy     
        pdc = 0.90;                                     % Constraint on probability of detection  
        sam_pts_th = 10000;                             % Samples for probability of detection  
        snr_rcvd = 10.^([0]/10);                        % snr for the channel g_p1  
        
        rcvd_energy_bar = (snr_rcvd + 1) * noise_power; % true value of received energy
        
        pd = zeros(1, M);
        pd_pts = zeros(length(snr_rcvd),M);             % Buffer to store the C0 points for the distribution function
        
                
        for k = 1:length(tsen)
          %% Simulation
          %% Estimation of rcvd_energy is supposed to be non-central chi squared distributed 
            for i=1:M
                rcvd_energy(i) = mean(random('norm',...
                    random('norm',0, sqrt(snr_rcvd * noise_power),1, test),...
                    sqrt(noise_power), 1, test).^2);
            end

            tsen(k) = tsen(k) + test;
            % Determine the threshold 
            threshold = gammaincinv(pdc, tsen(k)/2, 'upper') * 2 * rcvd_energy_bar / tsen(k);
            
            %% Realization of the random variable
            pd = gammainc(tsen(k)/2 * threshold ./ rcvd_energy, tsen(k)/2, 'upper');              


            %% Emperical CDF 
            [CDF_pd_temp, pd_pts_temp,~,~,eid] = cdfcalc(pd);
            CDF_pd_temp = CDF_pd_temp(1:length(pd_pts_temp));
            pd_pts = [pd_pts_temp; (pd_pts_temp(end) * ones(1, M - length(pd_pts_temp)))'];
            CDF_pd_sim =  [CDF_pd_temp; (CDF_pd_temp(end) * ones(1, M - length(pd_pts_temp)))'];
            
            %% CDF Analytical
            cdf_pd_ana = 1 - gammainc( test * tsen * threshold./...
                (4 * rcvd_energy_bar * gammaincinv(pd_pts, tsen/2, 'upper')), test/2, 'upper'); 
            
            figure(1);
            plot(pd_pts, CDF_pd_sim, 'r');
            hold on,
            plot(pd_pts, cdf_pd_ana, 'b');
            
            
            %% pdf of pd is computed by plotting histogram  
            figure(2);
            bins = 1000;
            [f pd_pts] = hist(pd, bins);
            

            %% Anayltical PDF
            if 0 %% Numerical
                h = 1/bins;
                cdf_pd_ana = 1 - gammainc( test * tsen * threshold./...
                    (4 * rcvd_energy_bar * gammaincinv(pd_pts, tsen/2, 'upper')), test/2, 'upper');
                pdf_pd_ana = diff(cdf_pd_ana)/h;

            end
            pd_pts_ana = pd_pts; %linspace(1e-4,1 - 1e-4,bins);
            pdf_pd_ana = exp(gammaln(tsen/2) - gammaln(test/2) - test * tsen * threshold./...
                    (4 * rcvd_energy_bar * gammaincinv(pd_pts_ana, tsen/2, 'upper')) +...
                    gammaincinv(pd_pts_ana, tsen/2, 'upper') + test/2 * log(test * tsen * threshold./...
                    (4 * rcvd_energy_bar * gammaincinv(pd_pts_ana, tsen/2, 'upper'))) - tsen/2 *...
                    log(gammaincinv(pd_pts_ana, tsen/2, 'upper')));
                
            
            

            % Plotting curves
            bar(pd_pts,f/sum(f));        
            hold on,
            plot(pd_pts_ana, pdf_pd_ana * ((pd_pts_ana(3) - pd_pts_ana(2))), 'r', 'Linewidth',2.5);

            %plot(pd_pts(:,1:length(pdf_pd_ana)), pdf_pd_ana * ((pd_pts(3) - pd_pts(2))), 'r', 'Linewidth',2.5);
   
        end
        %save('results_CDF_pd_diff_tsen.mat');
end
%load('results_CDF_pd_diff_tsen.mat');
        
figure(1);

% Plotting Curves
diff = 10;
sim = 1:diff:M;

k = 1;
plot(pd_pts(k, :), CDF_pd_ana(k, :), 'c-', 'Linewidth',5);
hold on,
plot(pd_pts(k, sim), CDF_pd_sim(k, sim), '-', 'Linewidth',2);

hold on,
k = 2;
plot(pd_pts(k, :), CDF_pd_ana(k, :), 'c-', 'Linewidth',5);
hold on,
plot(pd_pts(k, sim), CDF_pd_sim(k, sim), '-', 'Linewidth',2);

hold on,
k = 3;
plot(pd_pts(k, :), CDF_pd_ana(k, :), 'c-', 'Linewidth',5);
hold on,
plot(pd_pts(k, sim), CDF_pd_sim(k, sim), '-', 'Linewidth',2);

Fontsize = 9;

grid on;
xlabel('$\pd$','FontSize',Fontsize);
ylabel('CDF','FontSize',Fontsize);
hl = legend('Theoretical', 'Simulated');
%set(hl, 'Location', 'NorthWest', 'FontSize', Fontsize);
set(hl, 'position',[0.14 0.77 0.34 0.14]);
set(gca,'FontSize',Fontsize);
laprint(1, '../figures/fig_CDF_pd_diff_tsen', 'options', 'factory', 'width', 7, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');        
