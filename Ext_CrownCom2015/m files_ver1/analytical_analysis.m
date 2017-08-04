%% Hereby, we would like to investigate 

%% 1) the distortion in the probability of detection. 
%% In this regard, we would like to characterize the distribution
%% function and/or the probability density functions P_rcvd.

%% 2) determine the performance of the system by capturing the expected value of the 
%% capacities C0 = log2(1 + snr_s) and C1 =  log2(1 + snr_s / (snr_p  + 1)).  Channel Estimation
%% is performed (Case 1: Matched filtering is employed for both the channels (PT-SR, ST-SR), 
%% Case 2: Here Matched filtering is peformed on ST-SR and energy detection is performed on PT-SR), 
%% hence we first characterize the densities of snr_s and snr_p and thereafter
%% characterieze the densities of C_0 and C_1
clear all;
close all;
clc;

%% 1) Analysis

if 0
   if 1 % Short term analysis
        % Case where the energy measurement consisiting of K samples and where K samples sees N 
        % different realization of the channel and 1 < N < K
        M = 5e4;                    % Num of realizations
        t_est =  1000;               % Num of samples used for estimation  
        t_sen =  1000;               % Num of samples used for sensing  

        noise_power = 1e-10;
        x = zeros(1,M);
        SNR = +50;
        snr = 10^(SNR/10);
        signal_power = noise_power * 10^(SNR/10);    
        P_d_d = 0.9;
        energy = zeros(1,M);
        %% Lets first demonstrate the effect to fixing to the wrong threshold
        %% Hence, we determine the threshold corresponding to the received energy
        threshold_true = (sqrt(2/t_sen) * qfuncinv(P_d_d) + 1) * (signal_power + noise_power);
        %threshold_true = 1.0115e-10
        % Test: To investigated the effect of the increasing the threshold
        % on the distribution --> Not so relevent for the final
        % investigation
        % threshold_true = (1 - 0.03) * threshold_true;
        % Now capture the distortion on the Probability of dectection
        for i=1:M
                % Estimate the received power P_reg with additive noise amplitude (zero mean and noise variance) 
                samples = sqrt(signal_power) * ones(1, t_est) + random('norm', 0, sqrt(noise_power), 1, t_est);
                energy(i) = mean(samples.^2);                                  
        end 
        %energy = normrnd(signal_power + noise_power,...
                %sqrt(2/t_est) * (signal_power + noise_power), 1 , M);
        P_d = qfunc((threshold_true - energy) ./ (sqrt(2/t_sen) * energy));
  %      P_d = marcumq( sqrt(t_sen * (energy - noise_power)), sqrt(t_sen * threshold_true), t_sen/2);

        
        bins = 200;
        [f x_pts] = hist(P_d, bins);    
        
        %% Verify the analytical expression
        e = threshold_true;
        p = signal_power + noise_power;
        ts = sqrt(2/t_sen);
        te = sqrt(2/t_est);
        %pdf_pd_ana = e * exp( (erfcinv(2 * x_pts)).^2 - (e./(1 + ts * sqrt(2) * erfcinv(2 * x_pts)) - p).^2 /(sqrt(2) * p * te)^2) ...
        %* ts ./ (p * te * ( 1 + sqrt(2) * ts * erfcinv(2 * x_pts)).^2);
        
        %% Back to Q function
        scale = 1;
        pdf_pd_ana = 1/scale * e * exp( (qfuncinv(x_pts * scale)/sqrt(2)).^2 - (e./(1 + ts * qfuncinv(x_pts * scale)) - p).^2 /(sqrt(2) * p * te)^2)...
        * ts ./ (p * te * ( 1 + ts * qfuncinv(x_pts * scale)).^2);
    
        %% Applying MarcumQ function
        %pdf_pd_ana = 
        
        func_exp_pd = @(t) t .* e .* exp( (qfuncinv(t)/sqrt(2)).^2 - (e./(1 + ts * qfuncinv(t)) - p).^2 /(sqrt(2) * p * te)^2)...
        .* ts ./ (p * te * ( 1 + ts * qfuncinv(t)).^2);
    
        %Exp_p_d = integral(func_exp_pd,0,1); 
        
    
        % Plotting curves
        bar(x_pts,f/sum(f));    
        hold on,
        plot(x_pts, pdf_pd_ana * ((x_pts(3) - x_pts(2))), 'r', 'Linewidth',2.5);
        grid on,
    end
    
    if 0 % Long term scenario
        % Case where the energy measurement consisiting of K samples and where K samples sees N 
        % different realization of the channel and 1 < N < K
        M = 1e5;                    % Num of realizations
        K =  2000;                  % Num of samples exprnd(1, 1 ,1)  
        noise_power = 1e-10;
        x = zeros(1,M);
        SNR = -10;
        snr = 10^(SNR/10);
        signal_power = noise_power * 10^(SNR/10);    
        energy = zeros(1, M);
        m = 4;    
        g = gamrnd(m, 1, 1 ,M);
        for i=1:M       
            energy(i) = mean((sqrt(g(i) * signal_power) *...
                ones(1, K) + normrnd(0,sqrt(noise_power), 1 ,K)).^2);
        end
        bins = 500;
        [f x_pts] = hist(energy, bins);    
        %% Moments determined from the
        mean_nc = mean(energy);
        var_nc = moment(energy, 2);
        third_moment_nc_ana_wo_ch = moment(energy,3);
        fourth_moment_nc_ana_wo_ch = moment(energy,4);
        
        %% Moments determined from MGF
        mean_nc_ana = 1/K * K * (1 + m * snr) * noise_power;
        var_nc_ana = 1/K^2 * K * (2 + 4 * m * snr + K * m * snr^2 ) * noise_power^2; 
        third_central_nc_ana = 1/K^3 * 2 * K * (4 + m * snr * (12 + 6 * snr * K + snr^2 * K^2)) * noise_power^3; 
        fourth_moment_nc_ana = 1/K^4 * 6 * K * (8 + m * snr * (32 + 24 * snr * K + 8 * snr^2 * K^2 + snr^3 * K^3)) * noise_power^4;
        
        %% Gamma approximation
        %% Using the first two moments, we determine the parameters of Gamma distribution
        b = var_nc_ana/mean_nc_ana;
        a = mean_nc_ana/b;
        pdf_approx_gam  = gampdf(x_pts, a, b);  

        
        %% To analyze the peformformance of te Gamma approximation, we perform curve fitting  
        %para_gamma_fit = gamfit(energy);
        %pdf_curvefit  = gampdf(x_pts, para_gamma_fit(1), para_gamma_fit(2));

        % Log Normal approximation
        %para_logn_fit = lognfit(energy);
        %pdf_curvefit = lognpdf(x_pts, para_logn_fit(1), para_logn_fit(2));
        
        %% F-distribution
        %para_F_fit = mle(x,'pdf',@fpdf,'start',1:2)
        
        %% Generalized Extreme Value distribution
        para_gevd_fit = gevfit(energy);
        pdf_curvefit_gevd = gevpdf(x_pts, para_gevd_fit(1), para_gevd_fit(2), para_gevd_fit(3));

        if 1 %test
            k_evd_t = para_gevd_fit(1);
            g_1_t = gamma(1 - 1 * k_evd_t);
            g_2_t = gamma(1 - 2 * k_evd_t);
            g_3_t = gamma(1 - 3 * k_evd_t);
            skewnes_t =  (g_3_t  - 3 * g_1_t * g_2_t + 2 * g_1_t^3)/(g_2_t - g_1_t^2)^1.5;  
            var_t = para_gevd_fit(2)^2* (g_2_t  - g_1_t^2)/(k_evd_t)^2;
            mean_t = para_gevd_fit(3) + para_gevd_fit(2) *...
                (g_1_t - 1)/k_evd_t;
        end
        
        % Exponential approximation
        %para_exp_fit = expfit(energy);
        %pdf_curvefit = exppdf(x_pts, para_exp_fit);
        
        %% Burr Distribution
        para_burr = fitdist(energy','burr');
        pdf_curvefit_burr = pdf('burr',x_pts, para_burr.alpha, para_burr.c, para_burr.k);

        %% Inverse Gaussian distribution
        para_invgauss = fitdist(energy', 'InverseGaussian');
        pdf_curvefit_invgauss = pdf('InverseGaussian',x_pts,...
            para_invgauss.mu, para_invgauss.lambda);

        
        % Plotting curves
        bar(x_pts,f/sum(f));    
        hold on,
        plot(x_pts, pdf_approx_gam * ((x_pts(3) -x_pts(2))), 'r', 'Linewidth',2.5);
        grid on,
    end
  
end

%% 2) Analysis to determine the densities of C_0 and C_1 using Gamma Aprroximation
if 0
    %% Case 1: The case when the channel estimation is enabled by matched filtering only
    if 0       
        M = 1e5;                       % Number of realizations
        rho_s = 0.8;                   % Correlation coefficient for channel h_s
        rho_p2 = 0.8;                  % Correlation coefficient for channel h_s

        snr_s = 10^(-10/10);           % snr for the channel g_s  
        snr_p2 = 10^(-10/10);          % snr for the channel g_p2 

        g_s = zeros(1,M);              % power gain for the channel g_s 
        g_p2 = zeros(1,M);             % power gain for the channel g_p2  

        dof_s = 1;                     % degree of freedom for the channel g_s
        dof_p2 = 1;                    % degree of freedom for the channel g_s

        %%  Simulation
        %% g_s is supposed to be non-central chi squared distributed with 1 degree of freedom 
        for i=1:M
            g_s(i) = sum( normrnd(rho_s, sqrt(snr_s) * sqrt(1 - rho_s^2), 1 ,dof_s).^2);

            g_p2(i) = sum(normrnd(rho_p2, sqrt(snr_p2) * sqrt(1 - rho_p2^2), 1 ,dof_p2).^2);
        end

        bins = 200;

        %% density for SNR_s
        %% Gamma Aprroximation
        lambda_s =  dof_s *  rho_s^2/(snr_s * (1 - rho_s^2));
        mean_ana_s = snr_s * (1 - rho_s^2) * (dof_s + lambda_s);
        var_ana_s = snr_s^2 * (1 - rho_s^2)^2 * (2 * dof_s + 4 * lambda_s);
        mean_sim_s = mean(g_s);
        var_sim_s = var(g_s);

        b_gamma_s = var_ana_s/mean_ana_s;
        a_gamma_s = mean_ana_s/b_gamma_s;

        %  Plotting Curves    
        [f x_pts] = hist(g_s, bins);      
        figure(1);
        bar(x_pts,f/sum(f)); 
        pdf_approx_gamma_s  = 1/b_gamma_s^(a_gamma_s) *...
            (x_pts).^(a_gamma_s - 1)/gamma(a_gamma_s) .* exp(- x_pts/b_gamma_s);
        hold on,
        plot(x_pts, pdf_approx_gamma_s * ((x_pts(3) -x_pts(2))), 'r', 'Linewidth',2.5);


        %% density of Capacity C_0 = log2(1 + SNR_s)
        C_0 = log2(ones(1,M) +  g_s);
        [f x_pts] = hist(C_0, bins); 

        pdf_gamma_C_0 = 2.^(x_pts) * log(2) .* 1/b_gamma_s^(a_gamma_s)...
            .* (2.^x_pts - 1).^(a_gamma_s - 1)/gamma(a_gamma_s)...
            .* exp(- (2.^x_pts - 1)/b_gamma_s);

        %  Plotting Curves    
        figure(2);
        bar(x_pts,f/sum(f)); 
        hold on,
        plot(x_pts, pdf_gamma_C_0 * ((x_pts(3) -x_pts(2))), 'r', 'Linewidth',2.5);

        %% Determining the expectation
        % Simulated
        mean_C_0_sim = mean(C_0)

        % Numerically evalauted
        func_exp_C_0 = @(t) t .* 2.^(t) * log(2) .* 1/b_gamma_s^(a_gamma_s)...
            .* (2.^t - 1).^(a_gamma_s - 1)/gamma(a_gamma_s)...
            .* exp(- (2.^t - 1)/b_gamma_s);
        mean_C_O_ni = integral(func_exp_C_0,0, 100*max(x_pts))

        %% density for 1 + SNR_p
        g_p2 = (ones(1, M) + g_p2);

        lambda_p2 =  dof_p2 *  rho_p2^2/(snr_p2 * (1 - rho_p2^2));
        mean_ana_p2 = snr_p2 * (1 - rho_p2^2) * (dof_p2 + lambda_p2);
        var_ana_p2 = snr_p2^2 * (1 - rho_p2^2)^2 * (2 * dof_p2 + 4 * lambda_p2);
        mean_sim_p2 = mean(g_p2);
        var_sim_p2 = var(g_p2);    

        b_gamma_p2 = var_ana_p2/mean_ana_p2;
        a_gamma_p2 = mean_ana_p2/b_gamma_p2;          


        % Plotting curves
        [f x_pts] = hist(g_p2, bins);

        figure(3);
        bar(x_pts,f/sum(f));       
        pdf_approx_gamma_p2  = 1/b_gamma_p2^(a_gamma_p2) *...
            (x_pts - 1).^(a_gamma_p2 - 1)/gamma(a_gamma_p2) .* exp(- (x_pts - 1)/b_gamma_p2);
        hold on,
        plot(x_pts, pdf_approx_gamma_p2 * ((x_pts(3) -x_pts(2))), 'r', 'Linewidth',2.5);


        %% density of SNR_s/(1 + SNR_p)   
        [f x_pts] = hist(g_s./g_p2, bins);     

        figure(4);
        bar(x_pts,f/sum(f));

        %% Using Numerical Integration -- Analytical Solution is working
        if 0
            for i = 1:bins
                z = x_pts(i);
                if  1
                   f_density_SNR_s = @(y) 1/(gamma(a_gamma_s) * b_gamma_s^a_gamma_s)...
                       * 1/(gamma(a_gamma_p2) * b_gamma_p2^a_gamma_p2) .* y .* (z * y).^(a_gamma_s-1)...
                       .* exp(- z  * y/ b_gamma_s) .* (y - 1).^(a_gamma_p2 - 1) .*...
                       exp(-(y-1)/b_gamma_p2);
                   pdf_approx_SNR_s_ni(i) = integral(f_density_SNR_s, 1, Inf);            
                end
            end
        end    

        %% Analytical Expression
        pdf_approx_SNR_s_ana = 1/(b_gamma_p2^a_gamma_p2 * b_gamma_s^a_gamma_s) * exp(1/b_gamma_p2) .* x_pts.^(a_gamma_s -1) .*...
            ((1/b_gamma_p2 + x_pts/b_gamma_s).^(-a_gamma_s - a_gamma_p2) * gamma(a_gamma_p2 + a_gamma_s)...
            .* hypergeom(1 - a_gamma_p2, 1 - a_gamma_p2 - a_gamma_s, -(b_gamma_s + b_gamma_p2 * x_pts)/(b_gamma_s * b_gamma_p2))...
            + gamma(-a_gamma_p2 - a_gamma_s) * 1/gamma(-a_gamma_s) * gamma(a_gamma_p2) *...
            hypergeom(1 + a_gamma_s, 1 + a_gamma_p2 + a_gamma_s, -(b_gamma_s + b_gamma_p2 * x_pts)/(b_gamma_s * b_gamma_p2))...
            )/ (gamma(a_gamma_s) * gamma(a_gamma_p2)); 


        hold on,
        %plot(x_pts, pdf_approx_SNR_s_ni * ((x_pts(3) - x_pts(2))), 'r', 'Linewidth',2.5);
        %hold on,
        plot(x_pts, abs(pdf_approx_SNR_s_ana) * ((x_pts(3) - x_pts(2))), 'b', 'Linewidth',2.5);

        %% density of Capacity C_1 = log(1 + SNR_s/(SNR_p + 1))
        C_1 = log2(1  + g_s./g_p2);
        [f x_pts] = hist(C_1, bins);

        %% Using Numerical Integration
        if 0
            for i = 1:bins
                z = 2^x_pts(i) - 1;
                if  1
                   f_density_C_1 = @(y) (z + 1) * log(2) * 1/(gamma(a_gamma_s) * b_gamma_s^a_gamma_s)...
                       * 1/(gamma(a_gamma_p2) * b_gamma_p2^a_gamma_p2) .* y .* (z * y).^(a_gamma_s-1)...
                       .* exp(- z  * y/ b_gamma_s) .* (y - 1).^(a_gamma_p2 - 1) .*...
                       exp(-(y-1)/b_gamma_p2);
                   pdf_approx_C_1_sim(i) = integral(f_density_C_1, 1, Inf);            
                end
            end
        end    

        %% Analytical Expression
        pdf_approx_C_1_ana = (2.^x_pts) * log(2) * 1/(b_gamma_p2^a_gamma_p2 * b_gamma_s^a_gamma_s) * exp(1/b_gamma_p2) .* (2.^x_pts - 1).^(a_gamma_s -1) .*...
            ((1/b_gamma_p2 + (2.^x_pts - 1)/b_gamma_s).^(-a_gamma_s - a_gamma_p2) * gamma(a_gamma_p2 + a_gamma_s)...
            .* hypergeom(1 - a_gamma_p2, 1 - a_gamma_p2 - a_gamma_s, -(b_gamma_s + b_gamma_p2 * (2.^x_pts - 1))/(b_gamma_s * b_gamma_p2))...
            + gamma(-a_gamma_p2 - a_gamma_s) * 1/gamma(-a_gamma_s) * gamma(a_gamma_p2) *...
            hypergeom(1 + a_gamma_s, 1 + a_gamma_p2 + a_gamma_s, -(b_gamma_s + b_gamma_p2 * (2.^x_pts - 1))/(b_gamma_s * b_gamma_p2))...
            )/ (gamma(a_gamma_s) * gamma(a_gamma_p2));

        figure(5);
        bar(x_pts,f/sum(f));
        %hold on,
        %plot(x_pts, pdf_approx_C_1_sim * ((x_pts(3) - x_pts(2))), 'r', 'Linewidth',2.5);
        hold on,
        plot(x_pts, pdf_approx_C_1_ana * ((x_pts(3) - x_pts(2))), 'b', 'Linewidth',2.5);


        %% Expected capacity C_1
        mean_C_1_sim = mean(C_1)  


        f_exp_C_1_s = @(y,z) z .* (2.^z) .* log(2) * 1/(gamma(a_gamma_s) * b_gamma_s^a_gamma_s)...
                       * 1/(gamma(a_gamma_p2) * b_gamma_p2^a_gamma_p2) .* y .* ((2.^z - 1) .* y).^(a_gamma_s-1)...
                       .* exp(- (2.^z- 1)  .* y/ b_gamma_s) .* (y - 1).^(a_gamma_p2 - 1) .*...
                       exp(-(y-1)/b_gamma_p2);
        mean_C_1_ni = integral2(f_exp_C_1_s, 1, Inf, 0, 50*max(x_pts))

        f_exp_C_1_s = @(z) z .* (2.^z) * log(2) * 1/(b_gamma_p2^a_gamma_p2 * b_gamma_s^a_gamma_s) *...
            exp(1/b_gamma_p2) .* (2.^z - 1).^(a_gamma_s -1) .* ((1/b_gamma_p2 +...
            (2.^z - 1)/b_gamma_s).^(-a_gamma_s - a_gamma_p2) * gamma(a_gamma_p2 + a_gamma_s)...
            .* hypergeom(1 - a_gamma_p2, 1 - a_gamma_p2 - a_gamma_s,...
            -(b_gamma_s + b_gamma_p2 * (2.^z - 1))/(b_gamma_s * b_gamma_p2))...
            + gamma(-a_gamma_p2 - a_gamma_s) * 1/gamma(-a_gamma_s) * gamma(a_gamma_p2) *...
            hypergeom(1 + a_gamma_s, 1 + a_gamma_p2 + a_gamma_s,...
            -(b_gamma_s + b_gamma_p2 * (2.^z - 1))/(b_gamma_s * b_gamma_p2)))/...
            (gamma(a_gamma_s) * gamma(a_gamma_p2));
       %f_exp_C_1_s = @(z) hypergeom(1 - a_gamma_p2, 1 - a_gamma_p2 - a_gamma_s,...
       %     -(b_gamma_s + b_gamma_p2 * (2.^z - 1))/(b_gamma_s * b_gamma_p2));
       mean_C_1_ni = integral(f_exp_C_1_s, 0, 20)

       %f_exp_C_1_s = @(z) hypergeom(1 + a_gamma_s, 1 + a_gamma_p2 + a_gamma_s,...
       %      -(b_gamma_s + b_gamma_p2 * (2.^z - 1))/(b_gamma_s * b_gamma_p2));
       %mean_C_1_ni = integral(f_exp_C_1_s, 0, 25)

    end

    %% Case 2: The case when the channel estimation is enabled by matched filtering and energy detection (Model 1) 
    if 0
        M = 5e3;                       % Number of realizations
        K = 100;                       % Number of samples used for energy detection 
        rho_p2 = 0.8;                  % Correlation coefficient for channel h_s

        snr_s = 10^(-10/10);           % snr for the channel g_s  
        snr_p2 = 10^(-10/10);          % snr for the channel g_p2  
        
        g_s = zeros(1,M);              % Power gain for the channel g_s 
        snr_p_ = zeros(1,M);           % snr_p_prime = 1 + snr_p  
        
        noise_power = 10^(-100/10);    % Noise power at SR 
        dof_s = 1;                     % degree of freedom for the channel g_s

       %%  Simulation
       %% g_s is supposed to be non-central chi squared distributed with 1 degree of freedom 
        for i=1:M
            g_s(i) = sum( normrnd(sqrt(snr_s), 1 ,dof_s).^2);
            snr_p_(i) = mean((sqrt(snr_p2  * noise_power) * ones(1, K)...
                + random('norm', 0, sqrt(noise_power), 1, K)).^2);
        end

        bins = 200;

       %% density for SNR_s
       %% Gamma Aprroximation
        lambda_s =  dof_s *  rho_s^2/(snr_s * (1 - rho_s^2));
        mean_ana_s = snr_s * (1 - rho_s^2) * (dof_s + lambda_s);
        var_ana_s = snr_s^2 * (1 - rho_s^2)^2 * (2 * dof_s + 4 * lambda_s);
        mean_sim_s = mean(g_s);
        var_sim_s = var(g_s);

        b_gamma_s = var_ana_s/mean_ana_s;
        a_gamma_s = mean_ana_s/b_gamma_s;

        %  Plotting Curves    
        [f x_pts] = hist(g_s, bins);      
        figure(1);
        bar(x_pts,f/sum(f)); 
        pdf_approx_gamma_s  = 1/b_gamma_s^(a_gamma_s) *...
            (x_pts).^(a_gamma_s - 1)/gamma(a_gamma_s) .* exp(- x_pts/b_gamma_s);
        hold on,
        plot(x_pts, pdf_approx_gamma_s * ((x_pts(3) -x_pts(2))), 'r', 'Linewidth',2.5);


        %% density of Capacity C_0 = log2(1 + SNR_s)
        C_0 = log2(ones(1,M) +  g_s);
        [f x_pts] = hist(C_0, bins); 

        pdf_gamma_C_0 = 2.^(x_pts) * log(2) .* 1/b_gamma_s^(a_gamma_s)...
            .* (2.^x_pts - 1).^(a_gamma_s - 1)/gamma(a_gamma_s)...
            .* exp(- (2.^x_pts - 1)/b_gamma_s);

        % Plotting Curves    
        figure(2);
        bar(x_pts,f/sum(f)); 
        hold on,
        plot(x_pts, pdf_gamma_C_0 * ((x_pts(3) -x_pts(2))), 'r', 'Linewidth',2.5);

       %% Determining the expectation
        % Simulated
        mean_C_0_sim = mean(C_0)

       %% Numerically evalauted
        func_exp_C_0 = @(t) t .* 2.^(t) * log(2) .* 1/b_gamma_s^(a_gamma_s)...
            .* (2.^t - 1).^(a_gamma_s - 1)/gamma(a_gamma_s)...
            .* exp(- (2.^t - 1)/b_gamma_s);
        mean_C_O_ni = integral(func_exp_C_0,0, 100*max(x_pts))
        
       %% Analytical Expression
        if 0
            mean_C_O_ana = 1/( gamma(a_gamma_s) * log(2)) * (-(1/b_gamma_s)^(a_gamma_s) * pi *...
                csc(a_gamma_s * pi) * gamma(a_gamma_s, - 1/b_gamma_s) + 1/b_gamma_s *...
                (-1/b_gamma_s)^(a_gamma_s) * gamma(a_gamma_s - 1) *...
                hypergeom([1, 1], [2, 2 - a_gamma_s], 1/b_gamma_s) + gamma(a_gamma_s) *....
                ( (1/b_gamma_s)^a_gamma_s * pi * csc(a_gamma_s * pi) - (- 1/b_gamma_s)^(a_gamma_s) *...
                log(1/b_gamma_s) + (- 1/b_gamma_s)^(a_gamma_s) * psi(0, a_gamma_s)));  
        end

       %% density for SNR_p_
        mean_ana_p_ = noise_power * ( 1  + snr_p2);
        var_ana_p_ = noise_power^2/K * ( 2  + 4 * snr_p2);
        mean_sim_p_ = mean(snr_p_);
        var_sim_p_ = var(snr_p_);    

        b_gamma_p_ = var_ana_p_/mean_sim_p_;
        a_gamma_p_ = mean_sim_p_/b_gamma_p_;          


        % Plotting curves
        [f x_pts] = hist(snr_p_, bins);

        figure(3);
        bar(x_pts,f/sum(f));       
        pdf_approx_gamma_p_  = exp(a_gamma_p_ * log(1/b_gamma_p_) +...
            (a_gamma_p_ - 1) * log(x_pts) - gammaln(a_gamma_p_) - x_pts/b_gamma_p_);
        hold on,
        plot(x_pts, pdf_approx_gamma_p_ * ((x_pts(3) -x_pts(2))), 'r', 'Linewidth',2.5);


       %% density of SNR =  SNR_s/(1 + SNR_p)  =  SNR_s/SNR_p_ 
        [f x_pts] = hist(g_s./snr_p_, bins);     

        figure(4);
        bar(x_pts,f/sum(f));

       %% Analytical Expression of the density
        pdf_approx_SNR_ana = exp((a_gamma_s - 1) * log(x_pts) - (a_gamma_s + a_gamma_p_) *...
            log(1/b_gamma_p_ + x_pts/b_gamma_s) + gammaln(a_gamma_p_ + a_gamma_s) -...
            gammaln(a_gamma_p_) - gammaln(a_gamma_s) - a_gamma_p_ * log(b_gamma_p_) -... 
            a_gamma_s * log(b_gamma_s)); 


        hold on,
        plot(x_pts, abs(pdf_approx_SNR_ana) * ((x_pts(3) - x_pts(2))), 'b', 'Linewidth',2.5);

       %% density of Capacity C_1 = log(1 + SNR)
        C_1 = log2(1  + g_s./snr_p_);
        [f x_pts] = hist(C_1, bins);

       %% Analytical Expression
        pdf_approx_C_1_ana = log(2) * exp(log(2.^x_pts)  + (a_gamma_s - 1) * log(2.^x_pts - 1)...
            - (a_gamma_s + a_gamma_p_) * log(1/b_gamma_p_ + (2.^x_pts - 1)/b_gamma_s)...
            + gammaln(a_gamma_p_ + a_gamma_s) - gammaln(a_gamma_p_) - gammaln(a_gamma_s)...
            - a_gamma_p_ * log(b_gamma_p_) - a_gamma_s * log(b_gamma_s));

        figure(5);
        bar(x_pts,f/sum(f));
        %hold on,
        %plot(x_pts, pdf_approx_C_1_sim * ((x_pts(3) - x_pts(2))), 'r', 'Linewidth',2.5);
        hold on,
        plot(x_pts, pdf_approx_C_1_ana * ((x_pts(3) - x_pts(2))), 'b', 'Linewidth',2.5);


       %% Expected capacity C_1
        mean_C_1_sim = mean(C_1)  
        
       %% Numerically evaluated
        f_exp_C_1_s = @(z) z .* log(2) .* exp(log(2.^z)  + (a_gamma_s - 1) * log(2.^z - 1)...
            - (a_gamma_s + a_gamma_p_) * log(1/b_gamma_p_ + (2.^z - 1)/b_gamma_s)...
            + gammaln(a_gamma_p_ + a_gamma_s) - gammaln(a_gamma_p_) - gammaln(a_gamma_s)...
            - a_gamma_p_ * log(b_gamma_p_) - a_gamma_s * log(b_gamma_s));
        mean_C_1_ni = integral(f_exp_C_1_s, 0, 50*max(x_pts))

        %% Analytical Expresssion of the expected capacity C_1
        
        
    end

    %% Case 3: The case when the channel estimation is enabled by matched filtering and energy detection (Model 2)    
    if 1
        M = 5e4;                             % Number of realizations
        K = 100;                             % Number of samples used for energy detection 
        rho_p2 = 0.8;                        % Correlation coefficient for channel h_s

        snr_s = 10^(-10/10);                 % received snr for the channel g_s  
        snr_p2 = 10^(-10/10);                % snr for the channel g_p2  
        noise_power = 10^(-100/10);          % Noise power at SR 
        g_s = snr_s * noise_power;           % Power gain for the channel g_s         
        tran_snr_ST = snr_s / g_s;           % Transmitted SNR = Transmitted power at PT/Noise at SR  = received SNR / g_s
        E_s = 1;                             % Pilot Energy 
        N_s = 10;                            % Number of Pilot Symobols 
        sigma_est_error = noise_power/N_s;   % Variance of the Estimation error for h_s 
        dof_s = 1;
        
        g_s_hat = zeros(1,M);                % Estimated power gain for the channel g_s 
        snr_p_ = zeros(1,M);                 % snr_p_prime = 1 + snr_p  
     
        
       %%  Simulation
       %% g_hat_s is supposed to be non-central chi squared distributed with 1 degree of freedom 
        for i=1:M
            g_s_hat(i) = sum(normrnd(sqrt(g_s), sqrt(sigma_est_error), 1, dof_s).^2);
            snr_p_(i) = mean((sqrt(snr_p2  * noise_power) * ones(1, K)...
                + random('norm', 0, sqrt(noise_power), 1, K)).^2);
        end

        bins = 200;

       %% density for SNR_s
       %% Gamma Aprroximation
        if 0    % Test         

            lambda_s = dof_s *  g_s / sigma_est_error;
            mean_ana_s = sigma_est_error * (dof_s + lambda_s);
            var_ana_s =  (sigma_est_error)^2 * (2 * dof_s + 4 * lambda_s);
            mean_sim_s = mean(g_s_hat);
            var_sim_s = var(g_s_hat);

            b_gamma_s = var_ana_s/mean_ana_s;
            a_gamma_s = mean_ana_s/b_gamma_s;

            %  Plotting Curves    
            [f x_pts] = hist(g_s_hat, bins);      
            figure(1);
            bar(x_pts,f/sum(f)); 
            pdf_approx_gamma_s = 1/b_gamma_s^(a_gamma_s) *...
                (x_pts).^(a_gamma_s - 1)/gamma(a_gamma_s) .* exp(- x_pts/b_gamma_s);
            hold on,
            plot(x_pts, pdf_approx_gamma_s * ((x_pts(3) -x_pts(2))), 'r', 'Linewidth',2.5);
        end

        
        %% Density of the SNR_s
        lambda_s = dof_s *  g_s / sigma_est_error;
        snr_s_hat = tran_snr_ST * (g_s_hat);
        mean_ana_s = sigma_est_error * tran_snr_ST * (dof_s + lambda_s);
        var_ana_s = (sigma_est_error)^2 * tran_snr_ST^2 * (2 * dof_s + 4 * lambda_s);
        mean_sim_s = mean(snr_s_hat);
        var_sim_s = var(snr_s_hat);        
                
        b_gamma_s = var_ana_s/mean_ana_s;
        a_gamma_s = mean_ana_s/b_gamma_s;
        
        [f x_pts] = hist(snr_s_hat, bins);      
        figure(1);
        bar(x_pts,f/sum(f)); 
        pdf_approx_gamma_s = 1/b_gamma_s^(a_gamma_s) *...
            (x_pts).^(a_gamma_s - 1)/gamma(a_gamma_s) .* exp(- x_pts/b_gamma_s);
        hold on,
        plot(x_pts, pdf_approx_gamma_s * ((x_pts(3) -x_pts(2))), 'r', 'Linewidth',2.5);
        

        %% density of Capacity C_0 = log2(1 + SNR_s)
        C_0 = log2(ones(1,M) +  snr_s_hat);
        [f x_pts] = hist(C_0, bins); 

        pdf_gamma_C_0 = 2.^(x_pts) * log(2) .* 1/b_gamma_s^(a_gamma_s)...
            .* (2.^x_pts - 1).^(a_gamma_s - 1)/gamma(a_gamma_s)...
            .* exp(- (2.^x_pts - 1)/b_gamma_s);

        % Plotting Curves    
        figure(2);
        bar(x_pts,f/sum(f)); 
        hold on,
        plot(x_pts, pdf_gamma_C_0 * ((x_pts(3) -x_pts(2))), 'r', 'Linewidth',2.5);

       %% Determining the expectation
        % Simulated
        mean_C_0_sim = mean(C_0);

       %% Numerically evalauted
        func_exp_C_0 = @(t) t .* 2.^(t) * log(2) .* 1/b_gamma_s^(a_gamma_s)...
            .* (2.^t - 1).^(a_gamma_s - 1)/gamma(a_gamma_s)...
            .* exp(- (2.^t - 1)/b_gamma_s);
        mean_C_O_ni = integral(func_exp_C_0,0, max(x_pts))
        
       %% Analytical Expression
        if 0
            mean_C_O_ana = 1/( gamma(a_gamma_s) * log(2)) * (-(1/b_gamma_s)^(a_gamma_s) * pi *...
                csc(a_gamma_s * pi) * gamma(a_gamma_s, - 1/b_gamma_s) + 1/b_gamma_s *...
                (-1/b_gamma_s)^(a_gamma_s) * gamma(a_gamma_s - 1) *...
                hypergeom([1, 1], [2, 2 - a_gamma_s], 1/b_gamma_s) + gamma(a_gamma_s) *....
                ( (1/b_gamma_s)^a_gamma_s * pi * csc(a_gamma_s * pi) - (- 1/b_gamma_s)^(a_gamma_s) *...
                log(1/b_gamma_s) + (- 1/b_gamma_s)^(a_gamma_s) * psi(0, a_gamma_s)));  
        end

       %% density for SNR_p_
        mean_ana_p_ = noise_power * ( 1  + snr_p2);
        var_ana_p_ = noise_power^2/K * ( 2  + 4 * snr_p2);
        mean_sim_p_ = mean(snr_p_);
        var_sim_p_ = var(snr_p_);    

        b_gamma_p_ = var_ana_p_/mean_sim_p_;
        a_gamma_p_ = mean_sim_p_/b_gamma_p_;          


        % Plotting curves
        [f x_pts] = hist(snr_p_, bins);

        figure(3);
        bar(x_pts,f/sum(f));       
        pdf_approx_gamma_p_  = exp(a_gamma_p_ * log(1/b_gamma_p_) +...
            (a_gamma_p_ - 1) * log(x_pts) - gammaln(a_gamma_p_) - x_pts/b_gamma_p_);
        hold on,
        plot(x_pts, pdf_approx_gamma_p_ * ((x_pts(3) -x_pts(2))), 'r', 'Linewidth',2.5);


       %% density of SNR =  SNR_s/(1 + SNR_p)  =  SNR_s/SNR_p_ 
        [f x_pts] = hist(snr_s_hat./snr_p_, bins);     

        figure(4);
        bar(x_pts,f/sum(f));

       %% Analytical Expression of the density
        pdf_approx_SNR_ana = exp((a_gamma_s - 1) * log(x_pts) - (a_gamma_s + a_gamma_p_) *...
            log(1/b_gamma_p_ + x_pts/b_gamma_s) + gammaln(a_gamma_p_ + a_gamma_s) -...
            gammaln(a_gamma_p_) - gammaln(a_gamma_s) - a_gamma_p_ * log(b_gamma_p_) -... 
            a_gamma_s * log(b_gamma_s)); 


        hold on,
        plot(x_pts, abs(pdf_approx_SNR_ana) * ((x_pts(3) - x_pts(2))), 'b', 'Linewidth',2.5);

       %% density of Capacity C_1 = log(1 + SNR)
        C_1 = log2(1  + snr_s_hat./snr_p_);
        [f x_pts] = hist(C_1, bins);

       %% Analytical Expression
        pdf_approx_C_1_ana = log(2) * exp(log(2.^x_pts)  + (a_gamma_s - 1) * log(2.^x_pts - 1)...
            - (a_gamma_s + a_gamma_p_) * log(1/b_gamma_p_ + (2.^x_pts - 1)/b_gamma_s)...
            + gammaln(a_gamma_p_ + a_gamma_s) - gammaln(a_gamma_p_) - gammaln(a_gamma_s)...
            - a_gamma_p_ * log(b_gamma_p_) - a_gamma_s * log(b_gamma_s));

        figure(5);
        bar(x_pts,f/sum(f));
        %hold on,
        %plot(x_pts, pdf_approx_C_1_sim * ((x_pts(3) - x_pts(2))), 'r', 'Linewidth',2.5);
        hold on,
        plot(x_pts, pdf_approx_C_1_ana * ((x_pts(3) - x_pts(2))), 'b', 'Linewidth',2.5);


       %% Expected capacity C_1
        mean_C_1_sim = mean(C_1)  
        
       %% Numerically evaluated
        f_exp_C_1_s = @(z) z .* log(2) .* exp(log(2.^z)  + (a_gamma_s - 1) * log(2.^z - 1)...
            - (a_gamma_s + a_gamma_p_) * log(1/b_gamma_p_ + (2.^z - 1)/b_gamma_s)...
            + gammaln(a_gamma_p_ + a_gamma_s) - gammaln(a_gamma_p_) - gammaln(a_gamma_s)...
            - a_gamma_p_ * log(b_gamma_p_) - a_gamma_s * log(b_gamma_s));
        mean_C_1_ni = integral(f_exp_C_1_s, 0, 50)

        %% Analytical Expresssion of the expected capacity C_1
    end
end

%% 3) Validate the expection on the probability of detection, where 
%% transmitted source follows a Gaussian distribution
if 1   
    M = 1e4;                                        % Number of realizations
    test = 100;                                     % Estimation time
    tsen = 100;                                     % Sensing time  
    noise_power = 10^(-100/10);                     % Noise power        
    rcvd_energy = zeros(1, M);                      % Received Energy     
    pdc = 0.90;                                     % Constraint on probability of detection  
    snr_rcvd = 10.^([-00]/10);                      % snr for the channel g_p1  

    rcvd_energy_bar = (snr_rcvd + 1) * noise_power; % true value of received energy
    bins = 200;
      
    %% Simulation
    %% Estimation of rcvd_energy is supposed to be non-central chi squared distributed 
    for i=1:M
        rcvd_energy(i) = mean(random('norm',...
            random('norm',0, sqrt(snr_rcvd * noise_power),1, test),...
            sqrt(noise_power), 1, test).^2);
    end

    % Determine the threshold 
    threshold = gammaincinv(pdc, tsen/2,'upper') * 2 * rcvd_energy_bar / tsen;

    %% Realization of the random variable
    pd = gammainc(tsen/2, tsen/2 * threshold ./ rcvd_energy);  
    
    %% Expected value of pd
    % Simulated 
    exp_pd_sim = mean(pd);
    
    if 0
        t = linspace(min(rcvd_energy), max(rcvd_energy), 10000);
        func_pd = gammainc(tsen/2, tsen/2 * threshold ./ t) *...
            test/rcvd_energy_bar .* exp(-test/2 * log(2) - gammaln(test/2)  +...
            (test/2 - 1) * log(test * t/rcvd_energy_bar) + (-t * test/(2 * rcvd_energy_bar)));
        plot(t, func_pd); 
    end
    % Analytical -- Numerical Integration
    func_exp_pd = @(t) gammainc(tsen/2, tsen/2 * threshold ./ t) *...
            test/rcvd_energy_bar .* exp(-test/2 * log(2) - gammaln(test/2)  +...
            (test/2 - 1) * log(test * t/rcvd_energy_bar) + (-t * test/(2 * rcvd_energy_bar)));
    exp_pd_ana = integral(func_exp_pd, 0, threshold * 100);
end