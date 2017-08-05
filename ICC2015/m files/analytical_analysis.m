%% Investigating the distribution of the regulated power and power recieved at PR
%% Simulation of the Inv non central chi-squared distributed random variable
%% Analytical validatation of the density of the regulated power at ST and Power received at PR
clc,
close all,
clear all,
if 1
    %% AWGN or path loss case
    if 1        
        M = 1e4;                                          % Num of realizations
        K = 50;                
        pl = 10^(-100/10);                                % Path loss at link PR-ST
        pls = 10^(-100/10);                               % Path loss at link ST-SR
        np = 10^(-100/10);                                % Noise power at PR
        nps = 10^(-100/10);                               % Noise power at SR 
        tp = 10^(10/10);                                  % Transmit power at beacon
        lambda =  pl * K * tp / np;                       % non-cental chi-2 parameter
        it = 10^(-110/10);                                % Interference temperature
        nc = 1;                                           % normalization constant to compensate for
                                                          % scaling and shifting of the received power  
        %% Simulation 
        for i=1:M       
            ncx2_sim(i) = mean((sqrt(pl * ones(1, K) * tp) + normrnd(0, sqrt(np), 1 ,K)).^2);
        end       
        
        % Determine the scaling factor nc
        inv_ncx2_sim = it ./ ncx2_sim;                    
        scaled_inv_ncx2_sim = pl * inv_ncx2_sim;                
        nc = it/mean(scaled_inv_ncx2_sim);
        
        inv_ncx2_sim = it * nc./ ncx2_sim;               % Reguated power at ST (compensation included)  
        scaled_inv_ncx2_sim = pl .* inv_ncx2_sim;        % Power Received at PR 
        scaled2_inv_ncx2_sim = pls/nps * inv_ncx2_sim;   % SNR Recieved at SR
        cap_sim = log2(1 + scaled2_inv_ncx2_sim);        % Capacity at ST 
        
        
        %% Analytical validations        
        bins = 200;    
        
        %% Power Received at ST
        [f x_pts] = hist(ncx2_sim, bins);
        
        %% Final step: working
        besselln =  zeros(1, bins);        
        for i = 1:bins
            x_sym = sym(x_pts(i));
            besselln(i) = vpa(log(besseli(K/2 - 1, sqrt(pl * tp * K^2 * x_sym/np^2))));
        end       
        pdf_ncx2 = (x_pts(3) -x_pts(2))* K/np * 1/2 .* exp(- K/(2 * np) * (x_pts + (pl * tp))...
            + (K/4 - 1/2) * log((x_pts)./(pl * tp)) +besselln); 
        
        % Plot the curve
        figure(1);
        bar(x_pts,f/sum(f));
        hold on,
        plot(x_pts, pdf_ncx2, 'g', 'Linewidth',2.5);  
        title('Received Power at ST (P_{rcvd})');
    
        %% Power regulated at ST
        [f x_pts] = hist(inv_ncx2_sim, bins); 
        
        %% Final step: working
        besselln =  zeros(1, bins);  
        x_sym = sym(x_pts(i));
        for i = 1:bins
            x_sym = sym(x_pts(i));
            besselln(i) = vpa(log(besseli(K/2 - 1, sqrt(K^2 * pl * tp * it * nc./(x_sym * np^2)))));
        end
        pdf_inv_ncx2 = (x_pts(3) -x_pts(2)) * K/np * nc * it * (1./(x_pts)).^2 * 1/2 .*...
            exp(-K/(2 * np) * (it * nc./x_pts + pl * tp) +...
            (K/4 - 1/2) * log(it * nc./(x_pts * pl * tp)) + besselln);
        
        % Plot the curve
        figure(2);
        bar(x_pts,f/sum(f));
        hold on,
        plot(x_pts, pdf_inv_ncx2, 'g', 'Linewidth',2.5);
        title('Controlled Power at PR (P_c)');

        %% Power Received at PR
        [f x_pts] = hist(scaled_inv_ncx2_sim, bins); 
        
        %% Intermediate step: working
        %% without interfernce temp
        if 0
            besselln =  zeros(1, bins);
            for i = 1:bins
                x_sym = sym(x_pts(i));
                besselln(i) = vpa(log(besseli(K/2 - 1, sqrt(lambda * pl/(x_sym * np)))));
            end        
            pdf_scaled_inv_ncx2 = (x_pts(3) -x_pts(2)) * 1/np * 1/pl  *(pl./(x_pts)).^2 * 1/2 .*...
                exp(- (pl./(x_pts * np) + lambda)/2 + (K/4 - 1/2)* log(pl./(lambda * (x_pts * np)))...
                + besselln);
        end
        
        %% Final step: working
        %% including iterference temperature
        if 1
            if 0
                besselln =  zeros(1, bins);
                for i = 1:bins
                    x_sym = sym(x_pts(i));
                    besselln(i) = vpa(log(besseli(K/2 - 1, pl * K /np * sqrt(tp * it * nc/x_sym))));
                end
            end
            pdf_scaled_inv_ncx2 = (x_pts(3) -x_pts(2)) * (pl * it * K * nc)/np * (1./(x_pts)).^2 * 1/2 .*...
                exp(-K * pl/(2 * np) * (it * nc./x_pts + tp) +...
                (K/4 - 1/2)* log(it * nc./(tp * x_pts)) + besselln);
        end
        
        % Plot the curve
        figure(3);
        bar(x_pts,f/sum(f));
        hold on,
        plot(x_pts, pdf_scaled_inv_ncx2, 'g', 'Linewidth',2.5);
        title('Power Received at PR (P_p)');
        
        %% Received SNR at SR
        [f x_pts] = hist(scaled2_inv_ncx2_sim, bins); 
        
        %% Final step: working
        besselln =  zeros(1, bins);  
        for i = 1:bins
            x_sym = sym(x_pts(i));
            besselln(i) = vpa(log(besseli(K/2 - 1, sqrt(K^2 * pl * tp * it * nc * pls...
                ./(x_sym * np^2 * nps)))));
        end
        pdf_scaled2_inv_ncx2 = (x_pts(3) -x_pts(2)) * K * pls/(np * nps) * nc * it * (1./(x_pts)).^2 ...
        * 1/2 .* exp(-K/(2 * np) * (it * nc * pls ./(x_pts * nps)  + pl * tp) +...
            (K/4 - 1/2) * log(it * nc * pls./(x_pts * pl * tp * nps)) + besselln);
        
        % Plot the curve
        figure(4);
        bar(x_pts,f/sum(f));
        hold on,
        plot(x_pts, pdf_scaled2_inv_ncx2, 'g', 'Linewidth',2.5);
        title('Recieved SNR at SR (P_s)');
        
        %% Capacity at SR
        [f x_pts] = hist(cap_sim, bins);         
       
        %% Final step: working
        if 1
            besselln =  zeros(1, bins);
            for i = 1:bins
                x_sym = sym(x_pts(i));
                besselln(i) = vpa(log(besseli(K/2 - 1, sqrt(K^2 * pl * it * nc * tp * pls...
                    /((2^x_sym - 1) * nps * np^2)))));
            end 
            pdf_Cap_SR = (x_pts(3) -x_pts(2)) * (K * pls * nc * it)/ (2 * nps * np) *...
                2.^x_pts .* (1./(2.^x_pts - 1)).^2 .* log(2) .*... 
                exp(-K/(2 * np) * (it * nc * pls ./ ((2.^x_pts - 1) * nps) + pl * tp) +...
                (K/4 - 1/2) * log(it * nc * pls./(tp * pl * nps * (2.^x_pts - 1))) + besselln);
        end
        
        % Plot the curve
        figure(5);
        bar(x_pts,f/sum(f));
        hold on,
        plot(x_pts, pdf_Cap_SR, 'g', 'Linewidth',2.5);  
        title('Capacity at SR (R_s)');

    end

    if 0
        %% Fading channel case
        M = 1e4;                                          % Num of realizations
        K = 100;                                           % Number of samples used for energy detection
        tp = 10^(-15/10);                                 % Transmit power at beacon
        np = 10^(-100/10);                                % noise power 
        pl = 10^(-100/10);                                % distance dependent Path loss (PR - ST)
        it = 10^(-110/10);                                % Interference temperature
        nc = 1;                                           % normalization constant to compensate for
                                                          % scaling and shifting of the received power  
        snr = 1;
        lambda =  snr * K * pl / np;    
        bins = 100;    

        % Realizing channel
        g = snr * ones(1,M) .* gamrnd(1, 1, 1 ,M);
        for i=1:M       
            ncx2_sim(i) = mean((sqrt(g(i) * pl * tp) * ones(1, K) + normrnd(0, sqrt(np), 1 ,K)).^2);
        end       
        inv_ncx2_sim = it * nc./ ncx2_sim;
        scaled_inv_ncx2_sim = g * pl .* inv_ncx2_sim;
        
        %% To determine nc, calculate the distributions second time 
        nc = it/mean(scaled_inv_ncx2_sim);
        
        inv_ncx2_sim = it * nc./ ncx2_sim;
        scaled_inv_ncx2_sim = g * pl .* inv_ncx2_sim;

        [f x_pts] = hist(ncx2_sim, bins);  
        pdf_ncx2 = (x_pts(3) -x_pts(2)) * 1/2 .* 1/np * exp(- ((x_pts/np) + lambda)/2) .* ((x_pts/np)./lambda ).^(K/4 - 1/2)...
        .* besseli(K/2 - 1, sqrt(lambda * (x_pts/np))); 
    
        %% Regulated power at ST 
        [f x_pts] = hist(inv_ncx2_sim, bins);            

        %% Analytical Expresssion
        %% Bessel function case
        if 1
           %% version 1, where the snr is exponentially distributed
           %% intermediate Step: working 
            if 0
                pdf_inv_ncx2 = (x_pts(3) -x_pts(2)) * (K - 2) * snr * exp( -1./(2 * x_pts + snr * K * x_pts)...
                    + K/2 * log(2/snr + K) + (K - 2)/4 * log(1./(K * x_pts)) - (K - 2)/4 *...
                    log(K./x_pts) - log(2 * (2 + snr * K).^2 * x_pts.^2) - gammaln(K/2) + gammaln(K/2 - 1)) .*...
                    (1 - gammainc(K/2 - 1, snr * K./(4 * x_pts + 2 * snr * K * x_pts)));  
            end
           %% versrion 2, where the signal energy is exponentially distributed
           %% intermediate Step: working  
            if 0
                temp = pl * K./(4 * x_pts + 2 * pl * K * x_pts);
                pdf_inv_ncx2 = (x_pts(3) -x_pts(2)) * (K - 2) * exp( -1./(2 * x_pts + pl * K * x_pts)...
                    + (K/2  -2 ) * log(2 + K * pl) + (K - 2)/2 * log(1./(K * pl))... 
                    - log(2 * x_pts.^2) - gammaln(K/2) + gammaln(K/2 - 1)) .*...
                    (1 - gammainc(K/2 - 1, pl * K./(4 * x_pts + 2 * pl * K * x_pts))); 
            end
            %% versrion 2.5, where the signal energy is exponentially distributed and noise is included
            %% intermediate Step: working 
            if 0
                pdf_inv_ncx2 = zeros(1, bins);
                nterms = 500;
                for i = 1:bins
                    x_sym = sym(x_pts(i));
                    disp(strcat('pdf_inv_ncx2_25 bin: ',num2str(i)));
                    a_25 = -1./(2 * x_sym + K * pl * x_sym) +...
                        (K/2  -2 ) * log(2 + K * pl) +...
                        (K - 2)/2 * log(1./(K * pl))... 
                        - log(2 * x_sym^2) - gammaln(K/2) + gammaln(K/2 - 1);
                    temp_25 = pl * K./ (4 * x_sym + 2 * pl * K * x_sym);
                    b_25 = (log(1 - mygammaincln(K/2 - 1, temp_25, nterms)));
                    pdf_inv_ncx2(i) = (x_pts(3) -x_pts(2)) * (K - 2)...
                        * exp(double(a_25 + b_25));
                end                    
            end
           %% versrion 3, where the signal energy is exponentially distributed and noise is included
           %% intermediate Step: working
            if 0
                pdf_inv_ncx2 = zeros(1, bins);
                nterms = 500;
                for i = 1:bins
                    x_sym = sym(x_pts(i));
                     disp(strcat('pdf_inv_ncx2_30 bin: ',num2str(i)));
                    a_30 = -1./(2 * np * x_sym + K * pl * x_sym) +...
                        K/2 * log(2 + K * pl/np) + (K - 2)/4 * log(1./(K * pl * x_sym))... 
                        - (K/4 - 1/2) * log(K * pl ./ (np^2 * x_sym)) - ...
                        log(2 * (2 * np  + K * pl)^2 * x_sym^2) - gammaln(K/2) + gammaln(K/2 - 1);
                    temp_30 = pl * K./ (4 * np^2 * x_sym + 2 * np * pl * K * x_sym);
                    b_30 = log(1 - mygammaincln(K/2 - 1, temp_30, nterms));
                    pdf_inv_ncx2(i) = (x_pts(3) -x_pts(2)) * (K - 2) * np...
                        * exp(a_30 + b_30);
                end                    
            end
           %% versrion 4, where the signal energy (includes path loss, 
           %% noise power, iterference temperature, transmit power) is exponentially distributed
           %% Final Step: working
            if 0
                pdf_inv_ncx2 = zeros(1, bins);
                nterms = 200;
                for i = 1:bins
                    x_sym = sym(x_pts(i));
                     disp(strcat('pdf_inv_ncx2_40 bin: ',num2str(i)));
                    a_40 = -it * K * nc./(2 * np * x_sym + K * pl * tp * x_sym) +...
                        K/2 * log(2 + K * pl * tp/np) +...
                        (K + 2)/4 * log(it * nc./(pl * tp * x_sym))... 
                        - (K/4 - 1/2) * log(it  * K^2 * nc * pl * tp ./ (np^2 * x_sym)) - ...
                        log(2 * (2 * np  + K * pl * tp)^2 * x_sym) -...
                        gammaln(K/2) + gammaln(K/2 - 1);
                    temp_40 = it * K^2 * nc * pl * tp./...
                        (4 * np^2 * x_sym + 2 * K * np * pl * tp * x_sym);
                    b_40 = log(1 - mygammaincln(K/2 - 1, temp_40, nterms));
                    pdf_inv_ncx2(i) = (x_pts(3) -x_pts(2)) * (K - 2) * K * np * pl * tp...
                        * exp(a_40 + b_40);
                end                    
            end
        end        
        if 0
           %% Hypgeormetric case
           %% version 1, where the snr is exponentially distributed
           %% intermediate Step: working 
            pdf_inv_ncx2 = (x_pts(3) - x_pts(2)) * (K - 2) * K * snr * ...
                exp( -1./(2 * x_pts + snr * K * x_pts) + (K/2 + 2) * log(1./x_pts)...
                - K/2 * log(snr * K ./(2 * x_pts + snr * K * x_pts))... 
                - log(2 * (2 + snr * K).^2 ) - gammaln(K/2) + gammaln(K/2 - 1)) .*...
                (1 - gammainc(K/2 - 1, snr * K./(4 * x_pts + 2 * snr * K * x_pts)));
        end

        % Plotting curves
        figure(1);
        bar(x_pts,f/sum(f));
        %hold on,
        %plot(x_pts, pdf_inv_ncx2, 'g', 'Linewidth',2.5); 


        %% Power Received at PR
        [f x_pts] = hist(scaled_inv_ncx2_sim, bins);    
     
       
        %% Analytical Expresssion
        %% version 1, where the snr is exponentially distributed
        %% intermediate Step: working  
        if 0
            inter1 = (snr + 2 * x_pts + snr * K * x_pts).^2/snr^2;
            inter2 = sqrt((- 4 * K * x_pts + inter1)./x_pts.^2);
            pdf_scaled_inv_ncx2 = (x_pts(3) - x_pts(2)) * (K * exp((K/2 - 1)* log(2)... 
                + 2 * log(snr) + ((K - 2)/2) * log(1./(x_pts)) - K/2 * log(1./(snr * x_pts) .*...
                (snr + 2 * x_pts + snr * K * x_pts + snr * x_pts .* inter2))) .*...
                (8 - 4 * K  + inter1./x_pts + sqrt(inter1) .* inter2) ./...
                (snr * sqrt(inter1) .* sqrt(1 - 4 * K * x_pts./inter1) .*...
                (- 4 * snr^2 * K * x_pts + inter1 * snr^2))); 
        end
       %% versrion 2, where the signal energy (includes path loss) is exponentially distributed
       %% intermediate Step: working  
        if 0
            b_20 = zeros(1, bins);
            t_20 = zeros(1, bins);
            for i = 1:bins
                x_sym = sym(x_pts(i));
                t_20(i) = 4 * K * pl^2 * x_sym./(pl + 2 * x_sym + K* pl * x_sym).^2;
                b_20(i) = vpa(log(abs(hypergeom([(2 + K)/4, (4 + K) /4], K/2, t_20(i)))));
            end
            a_20 = (K - 2)/4 * log(1./(K * x_pts))  + (K + 2)/4 * log( K * pl^2 * x_pts ./ ...
                (pl + 2 * x_pts + K * pl * x_pts).^2);
           
            pdf_scaled_inv_ncx2_20 = (x_pts(3) - x_pts(2)) * 1./(pl * x_pts) .* exp(a_20 + b_20); 
        end
       %% versrion 2.5, where the signal energy (includes path loss) is exponentially distributed
       %% intermediate Step: working 
        if 0
            b_25 = zeros(1, bins);
            t_25 = zeros(1, bins);
            pl_np = pl/np; 
            for i = 1:bins
                x_sym = sym(x_pts(i));
                t_25(i) = 4 * K * pl_np^2 * x_sym./(pl_np + 2 * x_sym + K* pl_np * x_sym).^2;
                b_25(i) = vpa(log(abs(hypergeom([(2 + K)/4, (4 + K) /4], K/2, t_25(i)))));
            end
            a_25 = (K - 2)/4 * log(1./(K * x_pts))  + (K + 2)/4 * log( K * pl_np^2 * x_pts ./ ...
                (pl_np + 2 * x_pts + K * pl_np * x_pts).^2);
           
            pdf_scaled_inv_ncx2_25 = (x_pts(3) - x_pts(2)) * 1./(pl_np * x_pts) .* exp(a_25 + b_25); 
        end 
        
        %% versrion 3, where the signal energy (includes path loss and noise power) is exponentially distributed
        %% intermediate Step: working
        if 0
            b_30 = zeros(1, bins);
            t_30 = zeros(1, bins);
            for i = 1:bins
                x_sym = sym(x_pts(i));
                t_30(i) = 4 * K * pl^2 * x_sym/(pl + 2 * np * x_sym + K * pl * x_sym)^2;
                b_30(i) = vpa(log(abs(hypergeom([(2 + K)/4, (4 + K) /4], K/2, t_30(i)))));
            end
            a_30 = (K - 2)/4 * log(1./(K * x_pts))  + (K + 2)/4 * log( K * pl^2 * x_pts...
                ./(pl + 2 * np * x_pts + K * pl * x_pts).^2);
           
            pdf_scaled_inv_ncx2_30 = (x_pts(3) - x_pts(2)) * np./(pl * x_pts) .* exp(a_30 + b_30); 
        end  
        
        %% versrion 4, where the signal energy (includes path loss, 
        %% noise power, iterference temperature, transmit power) is exponentially distributed
        %% Final Step: working
        if 1
            b_40 = zeros(1, bins);
            t_40 = zeros(1, bins);
            for i = 1:bins
                x_sym = sym(x_pts(i));
                t_40(i) = 4 * it * K^2 * nc * pl^2 * tp * x_sym/...
                    (it * K * nc * pl + 2 * np * x_sym + K * pl * tp * x_sym)^2;
                b_40(i) = vpa(log(abs(hypergeom([(2 + K)/4, (4 + K) /4], K/2, t_40(i)))));
            end
            a_40 = (K + 2)/4 * log(it * nc./(tp * x_pts))  + (K + 2)/4 *...
                log( it * K^2 * nc * pl^2 * tp * x_pts...
                ./(it * K * nc * pl + 2 * np * x_pts + K * pl * tp * x_pts).^2);
           
            pdf_scaled_inv_ncx2_40 = (x_pts(3) - x_pts(2)) * np./(it * nc * pl)...
                .* exp(a_40 + b_40); 
        end
        
        % Plotting curves
        figure(2);
        bar(x_pts,f/sum(f));

        hold on,
        plot(x_pts, pdf_scaled_inv_ncx2_40, 'g', 'Linewidth',2.5);  
    end
end