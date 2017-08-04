%%
clc;
close all;
clear all;

% enable (simulation = 1) this to run the simulation, otherwise disable to plot the graphs 
simulation = 1;

if simulation == 1
    %% Parameters

    sim_len = 1e4;                      % Number of Simulation Steps
    K = 40e3;                           % Number of nodes 
    lambda_p = 1e-6;                    % Primary density   
    lambda_s = 1e-5 * (0.05:0.01:2.5);  % Secondary density   
    alpha = [2.5 3 3.5 4];              % Path loss exponent
    d_p = [50];                         % avg distance between transmitter and receiver
    N_p = 10.^(10/10);                  % Primary user constraint 10 dB
    P_p = 10.^(10/10);                  % 10 dBm

    p_a = 1;                            % Primary receiver access probability

    sir = zeros(1,sim_len);

    P_s_sim = zeros(length(alpha),length(lambda_s));
    P_s_th = zeros(length(alpha),length(lambda_s));
    theta = 0.05;
    
    for d = 1:length(d_p)
        for a = 1:length(alpha)
            for l=1:length(lambda_s) 
                disp(strcat('d_p      = ',num2str(d_p(d)))); 
                disp(strcat('alpha    = ',num2str(alpha(a)))); 
                disp(strcat('lambda_s = ',num2str(lambda_s(l)))); 
                for n = 1:sim_len
                    % Calculate the PPP region for the primary user
                    R = sqrt(K/lambda_p/pi);
                    [x,y] = generate_nodes(R,K);    
                    dist_p = sqrt(x.^2+y.^2);

                    % Calculate the PPP region for the secondary user
                    R = sqrt(K/lambda_s(l)/pi);
                    [x,y] = generate_nodes(R,K);
                    dist_s = sqrt(x.^2+y.^2);

                %    [val,ind] = min(dist);
                %    dist(ind) = [];     % Excluding the nearest base station 

                %     plot(x,y,'*');
                %     hold on
                %     plot(0,0,'or');
                %     plot(d, 0,'*r')

                    fading_0 = random('exp',1,1,1);
                    fading = random('exp',1,2,K);    

                    signal(n) = d_p(d)^(-alpha(a))*fading_0;
                    interference_p(n) = sum(dist_p.^(-alpha(a)).*fading(1,:));
                    interference_s(n) = sum(dist_s.^(-alpha(a)).*fading(2,:));    
                    % Sum the interference greater that N_p      

                end
                
                % success probability
                sir_p = P_p * signal./(P_p * interference_p);
                [CDF_sir_p, log_sir_p,~,~,eid] = cdfcalc(10*log10(sir_p));
                CDF_sir_p = CDF_sir_p(2:length(CDF_sir_p));

                
                % Finding the closed value to N_p in log_sir and corresponding CDF(closest)
                temp = abs(log_sir_p - 10* log10(N_p));
                [idx idx] =  min(temp);
                p_s_sim = 1 - CDF_sir_p(idx);
                p_s_th = exp(-2 * p_a * pi^2 * (lambda_p * d_p(d)^2 * N_p^(2/alpha(a)))/(alpha(a) * sin(2 * pi/alpha(a))));

                for P_s_dB = 50:-0.01:-100     % dB scale
                    % success probability
                    sir = P_p * signal./(P_p * interference_p + 10^(P_s_dB/10) * interference_s);
                    [CDF_sir, log_sir,~,~,eid] = cdfcalc(10*log10(sir));

                    % Finding the closed value to N_p in log_sir and corresponding CDF(closest)
                    temp = abs(log_sir - 10* log10(N_p));
                    [idx idx] =  min(temp);
                    CDF_sir(idx);
                    if CDF_sir(idx) <= 1 - (1- theta)*p_s_sim        
                        % time to stop
                        P_s_sim(a,l) = P_s_dB;
                        break;
                    end    
                end
                P_s_th(a,l) = 10 * log10(P_p/(N_p*d_p(d)^alpha(a))*(alpha(a) * sin (2* pi/alpha(a))/(2 * pi^2 * lambda_s(l)) *  log(1/(1 -theta)))^(alpha(a)/2));
                disp(strcat('P_s_th   = ',num2str(P_s_th(a,l)))); 
                disp(strcat('P_s_sim  = ',num2str(P_s_sim(a,l)))); 
            end
        end
    end
    save('results_PPP_CR_TransmitPowerVs_CR_lambda.mat');
    quit;
end

load('results_PPP_CR_TransmitPowerVs_CR_lambda.mat')
    
%% Plot the CDFs for the emperical value of Capacity
figure(1);
[CDF_sir, log_sir,~,~,eid] = cdfcalc(10*log10(sir));
CDF_sir = CDF_sir(1:length(log_sir));
plot(log_sir, CDF_sir, 'b', 'Linewidth', 2); 

[CDF_sir_p, log_sir_p,~,~,eid] = cdfcalc(10*log10(sir_p));
CDF_sir_p = CDF_sir_p(1:length(log_sir));
plot(log_sir_p, CDF_sir_p, 'r', 'Linewidth', 2); 

% Simulation
plot(lambda_s, P_s_sim, 'b', 'Linewidth', 2);
axis([-20 40 0 1]);

Fontsize = 18;
hold on;
% Theoretical
%CDF_sir_th = 1 - exp(-2 * pi^2 * (lambda * d^2 * (10.^(linspace(-20,40,20)/10)).^(2./alpha))./(alpha .* sin(2 * pi./alpha)));
%plot((linspace(-20,40,20)), CDF_sir_th,'*')
plot(lambda_s, P_s_th, 'or');
grid on;
ylabel('$P_s$ (dBm)','Interpreter','latex','FontSize', Fontsize + 2);
xlabel('$\lambda_s$','Interpreter','latex','FontSize', Fontsize + 2);
axis tight;
set(gca,'FontSize',Fontsize);
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize',[10 7.5]);
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
print(gcf,'-dpdf','../figures/fig_CR_Txp_vs_lambda_CR');
print(gcf,'-depsc','../figures/fig_CR_Txp_vs_lambda_CR');
