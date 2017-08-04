%%
clc;
close all;
clear all;

% enable (simulation = 1) this to run the simulation, otherwise disable to plot the graphs 
simulation = 1;

if simulation == 1
    %% Parameters

    sim_len = 1e3;                      % Number of Simulation Steps
    K = 1e5;                            % Number of nodes 
    lambda_p = 1e-6;                    % Primary density   
    lambda_s = 1e-5 * [1];              % Secondary density   
    alpha = [2.5 3 3.5 4];              % Path loss exponent
    d_p = [50];                         % avg distance between primary transmitter and receiver
    d_s = [20 30 40 50];                % avg distance between secondary transmitter and receiver
    N_p = 10.^(10/10);                  % Primary user interference constraint 10 dB
    R_s =  2;                           % Secondary user capacity constraint 2 bits/sec/Hz  
    P_p = 10.^(10/10);                  % Primary tx power 10 dBm
    
    
    p_a = 1;                            % Primary receiver access probability
    theta = (0.01:0.14:.99);

    sir = zeros(1,sim_len);

    th_pts = 10;               % Over the simulated range, a linsapce of 10 pts will be chosen     
    
    Cap_Out_sim = zeros(length(alpha), length(d_s), length(theta));
    Cap_Out_th = zeros(length(alpha), length(d_s), length(theta));
    epsilon_p = zeros(length(alpha), length(d_s), length(theta));

    
    for a = 1:length(alpha)
        for d = 1:length(d_s)
            for t = 1:length(theta)
                % Determine the P_s either theortically or through
                % simulations

                % Lets first determine P_s theoretically

                % Success Probability 
                p_s = exp(-2 * p_a * pi^2 * (lambda_p * d_p^2 * N_p^(2/alpha(a)))/(alpha(a) * sin(2 * pi/alpha(a))));
                epsilon_p(a,d,t) = (1 - (1- theta(t))*p_s);
                % Calculate the P_s for the given configuration of the primary 
                P_s = P_p/(N_p * d_p^(alpha(a))) * power((alpha(a)*sin(2*pi/alpha(a)))./(2 * pi^2 * lambda_s) * log(p_s/(1 - epsilon_p(a,d,t))),alpha(a)/2);

                % Simulate the PPP to determine the SIR CDF at the ID and compare it with the 
                % mentioned in the literature

                disp(strcat('d_s      = ',num2str(d_s(d)))); 
                disp(strcat('alpha    = ',num2str(alpha(a)))); 
                disp(strcat('theta    = ',num2str(theta(t)))); 
                for n = 1:sim_len
                    % Calculate the PPP region for the primary user
                    R = sqrt(K/lambda_p/pi);
                    [x,y] = generate_nodes(R,K);
                    dist_p = sqrt(x.^2+y.^2);

                    % Calculate the PPP region for the secondary user
                    R = sqrt(K/lambda_s/pi);
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
                    sir(n) = P_s * d_s(d)^(-alpha(a))*fading_0/(P_p * sum(dist_p.^(-alpha(a)).*fading(1,:)) + P_s * sum(dist_s.^(-alpha(a)).*fading(2,:)));    

                    Cap(n) = log2(1 + sir(n));
                end

                % Simulated CDF SIR ID 
                [F_Cap_ID, Cap_sim,~,~,eid] = cdfcalc(Cap);
                F_Cap_ID = F_Cap_ID(2:length(F_Cap_ID));
                
                % Finding the closest value to R_s in Cap_sim and corresponding CDF(closest)
                idx = min(find(Cap_sim >= R_s));
                if ~isempty(idx)
                    Cap_Out_sim(a,d,t) = interp1([Cap_sim(idx - 1) Cap_sim(idx)], [F_Cap_ID(idx - 1) F_Cap_ID(idx)], R_s);
                else
                    Cap_Out_sim(a,d,t) = 1;
                end
                                      
                % Theoretical Cap Outage ID
                Cap_Out_th(a,d,t) = 1 - power( p_s * (1-theta(t)), ...
                power(P_p/P_s * (2.^R_s  - 1)/N_p,  2/alpha(a)) * (d_s(d)/d_p)^2);
            
                disp(strcat('Cap_Out_sim   = ',num2str(Cap_Out_sim(a,d,t)))); 
                disp(strcat('Cap_Out_th    = ',num2str(Cap_Out_th(a,d,t)))); 
            end
        end
    end
    save('results_PPP_Cap_Out_ID_wkr_1e3.mat');
    quit;
end

load('results_PPP_Cap_Out_ID_wkr_1e3.mat')

theta_th = (0.001:0.001:.999);
Cap_Out_th = zeros(length(alpha), length(d_s), length(theta_th));
epsilon_p_th = zeros(length(alpha), length(d_s), length(theta));

% Theoretical Cap Outage ID
for a=1:length(alpha)
    for d=1:length(d_s)
        for t = 1:length(theta_th)
% Success Probability 
                p_s = exp(-2 * p_a * pi^2 * (lambda_p * d_p^2 * N_p^(2/alpha(a)))/(alpha(a) * sin(2 * pi/alpha(a))));
                epsilon_p_th(a,d,t) = (1 - (1- theta_th(t))*p_s);
                % Calculate the P_s for the given configuration of the primary 
                P_s = P_p/(N_p * d_p^(alpha(a))) * power((alpha(a)*sin(2*pi/alpha(a)))./(2 * pi^2 * lambda_s) * log(p_s/(1 - epsilon_p_th(a,d,t))),alpha(a)/2);
            Cap_Out_th(a,d,t) = 1 - power( p_s * (1-theta_th(t)), ...
            power(P_p/P_s * (2.^R_s  - 1)/N_p,  2/alpha(a)) * (d_s(d)/d_p)^2);
        end
    end
end

%% Plot the CDFs for the emperical value of Capacity
figure(1);

thin_factor_sim = 1;  % Describes the thining in the measurements
thin_factor_th = 1;  % Describes the thining in the measurements
a = 2
%for a = 1:length(alpha)
    for d = 1:length(d_s)
        % Thining 
        Cap_Out_sim_thin = reshape(Cap_Out_sim(a,d,:),1,length(theta));
        Cap_Out_sim_thin = Cap_Out_sim_thin(1:length(theta));
        Cap_Out_sim_thin = Cap_Out_sim_thin(1:thin_factor_sim:length(Cap_Out_sim_thin));
        
        epsilon_p_sim_thin = reshape(epsilon_p(a,d,:),1,length(theta));
        epsilon_p_sim_thin = epsilon_p_sim_thin(1:thin_factor_sim:length(epsilon_p_sim_thin));
        
        
        Cap_Out_th_thin = reshape(Cap_Out_th(a,d,:),1,length(theta_th));
        Cap_Out_th_thin = Cap_Out_th_thin(1:thin_factor_th:length(Cap_Out_th_thin));
        
        epsilon_p_th_thin = reshape(epsilon_p_th(a,d,:),1,length(theta_th));
        epsilon_p_th_thin = epsilon_p_th_thin(1:thin_factor_th:length(epsilon_p_th_thin));
        
        if d == 1            
            h2 = plot(epsilon_p_th_thin, Cap_Out_th_thin, '-','Linewidth',1.5);
            hold on,
            h1 = plot(epsilon_p_sim_thin, Cap_Out_sim_thin, 'r*');
            hold on,
        elseif d == 2            
            h2 = plot(epsilon_p_th_thin, Cap_Out_th_thin, '-','Linewidth',1.5);
            hold on,
            h1 = plot(epsilon_p_sim_thin, Cap_Out_sim_thin, 'r^');
            hold on,
        elseif d == 3            
            h2 = plot(epsilon_p_th_thin, Cap_Out_th_thin, '-','Linewidth',1.5);
            hold on,
            h1 = plot(epsilon_p_sim_thin, Cap_Out_sim_thin, 'rs');
            hold on,
        else            
            h2 = plot(epsilon_p_th_thin, Cap_Out_th_thin, '-','Linewidth',1.5);
            hold on,
            h1 = plot(epsilon_p_sim_thin, Cap_Out_sim_thin, 'ro');
            hold on,
        end
    end
%end  

    d = 2;
    for a = 1:length(alpha)
        % Thining 
        Cap_Out_sim_thin = reshape(Cap_Out_sim(a,d,:),1,sim_len);
        Cap_Out_sim_thin = Cap_Out_sim_thin(1:length(theta));
        Cap_Out_sim_thin = Cap_Out_sim_thin(1:thin_factor_sim:length(Cap_Out_sim_thin));
        
        epsilon_p_sim_thin = reshape(epsilon_p(a,d,:),1,length(theta));
        epsilon_p_sim_thin = epsilon_p_sim_thin(1:thin_factor_sim:length(epsilon_p_sim_thin));
        
        
        Cap_Out_th_thin = reshape(Cap_Out_th(a,d,:),1,length(theta_th));
        Cap_Out_th_thin = Cap_Out_th_thin(1:thin_factor_th:length(Cap_Out_th_thin));
        
        epsilon_p_th_thin = reshape(epsilon_p_th(a,d,:),1,length(theta_th));
        epsilon_p_th_thin = epsilon_p_th_thin(1:thin_factor_th:length(epsilon_p_th_thin));
        
        if a == 1            
            h2 = plot(epsilon_p_th_thin, Cap_Out_th_thin, '--','Linewidth',1.5);
            hold on,
            h1 = plot(epsilon_p_sim_thin, Cap_Out_sim_thin, 'r^');
            hold on,
        elseif a == 2            
%             h2 = plot(epsilon_p_th_thin, Cap_Out_th_thin, '-','Linewidth',1.5);
%             hold on,
%             h1 = plot(epsilon_p_sim_thin, Cap_Out_sim_thin, 'r^');
%             hold on,
        elseif a == 3            
            h2 = plot(epsilon_p_th_thin, Cap_Out_th_thin, '-.','Linewidth',1.5);
            hold on,
            h1 = plot(epsilon_p_sim_thin, Cap_Out_sim_thin, 'r^');
            hold on,
        else            
            h2 = plot(epsilon_p_th_thin, Cap_Out_th_thin, ':','Linewidth',1.5);
            hold on,
            h1 = plot(epsilon_p_sim_thin, Cap_Out_sim_thin, 'r^');
            hold on,
        end
    end


Fontsize = 18;
hold on;
% Theoretical
%CDF_sir_th = 1 - exp(-2 * pi^2 * (lambda * d^2 * (10.^(linspace(-20,40,20)/10)).^(2./alpha))./(alpha .* sin(2 * pi./alpha)));
%plot((linspace(-20,40,20)), CDF_sir_th,'*')
grid on;
ylabel('$\psi_s$','Interpreter','latex','FontSize', Fontsize + 2);
xlabel('$\epsilon_p$','Interpreter','latex','FontSize', Fontsize + 2);
hl = legend(strcat('the $d_s$ = ',num2str(d_s(1)),', $\alpha$ = ',num2str(alpha(2))),...
    strcat('sim $d_s$ = ',num2str(d_s(1)),', $\alpha$ = ',num2str(alpha(2))),... 
    strcat('the $d_s$ = ',num2str(d_s(2)),', $\alpha$ = ',num2str(alpha(2))),...
    strcat('sim $d_s$ = ',num2str(d_s(2)),', $\alpha$ = ',num2str(alpha(2))),...
    strcat('the $d_s$ = ',num2str(d_s(3)),', $\alpha$ = ',num2str(alpha(2))),...
    strcat('sim $d_s$ = ',num2str(d_s(3)),', $\alpha$ = ',num2str(alpha(2))),...
    strcat('sim $d_s$ = ',num2str(d_s(4)),', $\alpha$ = ',num2str(alpha(2))),... 
    strcat('the $d_s$ = ',num2str(d_s(4)),', $\alpha$ = ',num2str(alpha(2))),...
    strcat('sim $d_s$ = ',num2str(d_s(2)),', $\alpha$ = ',num2str(alpha(1))),...
    strcat('the $d_s$ = ',num2str(d_s(2)),', $\alpha$ = ',num2str(alpha(1))),...
    strcat('sim $d_s$ = ',num2str(d_s(2)),', $\alpha$ = ',num2str(alpha(3))),...
    strcat('sim $d_s$ = ',num2str(d_s(2)),', $\alpha$ = ',num2str(alpha(3))),...
    strcat('the $d_s$ = ',num2str(d_s(2)),', $\alpha$ = ',num2str(alpha(3))),...
    strcat('sim $d_s$ = ',num2str(d_s(2)),', $\alpha$ = ',num2str(alpha(4))),...
    'Location','NorthEast');
set(hl,'Interpreter','latex','FontSize', Fontsize + 2);
axis tight;
axis([0 1 0 1]);
set(gca,'FontSize',Fontsize);
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize',[10 7.5]);
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
print(gcf,'-dpdf','../figures/fig_ID_Cap_Out_vs_epsilon_wkr_1e3');
print(gcf,'-depsc','../figures/fig_ID_Cap_Out_vs_epsilon_wkr_1e3');
