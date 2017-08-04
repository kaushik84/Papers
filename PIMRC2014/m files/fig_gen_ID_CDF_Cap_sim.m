%%
clc;
close all;
clear all;

% enable (simulation = 1) this to run the simulation, otherwise disable to plot the graphs 
simulation = 1;

if simulation == 1
    %% Parameters

    sim_len = 1e3;                      % Number of Simulation Steps
    K = 40e3;                           % Number of nodes 
    lambda_p = 1e-6;                    % Primary density   
    lambda_s = 1e-5 * [1];              % Secondary density   
    alpha = [2.5];              % Path loss exponent
    d_p = [50];                         % avg distance between primary transmitter and receiver
    d_s = [25];                  % avg distance between secondary transmitter and receiver
    N_p = 10.^(10/10);                  % Primary user constraint 10 dB
    N_s =  10.^((50:-0.1:-20)/10);      % Secondary user constaint 10 dB  
    P_p = 10.^(10/10);                  % Primary tx power 10 dBm

    p_a = 1;                            % Primary receiver access probability
    
    sir = zeros(1,sim_len);

    th_pts = 10;               % Over the simulated range, a linsapce of 10 pts will be chosen     
    F_Cap_ID_sim = zeros(length(alpha), length(d_s), sim_len);
    F_Cap_ID_th = zeros(length(alpha), length(d_s), th_pts);   
    
    Cap_sim = zeros(length(alpha), length(d_s), sim_len);
    Cap_th = zeros(length(alpha), length(d_s), th_pts);
    theta = 0.05;
    
    for a = 1:length(alpha)
        for d = 1:length(d_s)
            % Determine the P_s either theortically or through
            % simulations

            % Lets first determine P_s theoretically

            % Success Probability 
            p_s = exp(-2 * p_a * pi^2 * (lambda_p * d_p^2 * N_p^(2/alpha(a)))/(alpha(a) * sin(2 * pi/alpha(a))));
            p_s
            epsilon_p = (1 - (1- theta)*p_s);
            % Calculate the P_s for the given configuration of the primary 
            P_s = P_p/(N_p * d_p^(alpha(a))) * power((alpha(a)*sin(2*pi/alpha(a)))./(2 * pi^2 * lambda_s) * log(p_s/(1 - epsilon_p)),alpha(a)/2);

            % Simulate the PPP to determine the SIR CDF at the ID and compare it with the 
            % mentioned in the literature

            disp(strcat('d_s      = ',num2str(d_s(d)))); 
            disp(strcat('alpha    = ',num2str(alpha(a)))); 
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

            % Sinmulated CDF SIR ID 
            [F_Cap_ID, Cap_sim(a,d,:),~,~,eid] = cdfcalc(Cap);
            F_Cap_ID_sim(a,d,:) = F_Cap_ID(2:length(F_Cap_ID));
            
            % Theoretical CDF SIR ID
            Cap_th(a,d,:) = linspace(min(Cap_sim(a,d,:)),max(Cap_sim(a,d,:)),th_pts);
            
            % Step 1
            %F_Cap_ID_th(a,d,:) = ones(1, 1, th_pts) - power( p_s * (1-theta), ...
            %power(P_p/P_s * (2.^Cap_th(a,d,:)  - 1)/N_p,  2/alpha(a)) * (d_s(d)/d_p)^2);
        
            % Step 2
            b = p_s * (1 - theta);
            m = P_p / (N_p * P_s) * (d_s(d)/d_p)^(alpha);
            
            % no Rule 
            %F_Cap_ID_th(a,d,:) =  ones(1, 1, th_pts) - b.^(((2.^Cap_th(a,d,:) - 1) * m).^(2/alpha(a)));

            
            % Rule e^(a*b) = (e^b)^a
            F_Cap_ID_th(a,d,:) =  ones(1, 1, th_pts) - b.^(2/alpha(a)).^((2.^Cap_th(a,d,:) - 1)*m);

            % Rule e^a^b  = e^b^a
            %c = (power(b^m, 2/alpha(a)))^2;
            %F_Cap_ID_th(a,d,:) =  ones(1, 1, th_pts) - ((b^m)).^(2/alpha(a)).^(2.^Cap_th(a,d,:) - 1);

            %n = power(b^m, 2/alpha(a));
            %F_Cap_ID_th(a,d,:) =  ones(1, 1, th_pts) - c.^Cap_th(a,d,:)/n;
        end
    end
    save('results_PPP_CDF_Cap_ID.mat');
%    quit;
end

load('results_PPP_CDF_Cap_ID.mat')
    
%% Plot the CDFs for the emperical value of Capacity
figure(1);
% Simulation
for a = 1:length(alpha)
    for d = 1:length(d_s)
        if a == 1
            plot(reshape(Cap_sim(a,d,:), 1, sim_len), reshape(F_Cap_ID_sim(a,d,:), 1, sim_len), '-','Linewidth',2);
            hold on,
            plot(reshape(Cap_th(a,d,:), 1, th_pts), reshape(F_Cap_ID_th(a,d,:), 1, th_pts), 'r*');
        elseif a == 2
            plot(reshape(Cap_sim(a,d,:), 1, sim_len), reshape(F_Cap_ID_sim(a,d,:), 1, sim_len), '-.','Linewidth',2);
            hold on,
            plot(reshape(Cap_th(a,d,:), 1, th_pts), reshape(F_Cap_ID_th(a,d,:), 1, th_pts), 'r^');
        elseif a == 3
            plot(reshape(Cap_sim(a,d,:), 1, sim_len), reshape(F_Cap_ID_sim(a,d,:), 1, sim_len), '--','Linewidth',2);
            hold on,
            plot(reshape(Cap_th(a,d,:), 1, th_pts), reshape(F_Cap_ID_th(a,d,:), 1, th_pts), 'rs');
        else
            plot(reshape(Cap_sim(a,d,:), 1, sim_len), reshape(F_Cap_ID_sim(a,d,:), 1, sim_len), ':','Linewidth',2);
            hold on,
            plot(reshape(Cap_th(a,d,:), 1, th_pts), reshape(F_Cap_ID_th(a,d,:), 1, th_pts), 'ro');
        end
    end
end

% Theory

Fontsize = 18;
hold on;
% Theoretical
%CDF_sir_th = 1 - exp(-2 * pi^2 * (lambda * d^2 * (10.^(linspace(-20,40,20)/10)).^(2./alpha))./(alpha .* sin(2 * pi./alpha)));
%plot((linspace(-20,40,20)), CDF_sir_th,'*')
grid on;
ylabel('CDF','Interpreter','latex','FontSize', Fontsize + 2);
xlabel('Capacity at ID','Interpreter','latex','FontSize', Fontsize + 2);
axis tight;
axis([0 12 0 1]);
set(gca,'FontSize',Fontsize);
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize',[10 7.5]);
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
print(gcf,'-dpdf','../figures/fig_ID_CDF_Cap');
print(gcf,'-depsc','../figures/fig_ID_CDF_Cap');
