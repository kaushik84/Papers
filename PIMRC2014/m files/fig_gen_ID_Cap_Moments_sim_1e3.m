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
    alpha = [2.5 3 3.5 4];                        % Path loss exponent
    d_p = [50];                         % avg distance between primary transmitter and receiver
    d_s_sim = [10:5:50];                    % avg distance between secondary transmitter and receiver
    N_p = 10.^(10/10);                  % Primary user constraint 10 dB
    N_s =  10.^((50:-0.1:-20)/10);      % Secondary user constaint 10 dB  
    P_p = 10.^(10/10);                  % Primary tx power 10 dBm

    p_a = 1;                            % Primary receiver access probability
    
    sir = zeros(1,sim_len);

    th_pts = 10;               % Over the simulated range, a linsapce of 10 pts will be chosen     
    theta = 0.05;
    
    Exp_Cap_sim = zeros(length(alpha), length(d_s_sim));
    var_Cap_sim = zeros(length(alpha), length(d_s_sim));

    
    for a = 1:length(alpha)
        for d = 1:length(d_s_sim)
            % Determine the P_s either theortically or through
            % simulations
            % Lets first determine P_s theoretically
            % Success Probability 
            p_s = exp(-2 * p_a * pi^2 * (lambda_p * d_p^2 * N_p^(2/alpha(a)))/(alpha(a) * sin(2 * pi/alpha(a))));
            epsilon_p = (1 - (1- theta)*p_s);
            % Calculate the P_s for the given configuration of the primary 
            P_s = P_p/(N_p * d_p^(alpha(a))) * power((alpha(a)*sin(2*pi/alpha(a)))./(2 * pi^2 * lambda_s) * log(p_s/(1 - epsilon_p)),alpha(a)/2);

            % Simulate the PPP to determine the SIR CDF at the ID and compare it with the 
            % mentioned in the literature

            disp(strcat('d_s      = ',num2str(d_s_sim(d)))); 
            disp(strcat('alpha    = ',num2str(alpha(a)))); 
            
            %%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            parfor n = 1:sim_len
                % Calculate the PPP region for the primary user
                R = sqrt(K/lambda_p/pi);
                [x,y] = generate_nodes(R,K);
                dist_p = sqrt(x.^2+y.^2);

                % Calculate the PPP region for the secondary user
                R = sqrt(K/lambda_s/pi);
                [x,y] = generate_nodes(R,K);
                dist_s = sqrt(x.^2+y.^2);

                fading_0 = random('exp',1,1,1);
                fading = random('exp',1,2,K);    
                sir(n) = P_s * d_s_sim(d)^(-alpha(a))*fading_0/(P_p * sum(dist_p.^(-alpha(a)).*fading(1,:)) + P_s * sum(dist_s.^(-alpha(a)).*fading(2,:)));    
                
                Cap(n) = log2(1 + sir(n));
            end

            % Simulated Exp Cap ID 
            Exp_Cap_sim(a,d) = mean(Cap);
            var_Cap_sim(a,d) = var(Cap);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Theoretical Evaluation just for comparisions with the
            % Simulations
            b = p_s*(1-theta);
            m = (P_p/(P_s * N_p))^(2/alpha(a)) * (d_s_sim(d)/d_p)^2;
            c = b^m;
            mu = log(1/c);
            % Expectation
            if (alpha(a) == 4)
                Exp_Cap_th = 2/log(2) * (- cosint(mu)*cos(mu) + (pi/2 - sinint(mu))*sin(mu));            
            else
                fun_exp_num = @(t) 1/log(2) * (1./(1+t)).*c.^(t.^(2/alpha(a)));
                Exp_Cap_th = integral(fun_exp_num,0,Inf);
            end
         
            % Variance
            fun_var = @(t) 2/(log(2)^2) * ((log(1+t))./(1+t)).*c.^(t.^(2/alpha(a)));
            var_Cap_th =  - Exp_Cap_th^2 + integral(fun_var,0,Inf);            
         
            disp(strcat('Exp_Cap_sim      = ',num2str(Exp_Cap_sim(a,d)))); 
            disp(strcat('var_Cap_sim      = ',num2str(var_Cap_sim(a,d)))); 
            disp(strcat('Exp_Cap_sim      = ',num2str(Exp_Cap_th))); 
            disp(strcat('var_Cap_sim      = ',num2str(var_Cap_th))); 
        
            
        end
    end
    save('results_PPP_Moments_Cap_ID_1e3.mat');
    quit;
end
load('results_PPP_Moments_Cap_ID_1e3.mat')
    
%% Plot the CDFs for the emperical value of Capacity
figure(1);
% Theoretical Analysis
d_s_th = [10:1:50];
var_Cap_th = zeros(length(alpha), length(d_s_th));
Exp_Cap_th = zeros(length(alpha), length(d_s_th));
for a = 1:length(alpha)
     for d = 1:length(d_s_th)
         p_s = exp(-2 * p_a * pi^2 * (lambda_p * d_p^2 * N_p^(2/alpha(a)))/(alpha(a) * sin(2 * pi/alpha(a))));
        epsilon_p = (1 - (1- theta)*p_s);
        % Calculate the P_s for the given configuration of the primary 
        P_s = P_p/(N_p * d_p^(alpha(a))) * power((alpha(a)*sin(2*pi/alpha(a)))./(2 * pi^2 * lambda_s) * log(p_s/(1 - epsilon_p)),alpha(a)/2);

         b = p_s*(1-theta);
         m = (P_p/(P_s * N_p))^(2/alpha(a)) * (d_s_th(d)/d_p)^2;
         c = b^m;
         mu = log(1/c);
         
         % Expectation
         if (alpha(a) == 4)
                Exp_Cap_th(a,d) = 2/log(2) * (- cosint(mu)*cos(mu) + (pi/2 - sinint(mu))*sin(mu));            
         else
                fun_exp_num = @(t) 1/log(2) * (1./(1+t)).*c.^(t.^(2/alpha(a)));
                Exp_Cap_th(a,d) = integral(fun_exp_num,0,Inf);
         end
         
         % Variance
         fun_var = @(t) 2/(log(2)^2) * ((log(1+t))./(1+t)).*c.^(t.^(2/alpha(a)));
         var_Cap_th(a,d) =  - Exp_Cap_th(a,d)^2 + integral(fun_var,0,Inf);            
     end
end


% Figures
for a = 1:length(alpha)
    if a == 1
        [ax h1 h2] = plotyy(d_s_sim, Exp_Cap_sim(a,:), d_s_sim, var_Cap_sim(a,:));
        set(h1,'LineStyle','*');
        set(h2,'LineStyle','*');
        set(ax,'NextPlot','add');
        plot(ax(1), d_s_th, Exp_Cap_th(a,:),'LineWidth',1.5);
        plot(ax(2), d_s_th, var_Cap_th(a,:),'LineWidth',1.5);
    elseif a == 2
        set(ax,'NextPlot','add');
        plot(ax(1), d_s_sim, Exp_Cap_sim(a,:),'LineStyle','s');
        plot(ax(2), d_s_sim, var_Cap_sim(a,:),'LineStyle','s');
        set(ax,'NextPlot','add');
        plot(ax(1), d_s_th, Exp_Cap_th(a,:),'LineWidth',1.5);
        plot(ax(2), d_s_th, var_Cap_th(a,:),'LineWidth',1.5);
    elseif a == 3
        set(ax,'NextPlot','add');
        plot(ax(1), d_s_sim, Exp_Cap_sim(a,:),'LineStyle','*');
        plot(ax(2), d_s_sim, var_Cap_sim(a,:),'LineStyle','*');
        set(ax,'NextPlot','add');
        plot(ax(1), d_s_th, Exp_Cap_th(a,:),'LineWidth',1.5);
        plot(ax(2), d_s_th, var_Cap_th(a,:),'LineWidth',1.5);
    else
        set(ax,'NextPlot','add');
        plot(ax(1), d_s_sim, Exp_Cap_sim(a,:),'LineStyle','^');
        plot(ax(2), d_s_sim, var_Cap_sim(a,:),'LineStyle','^');
        set(ax,'NextPlot','add');
        plot(ax(1), d_s_th, Exp_Cap_th(a,:),'LineWidth',1.5);
        plot(ax(2), d_s_th, var_Cap_th(a,:),'LineWidth',1.5);
    end
end

Fontsize = 18;
hold on;
grid on;
ylabel(ax(1),'Exp[Capacity]','Interpreter','latex','FontSize',Fontsize+2);
ylabel(ax(2),'var[Capacity]','Interpreter','latex','FontSize',Fontsize+2);
xlabel(ax(1),'$d_s$ = [m]','Interpreter','latex','FontSize',Fontsize+2);
xlabel(ax(2), '$d_s$ = [m]','Interpreter','latex','FontSize',Fontsize+2);
axis(ax(1), 'tight');
axis(ax(2), 'tight');
set(ax(1),'FontSize',Fontsize);
set(ax(2),'FontSize',Fontsize);
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize',[10 7.5]);
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
print(gcf,'-dpdf','../figures/fig_ID_Cap_Moments_vs_d_s_1e3');
print(gcf,'-depsc','../figures/fig_ID_Cap_Moments_vs_d_s_1e3');