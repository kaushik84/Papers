%%
clc;
close all;
clear all;

% enable (simulation = 1) this to run the simulation, otherwise disable to plot the graphs 
simulation = 1;

if simulation == 1
    %% Parameters

    sim_len = 1e5;                      % Number of Simulation Steps
    K = 1e5;                           % Number of nodes 
    lambda_p = 1e-6;                    % Primary density   
    lambda_s = 1e-5 * [1];              % Secondary density   
    alpha = [4];                        % Path loss exponent
    d_p = [50];                         % avg distance between primary transmitter and receiver
    d_s = [10:5:50];                    % avg distance between secondary transmitter and receiver
    N_p = 10.^(10/10);                  % Primary user constraint 10 dB
    N_s =  10.^((50:-0.1:-20)/10);      % Secondary user constaint 10 dB  
    P_p = 10.^(10/10);                  % Primary tx power 10 dBm

    p_a = 1;                            % Primary receiver access probability
    
    sir = zeros(1,sim_len);

    th_pts = 10;               % Over the simulated range, a linsapce of 10 pts will be chosen     
    
    theta = 0.05;
    
    Exp_Cap_sim = zeros(length(alpha), length(d_s));
    Exp_Cap_th = zeros(length(alpha), length(d_s));
    
    var_Cap_sim = zeros(length(alpha), length(d_s));
    var_Cap_sim = zeros(length(alpha), length(d_s));


    
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
            
%             %%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             for n = 1:sim_len
%                 % Calculate the PPP region for the primary user
%                 R = sqrt(K/lambda_p/pi);
%                 [x,y] = generate_nodes(R,K);
%                 dist_p = sqrt(x.^2+y.^2);
% 
%                 % Calculate the PPP region for the secondary user
%                 R = sqrt(K/lambda_s/pi);
%                 [x,y] = generate_nodes(R,K);
%                 dist_s = sqrt(x.^2+y.^2);
% 
%             %    [val,ind] = min(dist);
%             %    dist(ind) = [];     % Excluding the nearest base station 
% 
%             %     plot(x,y,'*');
%             %     hold on
%             %     plot(0,0,'or');
%             %     plot(d, 0,'*r')
% 
%                 fading_0 = random('exp',1,1,1);
%                 fading = random('exp',1,2,K);    
%                 sir(n) = P_s * d_s(d)^(-alpha(a))*fading_0/(P_p * sum(dist_p.^(-alpha(a)).*fading(1,:)) + P_s * sum(dist_s.^(-alpha(a)).*fading(2,:)));    
%                 
%                 Cap(n) = log2(1 + sir(n));
%             end
% 
%             % Simulated Exp Cap ID 
%             Exp_Cap_sim(a,d) = mean(Cap);
%             var_Cap_sim(a,d) = var(Cap);
            
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
            % Theoretical CDF SIR ID
            
            % Numeric
            b = p_s*(1-theta);
            m = (P_p/(P_s * N_p))^(2/alpha(a)) * (d_s(d)/d_p)^2;
            c = b^m;
            
            % Checkpoint 1 - passed
            fun_exp_test = @(t) 1/log(2) * (1./(1+t)).*c.^(t.^(2/alpha(a)));            
            temp = integral(fun_exp_test,0,Inf);
            
            % Attempt 2 failed
            %Exp_Cap_th(a,d) = - 1/(log(2)) * exp(- 2/alpha(a) * log(c)) * ei(2/alpha(a) * log(c));            
            
            % Checkpoint 2 - passed 
            %fun_exp =  @(t) -1/log(2) * log(c) * c.^t .* log( t.^(alpha(a)/2) + 1);   
            
            % Checkpoint 3 - passed
            %fun_exp =  @(t) 1/log(2) * log((log(c)^(alpha(a)/2) + log(t).^(alpha(a)/2))/log(c)^(alpha(a)/2));            
            
            % Checkpoint 4 - passed
            %fun_exp =  @(t) 1/log(2) * (exp(- log(1/c)* t.^(2/alpha(a))) .* 1./(1 + t));            
            
            % Checkpoint 5 - passed
            % fun_exp =  @(t) log(1/c)/log(2) * (exp(- log(1/c)* t) .* log(1 + t.^(alpha(a)/2)));            
            %Exp_Cap_th(a,d) = integral(fun_exp, 0, Inf);
            
            % Approximation from Checkpoint 4

            Exp_Cap_th(a,d) = -2/log(2) * (cosint(log(1/c))*cos(log(1/c)) + sinint(log(1/c))*sin(log(1/c)));            
            
            fun_B = @(t) log(t.^2+1)./(t.^2+1) .* 2 .* t .* exp(- log(1/c)* t);
            B = integral(fun_B, 0, Inf);
            
            fun_C = @(t) t .* exp(- log(1/c)* (exp(t) - 1 ).^0.5);
            C = integral(fun_C, 0, Inf);
            
            fun_D = @(t) log(1/c)/2 * log(t.^2 + 1).^2 .* exp(- log(1/c) .* t);
            D = integral(fun_D, 0, Inf);
            
            fun_E = @(t) 1/2 * log( log(t).^2/(log(1/c))^2 + 1).^2;
            E = integral(fun_E, 0, 1);
        
            fun_F = @(t) log(1/c)/4 * log(t + 1).^2  .* 1./(t) .* exp(- log(1/c) .* sqrt(t));
            F = integral(fun_D, 0, Inf);
            
            fun_R = @(t) 1/2 * ( log( log(t).^2 + log(1/c)^2 ) - log(log(1/c)^2));
            R = integral(fun_R, 0, 1);
            sum = 0;    
%             x  = 15;
%             for i=0:1000
%                 sum = sum  + ((-log(1/c) * x^0.5)^i)/factorial(i);
%             end
%             actual = exp(-log(1/c)* x^0.5);
%             approx = sum;
            
%             for k=0:0
%                 sum = sum  + ((-log(1/c))^k/factorial(k) * (pi /sin( (k + 2)/2) * ( exp(1) + psi( 1 - (k + 2)/2 ))));
%             end
%             B = 2/(log(2)^2) * sum;

            
            %fun_B = @(t) t .* exp(- log(1/c)* t);
            %B = integral(fun_B, 0, Inf);
            
            %fun_B =  @(t) 2*pi*log(2)* (exp(- log(1/c) * t) - log(1/c) * t .* exp(- log(1/c)* t));
            %B =  integral(fun_B, 0, Inf);
                     
            fun_var = @(t) 2/(log(2)^2) * ((log(1+t))./(1+t)).*c.^(t.^(2/alpha(a)));
            temp =  - Exp_Cap_th(a,d)^2 + integral(fun_var,0,Inf);
            var_Cap_th(a,d) = 2/(log(2))^2 * ( B ) - Exp_Cap_th(a,d)^2;

            
            disp(strcat('Exp_Cap_sim      = ',num2str(Exp_Cap_sim(a,d)))); 
            disp(strcat('Exp_Cap_th       = ',num2str(Exp_Cap_th(a,d)))); 
            
            disp(strcat('var_Cap_sim      = ',num2str(var_Cap_sim(a,d)))); 
            disp(strcat('var_Cap_th       = ',num2str(var_Cap_th(a,d))));
        end
    end
    save('results_PPP_Moments_Cap_ID.mat');
%    quit;
end

load('results_PPP_Moments_Cap_ID.mat')
    
%% Plot the CDFs for the emperical value of Capacity
figure(1);
% Simulation
for a = 1:length(alpha)
    if a == 1
%         plot(d_s, Exp_Cap_sim(a,:), '-','Linewidth',2);
        hold on,
        plot(d_s, Exp_Cap_th(a,:), '*','Linewidth',2);
    elseif a == 2
%         plot(d_s, Exp_Cap_th(a,:), '-.','Linewidth',2);
        hold on,
        plot(d_s, Exp_Cap_th(a,:), '^');
    elseif a == 3
%         plot(d_s, Exp_Cap_th(a,:), '--','Linewidth',2);
        hold on,
        plot(d_s, Exp_Cap_th(a,:), 's');
    else
%         plot(d_s, Exp_Cap_th(a,:), ':','Linewidth',2);
        hold on,
        plot(d_s, Exp_Cap_th(a,:), 'o');
    end
end
Fontsize = 18;
hold on;
grid on;
ylabel('Exp[Capacity]','Interpreter','latex','FontSize', Fontsize + 2);
xlabel('d_s = [m]','Interpreter','latex','FontSize', Fontsize + 2);
axis tight;
set(gca,'FontSize',Fontsize);
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize',[10 7.5]);
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
%print(gcf,'-dpdf','../figures/fig_ID_Exp_Cap');
%print(gcf,'-depsc','../figures/fig_ID_Exp_Cap');

figure(2);
% Simulation
for a = 1:length(alpha)
    if a == 1
       % plot(d_s, var_Cap_sim(a,:), '-','Linewidth',2);
        hold on,
        plot(d_s, var_Cap_th(a,:), '*','Linewidth',2);
    elseif a == 2
       % plot(d_s, var_Cap_sim(a,:), '-.','Linewidth',2);
        hold on,
        plot(d_s, var_Cap_th(a,:), '^');
    elseif a == 3
       % plot(d_s, var_Cap_sim(a,:), '--','Linewidth',2);
        hold on,
        plot(d_s, var_Cap_th(a,:), 's');
    else
       % plot(d_s, var_Cap_sim(a,:), ':','Linewidth',2);
        hold on,
        plot(d_s, var_Cap_th(a,:), 'o');
    end
end
Fontsize = 18;
hold on;
grid on;
ylabel('var[Capacity]','Interpreter','latex','FontSize', Fontsize + 2);
xlabel('d_s = [m]','Interpreter','latex','FontSize', Fontsize + 2);
axis tight;
set(gca,'FontSize',Fontsize);
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize',[10 7.5]);
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
%print(gcf,'-dpdf','../figures/fig_ID_Var_Cap');
%print(gcf,'-depsc','../figures/fig_ID_Var_Cap');

