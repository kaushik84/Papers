%% The script analyzes the transmit power of the Cognitive Relay over different parameters.
% 1) P_s = P_p * (2/lambda_CR 2 \sqrt(d_p^(-alpha)/N_P) ln(1 - epsilon_p) + lambda_PT \sqrt(N_p/d_p^(1-alpha)))^2  
% 2) P_s = P_p * (d_p^(alpha)/(2 * pi^2 * lambda_s^(alpha/2) N_p) * log (p_s/1 - epsilon_p)
% where p_s = exp^(-2 pi^2 * (lambda_p * d_p^2 * N_p^(2/alpha))/(alpha * sin(2 pi/alpha)))
% P_s = F(P_p, lambda_CR, lambda_PT, d_p, N_p)
% P_s = Secondary transmitter or Cognitive Relay tranmsit power
% P_p = Primary transmitter tranmsit power
% alpha = path loss constant
% d_p = distance from the reference primary transmitter to primary receiver
% N_p = SIR outage theshold
% epsilon = SIR outage probability
% lambda_PT = density of the primary transmitters
% lambda_CR = density of the secondary transmitters or CR
% p_a =  access probability of the primary user network
 
clc;
clear all;
close all;

%% Define Constants
d_p = [100];                
d_s = [25 50];

%% Define parameters
P_p = 10.^(10/10);                   % 10 dBm
N_p = 10.^((10:-10:10)/10);          % 10 dB
p_a = 1;                                  
lambda_p = 1e-5 * [1];
lambda_s = 1e-5 * (0.05:0.01:2.5);
alpha = [2.5 3 3.5 4];
theta = [0.05];
%% Analysing P_s vs lambda_CR, diff d_p, diff lambda_PR and N_p

%% P_s vs N_p 
P_s = zeros(1, length(lambda_s));
h = figure(1);
Fontsize = 17;
axes;
for l = 1:length(lambda_p)
    for d = 1:length(d_p)
        for a =1:length(alpha)
            % Success Probability 
            p_s = exp(-2 * p_a * pi^2 * (lambda_p(l) * d_p(d)^2 * N_p^(2/alpha(a)))/(alpha(a) * sin(2 * pi/alpha(a))));
            p_s 
            epsilon_p = (1 - (1- theta)*p_s);
            r = log(p_s/(1 - epsilon_p))
            if p_s < (1 - epsilon_p)
                display(strcat('check failed',' alpha = ',num2str(alpha(a)),' lambda_p = ',num2str(lambda_p(l))));
            end
            P_s = 10 * log10(P_p(1)/(N_p * d_p(d)^(alpha(a))) *power((alpha(a)*sin(2*pi/alpha(a)))./(2 * pi^2 * lambda_s) * log(p_s/(1 - epsilon_p)),alpha(a)/2));
            %P_s(i,:) = min(10,10 * log10(P_p(1) * (2./(pi^2 * lambda_CR) * sqrt(d_p^(-alpha)/N_p(i)) * log(1 - epsilon_p(e)) - lambda_PT(l)./lambda_CR).^2));
            if a == 1
                plot(lambda_s, P_s, 'k-','Linewidth',2);
            elseif a == 2
                plot(lambda_s, P_s, 'k-.','Linewidth',2);
            elseif a == 3
                plot(lambda_s, P_s, 'k--','Linewidth',2);
            else 
                plot(lambda_s, P_s, 'k:' ,'Linewidth',2);
            end
            hold on;
        end
    end
end
%axis([min(lambda_s) max(lambda_s) -35 12]);

grid on;
set(gca,'FontSize',Fontsize);
h = legend('$\alpha$ = 2.5',...
       '$\alpha$ = 3',...
       '$\alpha$ = 3.5',...
       '$\alpha$ = 4');
set(h,'Interpreter','latex');
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

