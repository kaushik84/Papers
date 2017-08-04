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
 
clc;
clear all;
close all;

%% Define Constants
d_p = [50];

%% Define parameters
P_p = 1e-3 * 10.^((30:-10:10)/10);  % 10 - 40 dBm
N_p = 10.^((15:-10:15)/10);          % -40 - 10 dB
lambda_p = 1e-5 * [0.01 0.05 0.1 0.5]; 
alpha = (2:0.01:4);

%% Analysing P_s vs lambda_CR, diff d_p, diff lambda_PR and N_p

%% P_s vs N_p 
p_s = zeros(1, length(alpha));
p_a = 1;
h = figure(1);
Fontsize = 18;
axes;

for l = 1:length(lambda_p)
    for n = 1:length(N_p)
        p_s = exp(-2 * p_a * pi^2 * (lambda_p(l) * d_p^2 * N_p(n).^(2./alpha))./(alpha .* sin(2 * pi./alpha)));
        %if p_s < (1 - epsilon_p)
        %    display(strcat('check failed',' alpha = ',num2str(alpha(a)),' lambda_p = ',num2str(lambda_p(l))));
        %end
        %P_s(i,:) = min(10,10 * log10(P_p(1) * (2./(pi^2 * lambda_CR) * sqrt(d_p^(-alpha)/N_p(i)) * log(1 - epsilon_p(e)) - lambda_PT(l)./lambda_CR).^2));
        if n == 1
            c = 'k'
        elseif n == 2
            c = 'g'
        end
            
        if l == 1
            plot(alpha, p_s, strcat(c,'-'),'Linewidth',2);
        elseif l == 2
            plot(alpha, p_s, strcat(c,'--'),'Linewidth',2);
        elseif l == 3
            plot(alpha, p_s, strcat(c,'-.'),'Linewidth',2);
        else 
            plot(alpha, p_s, strcat(c,':') ,'Linewidth',2)
        end
        hold on;
    end
end
%axis([min(lambda_s) max(lambda_s) -35 12]);
axis tight;
grid on;
set(gca,'FontSize',Fontsize);
ylabel('$p_s$','Interpreter','latex','FontSize', 18);
xlabel('$\alpha$','Interpreter','latex','FontSize', 18);
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize',[10 7.5]);
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
print(gcf,'-dpdf','../figures/fig_PT_SucProb_vs_alpha');
print(gcf,'-depsc','../figures/fig_PT_SucProb_vs_alpha');

