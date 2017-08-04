%% The script analyzes the transmit power of the Cognitive Relay over different parameters.
% lambda_CR = sqrt(P_p/P_s) * (2/lambda_CR \sqrt(d_p^(-alpha)/N_P) ln(1 - epsilon_p) + lambda_PT sqrt(N_p/d_p^(1-alpha))) 
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
alpha = 4;
epsilon_p = [0.01]; 
d_p = [50];

%% Define parameters
P_p = 1e-3 * 10.^((10:-10:10)/10);  % 10 - 40 dBm
N_p = 10.^((10:-10:0)/10);       % -40 - 10 dB
lambda_PT = 1e-5 * [0.5 1 2 4]; 
P_s = 1e-3 * 10.^((-20:1:0)/10);

%% Analysing P_s vs lambda_CR, diff d_p, diff lambda_PR and N_p

%% P_s vs N_p 
lambda_CR = zeros(length(N_p), length(P_s));
h = figure(1);
Fontsize = 18;
axes;
for l = 1:length(lambda_PT)
    for e = 1:length(epsilon_p)
        for i=1:length(N_p)
            lambda_CR(i,:) = abs(sqrt(P_p(1)./P_s) * (2./(pi^2) * sqrt(d_p^(-alpha)/N_p(i)) * log(1 - epsilon_p(e)) - lambda_PT(l)));
            if i == 1
                plot(10 * log10(P_s/10^-3), lambda_CR(i,:),  'k-','Linewidth',2);
            elseif i == 2
                plot(10 * log10(P_s/10^-3), lambda_CR(i,:),  'k-.','Linewidth',2);
            else
                plot(10 * log10(P_s/10^-3), lambda_CR(i,:),  'k--','Linewidth',2);
            end
            hold on;
        end
    end
end
axis tight;
grid on;
set(gca,'FontSize',Fontsize);
xlabel('$P_s$ (dBm)','Interpreter','latex','FontSize', 18);
ylabel('$\lambda_{CR}$','Interpreter','latex','FontSize', 18);
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize',[10 7.5]);
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
print(gcf,'-dpdf','../figures/fig_CR_Txp_vs_lambda_CR');
print(gcf,'-depsc','../figures/fig_CR_Txp_vs_lambda_CR');

