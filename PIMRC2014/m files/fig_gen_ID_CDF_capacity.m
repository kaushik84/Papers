%% The script analyzes the CDF of Capacity at the Indoor device over different parameters.
% F_SIR_ID = 1 - e^(-a * sqrt(2^R_s -1)) 
% where a =  pi^2 / (2 sqrt(d_s^(-alpha))*((sqrt(P_p)* lambda_PT + sqrt(P_s) * lambda_CR)/sqrt(P_s))
% F_SIR_ID = F(R_s, lambda_CR, lambda_PT, d_s)
% P_s = Secondary transmitter or Cognitive Relay tranmsit power
% P_p = Primary transmitter tranmsit power
% alpha = path loss constant
% d_p = distance from the reference primary transmitter to primary receiver
% d_s = distance from the reference secodary transmitter (CR) to secondar receiver (ID)
% lambda_PT = density of the primary transmitters
% lambda_CR = density of the secondary transmitters or CR
 
clc;
clear all;
close all;

%% Define Constants
alpha = 4;
epsilon_p = [0.01]; 
d_p = [25];
d_s = [12];

%% Define parameters
P_p = 1e-3 * 10.^((10:-10:10)/10);  % 10 - 40 dBm
N_p = 10.^((5:-5:5)/10);       % -40 - 10 dB
lambda_PT = 1e-5 * [0.1 0.5 1]; 
lambda_CR = 1e-5 * [0.1:0.4:1];
R_s = (0:0.1:20);

%% Analysing P_s vs lambda_CR, diff d_p, diff lambda_PR and N_p

%% F_SIR_ID vs R_s for diff lambda_PR  
F_SIR_ID = zeros(1, length(R_s));
h = figure(1);
Fontsize = 18;
axes;
for l_CR = 1:length(lambda_CR)
    for l_PT = 1:length(lambda_PT)
        for e = 1:length(epsilon_p)
            for i=1:length(N_p)
                % Calculate the P_s for the given configuration of the primary 
                % A quick Check up 
                if d_p > pi * power(N_p(i)/(4 * lambda_PT(l_PT)^2) * (log(1/(1 - epsilon_p(e))))^2 ,.25)
                    display(strcat('Check failed with parameters:',' N_p = ',num2str(N_p(i)),' epsilon_p = ', epsilon_p(e),' lambda_PT = ',num2str(lambda_PT(l_PT))));
                end
                P_s = min(0.001, P_p(1) * (2./(pi^2 * lambda_CR(l_CR)) * sqrt(d_p^(-alpha)/N_p(i)) * log(1/(1 - epsilon_p(e))) - lambda_PT(l_PT)./lambda_CR(l_CR)).^2);
                P_s
                % Calculate F_SIR_ID
                a = pi^2 / (2 * sqrt(d_s^(-alpha)))*((sqrt(P_p(1))* lambda_PT(l_PT) + sqrt(P_s) * lambda_CR(l_CR))/sqrt(P_s));
                F_SIR_ID = 1 - exp(-a * sqrt(2.^R_s - 1));
                if P_s == 0.1
                    display(strcat('Limit Reachad, lambda_CR = ',num2str(lambda_CR(l_CR))));
                    c = 'r';
                else
                    c = 'k';
                end
                if l_PT == 1
                    plot(R_s, F_SIR_ID, strcat(c,'-'),'Linewidth',2);
                elseif l_PT == 2
                    plot(R_s, F_SIR_ID, strcat(c,'-.'),'Linewidth',2);
                else
                    plot(R_s, F_SIR_ID, strcat(c,'--'),'Linewidth',2);
                end
                hold on;
            end
        end
    end
end
axis tight;
grid on;
set(gca,'FontSize',Fontsize);
ylabel('CDF','Interpreter','latex','FontSize', 18);
xlabel('$R_s$','Interpreter','latex','FontSize', 18);
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize',[10 7.5]);
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
print(gcf,'-dpdf','../figures/fig_ID_CDF_capacity');
print(gcf,'-depsc','../figures/fig_ID_CDF_capacity');

