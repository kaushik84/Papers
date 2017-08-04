%% The script analyzes the CDF of Capacity at the Indoor device over different parameters.
% Exp_C = - 2/log(2) (ci(a)cos(a) + si(a)sin(a)) 
% where a =  pi^2 / (2 sqrt(d_s^(-alpha))*((sqrt(P_p)* lambda_PT + sqrt(P_s) * lambda_CR)/sqrt(P_s))
% Exp_C = F(R_s, lambda_CR, lambda_PT, d_s)
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
epsilon_p = [0.05]; 
d_p = [25];
d_s = [6 12];

%% Define parameters
P_p = 1e-3 * 10.^((10:-10:10)/10);  % 10 - 40 dBm
N_p = 10.^((5:-10:5)/10);       % -40 - 10 dB
lambda_PT = 1e-5 * [0.1 0.5 2.5]; 
lambda_CR = 1e-5 * (0.1:0.01:2);

%% Analysing P_s vs lambda_CR, diff d_p, diff lambda_PR and N_p

%% F_SIR_ID vs R_s for diff lambda_PR  
Exp_C = zeros(1, length(lambda_CR));
h = figure(1);
Fontsize = 18;
axes;
for l_PT = 1:length(lambda_PT)
    for e = 1:length(epsilon_p)
        for i=1:length(N_p)
            for d = 1:length(d_s)
                if d_p > power(N_p(i)/(4 * lambda_PT(l_PT)^2) * (log(1/(1 - epsilon_p(e))))^2 ,.25)
                    power(N_p(i)/(4 * lambda_PT(l_PT)^2) * (log(1/(1 - epsilon_p(e))))^2 ,.25)
                    display(strcat('Check failed with parameters:',' N_p = ',num2str(N_p(i)),' epsilon_p = ', epsilon_p(e),' lambda_PT = ',num2str(lambda_PT(l_PT))));
                end
                % Calculate the P_s for the given configuration of the primary 
                P_s = P_p(1) * (2./(pi^2 * lambda_CR) * sqrt(d_p^(-alpha)/N_p(i)) * log(1/(1 - epsilon_p(e))) - lambda_PT(l_PT)./lambda_CR).^2;

                % Calculate a
                a = pi^2 / (2 * sqrt(d_s(d)^(-alpha)))*((sqrt(P_p(1))* lambda_PT(l_PT) + sqrt(P_s) .* lambda_CR)./sqrt(P_s));
                % Calculate expected Capacity 
                Exp_C = - 2/log(2)*(cosint(a).*cos(a) + sinint(a).*sin(a));
                if e == 1
                    c = 'k';
                else 
                    c = 'r';
                end
                if d == 1
                    plot(lambda_CR, Exp_C, strcat(c,'-'),'Linewidth',2);
                elseif d == 2
                    plot(lambda_CR, Exp_C, strcat(c,'-.'),'Linewidth',2);
                else
                    plot(lambda_CR, Exp_C, strcat(c,'--'),'Linewidth',2);
                end
                hold on;
            end
        end
    end
end
axis tight;
grid on;
set(gca,'FontSize',Fontsize);
ylabel('C (bits/sec/Hz)','Interpreter','latex','FontSize', 18);
xlabel('$\lambda_{CR}$','Interpreter','latex','FontSize', 18);
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize',[10 7.5]);
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
print(gcf,'-dpdf','../figures/fig_ID_Exp_capacity');
print(gcf,'-depsc','../figures/fig_ID_Exp_capacity');
