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
d_p = [100];
d_s = [25 50];

%% Define parameters
P_p = 1e-3 * 10.^((30:-5:10)/10);    % 30 - 10 dBm
N_p = 10.^((10:-10:10)/10);         % -40 - 10 dB
p_a = 1; 
N_s =  10.^((50:-0.1:-20)/10);           
lambda_p = 1e-7 * [1];
lambda_s = 1e-6 * [1];
alpha = [2.5 3 3.5 4];
theta = [0.05];

%% F_SIR_ID vs N_s for diff lambda_PR  
F_SIR_ID = zeros(1, length(N_s));
h = figure(1);
Fontsize = 18;
axes;
for l = 1:length(lambda_s)
    for ddp = 1:length(d_p)
        for dds = 1:length(d_s)
            for a =1:length(alpha)
                % Success Probability 
                p_s = exp(-2 * p_a * pi^2 * (lambda_p * d_p(ddp)^2 * N_p^(2/alpha(a)))/(alpha(a) * sin(2 * pi/alpha(a))));
                p_s
                epsilon_p = (1 - (1- theta)*p_s);
                % Calculate the P_s for the given configuration of the primary 
                P_s = P_p(1)/(N_p * d_p(ddp)^(alpha(a))) * power((alpha(a)*sin(2*pi/alpha(a)))./(2 * pi^2 * lambda_s(l)) * log(p_s/(1 - epsilon_p)),alpha(a)/2);

                temp = p_s * (1 - theta);
                % Calculate F_SIR_ID, inserting P_s 
                F_SIR_ID = ones(1, length(N_s)) - power( power(p_s, power(P_p/P_s, 2/alpha(a))) * P_s/P_p * 1/(1-theta), power(N_s/N_p, 2/alpha(a)) * (d_s(dds)/d_p(ddp))^2);
                if l == 1
                    c = 'k';    
                else
                    c = 'k';
                end
                if a == 1
                    plot(10*log10(N_s), F_SIR_ID, strcat(c,'-'),'Linewidth',2);
                elseif a == 2
                    plot(10*log10(N_s), F_SIR_ID, strcat(c,'-.'),'Linewidth',2);
                elseif a == 3
                    plot(10*log10(N_s), F_SIR_ID, strcat(c,'--'),'Linewidth',2);
                else
                    plot(10*log10(N_s), F_SIR_ID, strcat(c,':'),'Linewidth',2);
                end
                hold on;
            end
        end
    end
end

axis([min(10 * log10(N_s)) max(10 * log10(N_s)) 0 1]);
grid on;
%h = legend('$\alpha$ = 2.5','$\alpha$ = 3','$\alpha$ = 3.5','$\alpha$ = 4');
%set(h,'Interpreter','latex','location','NorthWest');
set(gca,'FontSize',Fontsize);
ylabel('CDF','Interpreter','latex','FontSize', 18);
xlabel('SIR','Interpreter','latex','FontSize', 18);
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize',[10 7.5]);
set(gcf, 'PaperPositionMode','manual');
set(gcf, 'PaperPosition',[ 0 0 10 7.5]);
print(gcf,'-dpdf','../figures/fig_ID_CDF_SIR');
print(gcf,'-depsc','../figures/fig_ID_CDF_SIR');

