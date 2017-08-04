close all;
clc;
clear all;

m = 1;
v = 2;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
X = lognrnd(mu,sigma,1,1e6);
%Y_sam = X;
Y_sam = 5.^X;
lognfit(Y_sam)
figure(1)
[F_Y, Y,~,~,eid] = cdfcalc(Y_sam);
F_Y = F_Y(1:length(Y_sam));

plot(Y, F_Y,'k', 'Linewidth',2);
axis([0 400 0 1]);
%% Determine the distribution parmaters using MLE from the emperical path loss values
para_ln = fitdist(Y, 'lognormal');



