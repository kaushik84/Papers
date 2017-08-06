
close all;
h = openfig('fbmc_ofdm_ota_combined3.fig');
pbaspect([1 0.62 1]);
Fontsize = 9;
title('');
xlabel('f [MHz]','FontSize',Fontsize);
ylabel('Normalized PSD [dB/Hz]','FontSize',Fontsize);
hl = legend('OFDM', 'FBMC');
set(hl, 'Location', 'NorthEast', 'FontSize', Fontsize);
set(hl, 'position',[0.665 0.685 0.24 0.15]);
set(gca,'FontSize',Fontsize);
grid on;

laprint(1, 'OFDMvsFBMC', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');