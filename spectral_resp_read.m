clear;
clc;
close all;


if isunix
    addpath('/home/felipe/Dropbox/VLC_exp/data/spectral_resp/');
else
    addpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\spectral_resp');
end



spect = dlmread('white_LED_diffused.csv');


plot(spect(:,1),spect(:,2)./max(spect(:,2)));
ylabel('Normalized Luminous Intensity [cd]','interpreter','latex');
xlabel('Wavelength [nm]','interpreter','latex');
xlim([min(spect(:,1)) 700]);

if isunix
    rmpath('/home/felipe/Dropbox/VLC_exp/data/spectral_resp/');
else
    rmpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\spectral_resp');
end