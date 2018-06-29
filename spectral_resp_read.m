clear;
clc;
close all;

if isunix
    addpath('/home/felipe/Dropbox/VLC_exp/data/spectral_resp/');
else
    addpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\spectral_resp');
end

addpath(['.' filesep 'Misc' filesep]); 

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]); 

spect = dlmread('white_LED_diffused.csv');

plot(spect(:,1),spect(:,2)./max(spect(:,2)));
ylabel('Normalized Luminous Intensity','interpreter','latex');
% xlabel('Wavelength [nm]','interpreter','latex');
xlim([min(spect(:,1)) 700]);

set(gca,'xtick',[])
set(gca,'ytick',[])
box off
formatFig( gcf ,['.' filesep 'figs' filesep 'spectral_response'], 'en' , figProp );

if isunix
    rmpath('/home/felipe/Dropbox/VLC_exp/data/spectral_resp/');
else
    rmpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\spectral_resp');
end

rmpath(['.' filesep 'Misc' filesep]); 