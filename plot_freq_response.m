clear;
clc;
close all

if isunix
    addpath('/home/felipe/Dropbox/VLC_exp/data/freq_resp');
    addpath('/home/felipe.silva/Dropbox/VLC_exp/data/freq_resp');
else
    addpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\freq_resp');
end
%MI_1 = 0.05
%MI_2 = 0.075
%MI_3 = 0.1

%VDC_1 = 3.4
%VDC_2 = 3.5
%VDC_3 = 3.6
%VDC_4 = 3.1

data_MI_1_VDC_1 = load('MI0.05_340mvpp_3.4V.mat');
data_MI_1_VDC_2 = load('MI0.05_510mvpp_3.5V.mat');
data_MI_1_VDC_3 = load('MI0.05_680mvpp_3.6V.mat');
data_MI_1_VDC_4 = load('MI0.05_465mvpp_3.1V.mat');


data_MI_2_VDC_1 = load('MI0.075_350mvpp_3.4V.mat');
data_MI_2_VDC_2 = load('MI0.075_525mvpp_3.5V.mat');
data_MI_2_VDC_3 = load('MI0.075_700mvpp_3.6V.mat');
data_MI_2_VDC_4 = load('MI0.075_620mvpp_3.1V.mat');


data_MI_3_VDC_1 = load('MI0.1_360mvpp_3.4V.mat');
data_MI_3_VDC_2 = load('MI0.1_540mvpp_3.5V.mat');
data_MI_3_VDC_3 = load('MI0.1_720mvpp_3.6V.mat');
data_MI_3_VDC_4 = load('MI0.1_310mvpp_3.1V.mat');



w = linspace(200e3,10e6,1e3);

w1 = 15.5e6/(2*pi);
w2 = 3.26e6/(2*pi);
w3 = 10.86e6/(2*pi);
wc = 1e6/(2*pi);

Gb = @(wAux) exp(-wAux/w1);

figure
plot(w,(Gb(w)));
hold on;

wAux2 = w(w<wc);

Gw = exp(-wAux2/w2);

GwLength = length(Gw);

wAux3 = w(w>wc);

Gw2 = exp(-wc/w2)*exp(wc/w3)*exp(-wAux3/w3);

maxLength = max(GwLength,length(Gw2));

Gw = [Gw zeros(1,length(w) - GwLength)];

Gw(1,GwLength+1:end) = Gw2;

plot(w,(Gw),'r');

% refline(0,-3)
% 
% ylim([-10 0])
% 
% xlim([0 20e6])

grid on;

xlabel('Frequency (MHz)','interpreter','latex');
ylabel('Amplitude (dB)','interpreter','latex');

H = legend('Blue Channel','White Channel');
set(H,'interpreter','latex')

f = data_MI_1_VDC_1.freq;


data_MI_1_VDC_1_aux = interp1(f, data_MI_1_VDC_1.Amplitude, w);%interpolated frequency response
data_MI_1_VDC_2_aux = interp1(f, data_MI_1_VDC_2.Amplitude, w);%interpolated frequency response
data_MI_1_VDC_3_aux = interp1(f, data_MI_1_VDC_3.Amplitude, w);%interpolated frequency response
data_MI_1_VDC_4_aux = interp1(f, data_MI_1_VDC_4.Amplitude, w);%interpolated frequency response

data_MI_2_VDC_1_aux = interp1(f, data_MI_2_VDC_1.Amplitude, w);%interpolated frequency response
data_MI_2_VDC_2_aux = interp1(f, data_MI_2_VDC_2.Amplitude, w);%interpolated frequency response
data_MI_2_VDC_3_aux = interp1(f, data_MI_2_VDC_3.Amplitude, w);%interpolated frequency response
data_MI_2_VDC_4_aux = interp1(f, data_MI_2_VDC_4.Amplitude, w);%interpolated frequency response

data_MI_3_VDC_1_aux = interp1(f, data_MI_3_VDC_1.Amplitude, w);%interpolated frequency response
data_MI_3_VDC_2_aux = interp1(f, data_MI_3_VDC_2.Amplitude, w);%interpolated frequency response
data_MI_3_VDC_3_aux = interp1(f, data_MI_3_VDC_3.Amplitude, w);%interpolated frequency response
data_MI_3_VDC_4_aux = interp1(f, data_MI_3_VDC_4.Amplitude, w);%interpolated frequency response


figure;
plot(f,(data_MI_1_VDC_4.Amplitude/max(data_MI_1_VDC_4.Amplitude)))
hold on
plot(f,(data_MI_2_VDC_4.Amplitude/max(data_MI_2_VDC_4.Amplitude)))
plot(f,(data_MI_3_VDC_4.Amplitude/max(data_MI_3_VDC_4.Amplitude)))
plot(w,(Gw));

H = legend('MI = 0.05, VDC = 3.1','MI = 0.075, VDC = 3.1',...
    'MI = 0.1, VDC = 3.1','Model');

set(H,'interpreter','latex')
xlim([min(w) max(w)]);



figure
plot(f,(data_MI_1_VDC_1.Amplitude/max(data_MI_1_VDC_1.Amplitude)))
hold on
plot(f,(data_MI_1_VDC_2.Amplitude/max(data_MI_1_VDC_2.Amplitude)))
plot(f,(data_MI_1_VDC_3.Amplitude/max(data_MI_1_VDC_3.Amplitude)))
plot(f,(data_MI_2_VDC_1.Amplitude/max(data_MI_2_VDC_1.Amplitude)))
plot(f,(data_MI_2_VDC_2.Amplitude/max(data_MI_2_VDC_2.Amplitude)))
plot(f,(data_MI_2_VDC_3.Amplitude/max(data_MI_2_VDC_3.Amplitude)))
plot(f,(data_MI_3_VDC_1.Amplitude/max(data_MI_3_VDC_1.Amplitude)))
plot(f,(data_MI_3_VDC_2.Amplitude/max(data_MI_3_VDC_2.Amplitude)))
plot(f,(data_MI_3_VDC_3.Amplitude/max(data_MI_3_VDC_3.Amplitude)))
plot(w,(Gw));

H = legend('MI = 0.05, VDC = 3.4','MI = 0.05, VDC = 3.5',...
    'MI = 0.05, VDC = 3.6','MI = 0.075, VDC = 3.4','MI = 0.075, VDC = 3.5'...
    ,'MI = 0.075, VDC = 3.6','MI = 0.1, VDC = 3.4','MI = 0.1, VDC = 3.5',...
    'MI = 0.1, VDC = 3.6','Model');

set(H,'interpreter','latex')
xlim([min(w) max(w)]);


figure
plot(w,(data_MI_1_VDC_1_aux/max(data_MI_1_VDC_1_aux)))
hold on
plot(w,(data_MI_1_VDC_2_aux/max(data_MI_1_VDC_2_aux)))
plot(w,(data_MI_1_VDC_3_aux/max(data_MI_1_VDC_3_aux)))
plot(w,(data_MI_2_VDC_1_aux/max(data_MI_2_VDC_1_aux)))
plot(w,(data_MI_2_VDC_2_aux/max(data_MI_2_VDC_2_aux)))
plot(w,(data_MI_2_VDC_3_aux/max(data_MI_2_VDC_3_aux)))
plot(w,(data_MI_3_VDC_1_aux/max(data_MI_3_VDC_1_aux)))
plot(w,(data_MI_3_VDC_2_aux/max(data_MI_3_VDC_2_aux)))
plot(w,(data_MI_3_VDC_3_aux/max(data_MI_3_VDC_3_aux)))
plot(w,(Gw));

H = legend('MI = 0.05, VDC = 3.4','MI = 0.05, VDC = 3.5',...
    'MI = 0.05, VDC = 3.6','MI = 0.075, VDC = 3.4','MI = 0.075, VDC = 3.5'...
    ,'MI = 0.075, VDC = 3.6','MI = 0.1, VDC = 3.4','MI = 0.1, VDC = 3.5',...
    'MI = 0.1, VDC = 3.6','Model');

set(H,'interpreter','latex')
xlim([min(w) max(w)]);

figure
plot(w,10*log10(data_MI_1_VDC_1_aux/max(data_MI_1_VDC_1_aux)))
hold on
plot(w,10*log10(data_MI_1_VDC_2_aux/max(data_MI_1_VDC_2_aux)))
plot(w,10*log10(data_MI_1_VDC_3_aux/max(data_MI_1_VDC_3_aux)))
plot(w,10*log10(data_MI_2_VDC_1_aux/max(data_MI_2_VDC_1_aux)))
plot(w,10*log10(data_MI_2_VDC_2_aux/max(data_MI_2_VDC_2_aux)))
plot(w,10*log10(data_MI_2_VDC_3_aux/max(data_MI_2_VDC_3_aux)))
plot(w,10*log10(data_MI_3_VDC_1_aux/max(data_MI_3_VDC_1_aux)))
plot(w,10*log10(data_MI_3_VDC_2_aux/max(data_MI_3_VDC_2_aux)))
plot(w,10*log10(data_MI_3_VDC_3_aux/max(data_MI_3_VDC_3_aux)))
plot(w,10*log10(Gw));

H = legend('MI = 0.05, VDC = 3.4','MI = 0.05, VDC = 3.5',...
    'MI = 0.05, VDC = 3.6','MI = 0.075, VDC = 3.4','MI = 0.075, VDC = 3.5'...
    ,'MI = 0.075, VDC = 3.6','MI = 0.1, VDC = 3.4','MI = 0.1, VDC = 3.5',...
    'MI = 0.1, VDC = 3.6','Model');

set(H,'interpreter','latex','Location','Best')
xlim([min(w) max(w)]);

if isunix
    rmpath('/home/felipe/Dropbox/VLC_exp/data/freq_resp');
    rmpath('/home/felipe.silva/Dropbox/VLC_exp/data/freq_resp');
else
    rmpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\freq_resp');
end




