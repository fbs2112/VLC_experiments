clear;
clc;
close all;

addpath(['.' filesep 'VLC_param' filesep]);
addpath(['.' filesep 'VLC_param' filesep 'LED Parameters' filesep]);
addpath(['.' filesep 'Utils' filesep]);
addpath(['.' filesep 'Misc' filesep]);

load whiteLED_334-15_Param.mat;
load whiteLED_334-15.mat;
load VLC_param01.mat;

if isunix
    addpath('/home/felipe/Dropbox/VLC_exp/data');
    addpath('/home/felipe.silva/Dropbox/VLC_exp/data');
else
    addpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data');
end

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

load VDCxDCxP.mat;
P = P*1e-6;

receiverDiameter = 3.2e-3;

area = 4*pi*(receiverDiameter/2)^2;

plot(VDC, DC*1e3)
hold on

I_model = I_V_Fun(VDC,VT,nLED,ISat);

plot(VDC,I_model*1e3);

xlabel('Voltage [V]','interpreter','latex');
ylabel('Current [mA]','interpreter','latex');

H = legend('Experimental','Model');
set(H, 'interpreter','latex','location','best');
grid on;

formatFig( gcf ,['.' filesep 'figs' filesep 'current_voltage'], 'en' , figProp );

opticalPowerConstant = (maxLuminousIntensityLED/1000)/max(P/area);
    
P = P *opticalPowerConstant;
figure
plot(DC,P/area)

figure
plot(DC.* VDC,P/area)
hold on
% maxLuminousIntensityLED = 10000;
eletricalPowerOutput = DC.*VDC;
% k = 
% opticalPowerOutput = Poptical(ledLuminousEfficacy,eletricalPowerOutput,k,maxLuminousIntensityLED);
% plot(eletricalPowerOutput,opticalPowerOutput)


%--Find best k that fits the experimental data
counter = 1;
k = 1:0.01:200;

for index = 1:length(k)
    x(index,:) = Poptical(ledLuminousEfficacy,eletricalPowerOutput,k(index),maxLuminousIntensityLED);
    error(index) = sum(abs(x(index,:) - P/area).^2);
end

[~,minErrorIndex] = min(error);

opticalPowerOutput = Poptical(ledLuminousEfficacy,eletricalPowerOutput,k(minErrorIndex),maxLuminousIntensityLED);
plot(eletricalPowerOutput, opticalPowerOutput)
H = legend('Experimental','Model', 'Best model');
set(H, 'interpreter','latex','location','best');
xlim([0 round(max(DC.* VDC),2)]);



modFun = @(vmax, vdc) (vmax - vdc)/vdc; 
vmaxFun = @(MI, vdc) kron(vdc,MI);

voltageVector = [3.1];

 %quando a tens�o caiu 4 vezes do valor inicial;
% 
% vdc = data(:,2);
MI = [0.05 0.075 0.1];
% 
vmax = vmaxFun(MI,voltageVector)*2000; %Vpp [mV]

rmpath(['.' filesep 'VLC_param' filesep]);
rmpath(['.' filesep 'VLC_param' filesep 'LED Parameters' filesep]);
rmpath(['.' filesep 'Utils' filesep]);

if isunix
    rmpath('/home/felipe/Dropbox/VLC_exp/data');
    rmpath('/home/felipe.silva/Dropbox/VLC_exp/data');
    rmpath(['.' filesep 'Utils' filesep]);
else
    rmpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data');
end

rmpath(['.' filesep 'Utils' filesep]);
rmpath(['.' filesep 'Misc' filesep]);

