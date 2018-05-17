clear;
clc;
close all;

% opticalPower = 1e-6*[6.98;29.68;58.82;71.31;87.24;99.8;110;119;126.6;135.1;142.5;149.6;155.5;158.2;160;160;162;159;150;147.8;145;145;141;132;122;111;93;89];
%     
% % current = 1e-3*[1 5:5:125];
% 
% voltage = [2.8;3;3.3;3.6;3.8;4.1;4.3;4.5;4.8;5;5.3;5.6;5.8;6.1;6.3;6.6;6.8;6.9;7;7.2;7.5;7.7;7.9;8.2;8.4;8.7;8.9;9];

% plot(current(1:end-2), voltage)
% 
% figure
% plot(current(1:end-2),opticalPower)
% 
% figure
% plot(current(1:end-2).*voltage.',opticalPower)


% data = [1	2.8	10
% 5	3	43.2
% 10	3.4	78
% 15	3.7	107
% 20	3.9	132.9
% 25	4.2	155.4
% 30	4.5	175.4
% 35	4.7	191.5
% 40	5	205.7
% 45	5.2	215.8
% 50	5.5	224.7
% 55	5.7	231.1
% 60	5.9	237.8
% 65	6.2	240.5];

load DCeVDC.mat;

opticalPower = [0.05
0.554
1.47
3
5.2
7.5
11.05
15.68
21.81
29.05
36.87
45.2
53.78
62.53
71.33
80.3
89.48
96.55
97.08
97.02
96.97
96.84
96.82
96.75
96.74
96.6
96.55
96.63
96.59
96.53] * 1e-6;
receiverDiameter = 3.2e-3;

area = pi*(receiverDiameter/2)^2;

plot(VDC, DC)

aux = find(round(VDC,2)==2.35);
figure
plot(DC(aux:length(opticalPower)+aux-1),opticalPower/area)

figure
plot(DC(aux:length(opticalPower)+aux-1).* VDC(aux:length(opticalPower)+aux-1),opticalPower/area)



    


load whiteLED_334-15_Param.mat;
load whiteLED_334-15.mat;
load VLC_param01.mat;

eletricalPowerOutput = 0:0.01:0.3;
opticalPowerOutput = Poptical(ledLuminousEfficacy,eletricalPowerOutput,kNonLinearity,maxLuminousIntensityLED);

figure

plot(eletricalPowerOutput, opticalPowerOutput)

modFun = @(vmax, vdc) (vmax - vdc)/vdc; 
vmaxFun = @(MI, vdc) (vdc.*MI);

voltageVector = [2.8 3];

 %quando a tens�o caiu 4 vezes do valor inicial;
% 
% vdc = data(:,2);
MI = [0.1];
% 
vmax = vmaxFun(MI,voltageVector)*2000; %Vpp [mV]




