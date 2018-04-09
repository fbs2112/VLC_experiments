clear;
clc;
close all;


addpath(['.' filesep 'VLC_param' filesep]); 
addpath(['.' filesep 'VLC_param' filesep 'Utils' filesep]); 
addpath(['.' filesep 'VLC_param' filesep 'LED Parameters' filesep]); 


load VLC_param01.mat;

SNR = pow2db(30);
numberOfBits = 2;
numberOfSymbols = 2^numberOfBits;
numberOfSymbolsInBlock = 100;
signalPower = 1;
modulationIndexVector = 0.05;

maxVoltage = VDC*(1+modulationIndexVector);

input = randi([0,numberOfSymbols-1], numberOfSymbolsInBlock, 1);

pilot = pammod(input, numberOfSymbols, 0, 'gray');

pilot2 = pilot.*sqrt(signalPower/var(pilot));

xAux = VLC_channel(pilot2, modulationIndexVector, maxVoltage, VDC, SNR);

rmpath(['.' filesep 'VLC_param' filesep]); 
rmpath(['.' filesep 'VLC_param' filesep 'Utils' filesep]); 
rmpath(['.' filesep 'VLC_param' filesep 'LED Parameters' filesep]); 