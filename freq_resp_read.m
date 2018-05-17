%Sine B2BE
clear;
clc;
close all;

addpath(['.' filesep 'freq_resp' filesep]);

load('sine_B2BO.mat');
transmittedSignal = signal;
figure;
plot(transmittedSignal)
title('Transmitted Signal')


receivedSignal = sinal_recebido;
figure;
plot(receivedSignal)


%plot Correlation

% figure;
[corrTx,lags] = xcorr(transmittedSignal);

% figure;
% plot(lags,corrTx)

[corrRx,lags] = xcorr(receivedSignal);
% figure;
% plot(lags,corrRx)

% figure;
[corrTx_Rx,lags] = xcorr(receivedSignal,transmittedSignal);

figure;
plot(lags,corrTx_Rx)

[~,peaksOrdered] = sort(corrTx_Rx,'descend');

peak = lags(peaksOrdered(1));

RxSignalAligned = receivedSignal(peak+1:length(transmittedSignal) + peak);

figure;
plot(RxSignalAligned);


figure;
plot(RxSignalAligned);
title('Received Signal')

varTx = var(transmittedSignal);

rescaledReceivedSignal = RxSignalAligned*sqrt(varTx/var(RxSignalAligned));

figure;
plot(rescaledReceivedSignal)
title('Received Signal rescaled')


rmpath(['.' filesep 'freq_resp' filesep]);
