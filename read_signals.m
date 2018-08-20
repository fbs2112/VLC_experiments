% %%
% %Sine B2BE
% clear;
% clc;
% close all;
% 
% % addpath('/home/felipe/Dropbox/VLC_exp/data/');
% 
% 
% if isunix
%     addpath('/home/felipe/Dropbox/VLC_exp/data/');
% else
%     addpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data');
% end
% 
% transmittedSignal = dlmread('Signal_Tx_SIN_B2BE.csv');
% figure;
% plot(transmittedSignal)
% title('Transmitted Signal')
% 
% 
% receivedSignal = dlmread('Signal_Rx_SIN_B2BE.csv');
% figure;
% plot(receivedSignal)
% 
% 
% %plot Correlation
% 
% figure;
% [corrTx,lags] = xcorr(transmittedSignal);
% 
% figure;
% plot(lags,corrTx)
% 
% [corrRx,lags] = xcorr(receivedSignal);
% figure;
% plot(lags,corrRx)
% 
% figure;
% [corrTx_Rx,lags] = xcorr(receivedSignal,transmittedSignal);
% 
% figure;
% plot(lags,corrTx_Rx)
% 
% [~,peaksOrdered] = sort(corrTx_Rx,'descend');
% 
% peak = lags(peaksOrdered(1));
% 
% RxSignalAligned = receivedSignal(peak+1:length(transmittedSignal) + peak);
% 
% figure;
% plot(RxSignalAligned);
% 
% varTx = var(transmittedSignal);
% 
% rescaledReceivedSignal = RxSignalAligned*sqrt(varTx/var(RxSignalAligned));
% 
% figure;
% plot(rescaledReceivedSignal)
% title('Received Signal')
% 
% 
% % rmpath('/home/felipe/Dropbox/VLC_exp/data/');
% 
% if isunix
%     rmpath('/home/felipe/Dropbox/VLC_exp/data/');
% else
%     rmpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data');
% end
% 
% %%
% %4-PAM B2BE
% clear;
% clc;
% close all;
% 
% % addpath('/home/felipe/Dropbox/VLC_exp/data/');
% if isunix
%     addpath('/home/felipe/Dropbox/VLC_exp/data/Aquisicao_24-04-18/Sinal com oversampling/');
% else
%     addpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\Aquisicao_24-04-18\Sinal com oversampling');
% end
% 
% 
% transmittedSignal = dlmread('4PAM_Tx_OVS_B2BE_aquisicoes_ 5_ (1).csv');
% % index = find(transmittedSignal);
% % transmittedSignal = transmittedSignal(index:end); 
% figure;
% scatterplot(transmittedSignal)
% title('Transmitted Signal')
% 
% 
% % receivedSignal = dlmread('Signal_Rx_4PAM_B2BE.csv');
% receivedSignal = dlmread('4PAM_Rx_OVS_B2BE_aquisicoes_5_ (1).csv');
% 
% figure;
% scatterplot(receivedSignal)
% N = 5;
% mu = 0.5;
% gamma = 1e-8;
% xAux = flipud(buffer(receivedSignal,N,N-1));
% 
% w = zeros(N,length(transmittedSignal),10);
% 
% for delayinSamples = 1:10
% % delayinSamples = 1;
% 
%     
%     for k = N + 10:length(transmittedSignal)
% 
%         x = receivedSignal(k:-1:k-N+1);
%     %     xTDLAux = zeros((N*N+N)/2,1);
%     %     
%     %     
%     %     for lIndex = 1:length(l1)
%     %         xTDLAux(lIndex,1) = x(l1(lIndex),1)*(x(l2(lIndex),1));
%     %     end
%     %     
%     %     
%     %     xConc = [x;xTDLAux];
% 
% 
%         d(k) = (transmittedSignal(-delayinSamples + k + 1));
% 
% 
%         e(k) = d(k) - w(:,k,delayinSamples)'*x;
% 
%         w(:,k+1,delayinSamples) = w(:,k,delayinSamples) + mu*x*((x'*x+gamma*eye(1))\eye(1))*conj(e(k));
% 
% 
%     end
% 
%     figure;
% 
%     plot(10*log10(abs(e).^2));
% end
% 
% 
% %plot Correlation
% 
% figure;
% [corrTx,lags] = xcorr(transmittedSignal);
% 
% figure;
% plot(lags,corrTx)
% 
% figure;
% [corrRx,lags] = xcorr(receivedSignal);
% 
% figure;
% plot(lags,corrRx)
% 
% [corrTx_Rx,lags] = gccPHATCorrelationOpt(receivedSignal,transmittedSignal);
% 
% figure;
% plot(lags,corrTx_Rx)
% 
% [~,peaksOrdered] = sort(corrTx_Rx,'descend');
% 
% peak = lags(peaksOrdered(1));
% 
% RxSignalAligned = receivedSignal(peak+1:length(transmittedSignal) + peak);
% 
% figure;
% 
% plot(RxSignalAligned);
% 
% varTx = var(transmittedSignal);
% 
% rescaledReceivedSignal = RxSignalAligned*sqrt(varTx/var(RxSignalAligned));
% 
% figure;
% plot(rescaledReceivedSignal)
% title('Received Signal')
% 
% decidedReceivedSignal = pamHardThreshold(rescaledReceivedSignal);
% 
% 
% if isunix
%     rmpath('/home/felipe/Dropbox/VLC_exp/data/Aquisicao_24-04-18/Sinal com oversampling/');
% else
%     rmpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\Aquisicao_24-04-18\Sinal com oversampling');
% end
% 
% %%
% %Sine LED
% clear;
% clc;
% close all;
% 
% % addpath('/home/felipe/Dropbox/VLC_exp/data/');
% 
% if isunix
%     addpath('/home/felipe/Dropbox/VLC_exp/data/');
% else
%     addpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data');
% end
% 
% transmittedSignal = dlmread('Signal_Tx_SIN_B2BO_450nm.csv');
% figure;
% plot(transmittedSignal)
% title('Transmitted Signal')
% 
% 
% receivedSignal = dlmread('Signal_Rx_SIN_B2BO_450nm.csv');
% figure;
% plot(receivedSignal)
% 
% 
% %plot Correlation
% 
% figure;
% [corrTx,lags] = xcorr(transmittedSignal);
% 
% figure;
% plot(lags,corrTx)
% 
% [corrRx,lags] = xcorr(receivedSignal);
% figure;
% plot(lags,corrRx)
% 
% figure;
% [corrTx_Rx,lags] = xcorr(receivedSignal,transmittedSignal);
% 
% figure;
% plot(lags,corrTx_Rx)
% 
% [~,peaksOrdered] = sort(corrTx_Rx,'descend');
% 
% peak = lags(peaksOrdered(1));
% 
% RxSignalAligned = receivedSignal(peak+1:length(transmittedSignal) + peak);
% 
% figure;
% plot(RxSignalAligned);
% 
% varTx = var(transmittedSignal);
% 
% rescaledReceivedSignal = RxSignalAligned*sqrt(varTx/var(RxSignalAligned));
% 
% figure;
% plot(rescaledReceivedSignal)
% title('Received Signal')
% 
% 
% 
% figure
% plot(transmittedSignal)
% hold on
% plot(rescaledReceivedSignal)
% H = legend('Transmitted signal','Received signal');
% set(H,'interpreter','latex');
% 
% if isunix
%     rmpath('/home/felipe/Dropbox/VLC_exp/data/');
% else
%     rmpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data');
% end
% 
% 
% freqz(rescaledReceivedSignal)
% %%
% %4-PAM B2BO
% clear;
% clc;
% close all;
% 
% addpath('/home/felipe/Dropbox/VLC_exp/data/');
% 
% 
% transmittedSignal = dlmread('Signal_Tx_4PAM_B2BO_450nm.csv');
% figure;
% scatterplot(transmittedSignal)
% title('Transmitted Signal')
% 
% 
% receivedSignal = dlmread('Signal_Rx_4PAM_B2BO_450nm.csv');
% figure;
% scatterplot(receivedSignal)
% N = 20;
% mu = 0.8;
% gamma = 1e-8;
% xAux = flipud(buffer(receivedSignal,N,N-1));
% 
% w = zeros(N,length(transmittedSignal),10);
% 
% for delayinSamples = 1:10
% % delayinSamples = 1;
% 
%     
%     for k = N + 10:length(transmittedSignal)
% 
%         x = xAux(k:-1:k-N+1);
%         xConc = x.';
%     %     xTDLAux = zeros((N*N+N)/2,1);
%     %     
%     %     
%     %     for lIndex = 1:length(l1)
%     %         xTDLAux(lIndex,1) = x(l1(lIndex),1)*(x(l2(lIndex),1));
%     %     end
%     %     
%     %     
%     %     xConc = [x;xTDLAux];
% 
% 
%         d(k) = (transmittedSignal(-delayinSamples + k + 1));
% 
% 
%         e(k) = d(k) - w(:,k,delayinSamples)'*xConc;
% 
%         w(:,k+1,delayinSamples) = w(:,k,delayinSamples) + mu*xConc*((xConc'*xConc+gamma*eye(1))\eye(1))*conj(e(k));
% 
% 
%     end
% 
%     figure;
% 
%     plot(10*log10(abs(e).^2));
% end
% 
% % 
% % %plot Correlation
% % 
% % figure;
% % [corrTx,lags] = xcorr(transmittedSignal);
% % 
% % figure;
% % plot(lags,corrTx)
% % 
% % figure;
% % [corrRx,lags] = xcorr(receivedSignal);
% % 
% % figure;
% % plot(lags,corrRx)
% % 
% % [corrTx_Rx,lags] = gccPHATCorrelationOpt(receivedSignal,transmittedSignal);
% % 
% % figure;
% % plot(lags,corrTx_Rx)
% % 
% % [~,peaksOrdered] = sort(corrTx_Rx,'descend');
% % 
% % peak = lags(peaksOrdered(1));
% % 
% % RxSignalAligned = receivedSignal(peak+1:length(transmittedSignal) + peak);
% % 
% % figure;
% % plot(RxSignalAligned);
% % 
% % varTx = var(transmittedSignal);
% % 
% % rescaledReceivedSignal = RxSignalAligned*sqrt(varTx/var(RxSignalAligned));
% % 
% % figure;
% % plot(rescaledReceivedSignal)
% % title('Received Signal')
% % 
% 
% rmpath('/home/felipe/Dropbox/VLC_exp/data/');
% 
% 
% %%
% 
% %Sine B2BE
% clear;
% clc;
% close all;
% 
% if isunix
%     addpath('/home/felipe/Dropbox/VLC_exp/data/PAM');
% else
%     addpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\data_aquisition_PAM');
% end
% load('Sinal_Tx_Rx_10symbs.mat');
% transmittedSignal = waveform;
% figure;
% stem(transmittedSignal)
% title('Transmitted Signal')
% 
% 
% receivedSignal = Signal_Rx;
% receivedSignal = receivedSignal - mean(receivedSignal);
% figure;
% stem(receivedSignal)
% 
% 
% %plot Correlation
% 
% % figure;
% [corrTx,lags] = xcorr(transmittedSignal);
% 
% % figure;
% % plot(lags,corrTx)
% 
% [corrRx,lags] = xcorr(receivedSignal);
% % figure;
% % plot(lags,corrRx)
% 
% % figure;
% [corrTx_Rx,lags] = xcorr(receivedSignal,transmittedSignal);
% 
% figure;
% plot(lags,corrTx_Rx)
% 
% [~,peaksOrdered] = sort(corrTx_Rx,'descend');
% 
% peak = lags(peaksOrdered(1));
% 
% RxSignalAligned = receivedSignal(peak+1:length(transmittedSignal) + peak);
% 
% figure;
% plot(RxSignalAligned);
% 
% 
% figure;
% plot(RxSignalAligned);
% title('Received Signal')
% 
% varTx = var(transmittedSignal);
% 
% rescaledReceivedSignal = RxSignalAligned*sqrt(varTx/var(RxSignalAligned));
% 
% figure;
% plot(rescaledReceivedSignal)
% title('Received Signal rescaled')
% 
% 
% % rmpath('/home/felipe/Dropbox/VLC_exp/data/');
% 
% if isunix
%     rmpath('/home/felipe/Dropbox/VLC_exp/data/');
% else
%     rmpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\data_aquisition_PAM');
% end


%%
%4-PAM B2BE Novo
clear;
clc;
close all;


if isunix
    addpath('/home/felipe/Dropbox/VLC_exp/data/2018-08-16');
else
    addpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\2018-08-16');
end

addpath(['.' filesep 'Algorithms' filesep]);

% load B2BO_IM_0.005_VDC_3.4V_Amp_340mVpp_B_2MHz_Beff_2.8MHz_Fs_50MSa.mat;
load B2BO_IM_0.005_VDC_3.5V_Amp_510mVpp_B_2MHz_Beff_2.8MHz_Fs_50MSa.mat;
% load B2BO_IM_0.005_VDC_3.6V_Amp_680mVpp_B_2MHz_Beff_2.8MHz_Fs_50MSa.mat;

% transmittedSignal = waveform;
reference = symbs_Tx;
receivedSignal = symbs_Rx;


if size(Signal_Rx,2) > 1
    sinal_recebido = Signal_Rx.';
end

figure;

stem(reference);
title('Transmitted signal' );
% index = find(transmittedSignal);
% transmittedSignal = transmittedSignal(index:end); 
scatterplot(reference)
title('Transmitted Signal Constellation')


% receivedSignal = dlmread('Signal_Rx_4PAM_B2BE.csv');
% receivedSignal = sinal_recebido;
% 
% receivedSignal_2 = receivedSignal.*sqrt(var(transmittedSignal)/var(receivedSignal));
% receivedSignal_2 = receivedSignal_2 - mean(receivedSignal_2);

figure;
stem(receivedSignal)
title('received')
scatterplot(receivedSignal)

N = 10;
mu = 0.3;
lambda = 0.98;
stopTrainingIt = 300;
gamma = 1e-8;
delayVector = N/2;
w = zeros(N,length(reference),length(delayVector));

for idx = 1:length(delayVector)
    e = zeros(length(reference),1);
    Sd = zeros(N,N,length(reference));
    Sd(:,:,N + 10) = eye(N)*(1-lambda)/1e-1; 
% delayinSamples = 1;
    for k = N + 10:length(reference)
        x = receivedSignal(k:-1:k-N+1);
        d(k) = (reference(-delayVector(idx) + k + 1));
        y(k,idx) = w(:,k,idx)'*x;
        [yy(k,idx),~] = pamHardThreshold2(y(k,idx));
        e(k) = conj(d(k) - y(k,idx));
        if k < stopTrainingIt
            [w(:,k+1,idx),Sd(:,:,k+1),~] = RLS_fun(w(:,k),Sd(:,:,k),x,lambda,e(k));
        else
            w(:,k+1,idx) = w(:,k,idx);
        end
    end
    figure;
    i = find(e,1);
    plot(10*log10(abs(e(i:end)).^2));
    ylim([-60 10]);
    
    [corrTx_Rx,lags] = xcorr(yy(:,idx),reference);
    [~,peaksOrdered] = sort(corrTx_Rx,'descend');
    peak = lags(peaksOrdered(1));
    error = yy(peak+1:end,idx) ~= reference(1:end - peak);
    symbError(idx) = sum(error(stopTrainingIt:end));
end

N = 8;
l1 = cell(length(N),1);
l2 = cell(length(N),1);
adapFiltLength = zeros(length(N),1);
for i = 1:length(N)
    auxMatrix = triu(ones(N(i)));
    [l1{i},l2{i}] = find(auxMatrix);
    adapFiltLength(i) = (N(i)^2+N(i))/2 + N(i);
end

mu = 0.4;
stopTrainingIt = 600;
gamma = 1e-8;
delayVector = round(N/2);
w = zeros(adapFiltLength,length(reference),length(delayVector));

for idx = 1:length(delayVector)
    e = zeros(length(reference),1);
    Sd = zeros(adapFiltLength,adapFiltLength,length(reference));
    Sd(:,:,N + 20) = eye(adapFiltLength)*(1-lambda)/1e-1; 
% delayinSamples = 1;
    for k = N + 20:length(reference)
        xLin = receivedSignal(k:-1:k-N+1);
        xNonLin = zeros(length(l1{1}),1);
        for lIndex = 1:length(l1{1})
            xNonLin(lIndex,1) = xLin(l1{1}(lIndex))*(xLin(l2{1}(lIndex)));
        end
        
        x = [xLin;xNonLin];
        d(k) = (reference(-delayVector(idx) + k + 1));
        y(k,idx) = w(:,k,idx)'*x;
        [yy(k,idx),~] = pamHardThreshold2(y(k,idx));
        e(k) = d(k) - y(k,idx);
        if k < stopTrainingIt
            [w(:,k+1,idx),Sd(:,:,k+1),~] = RLS_fun(w(:,k),Sd(:,:,k),x,lambda,e(k));
        else
            w(:,k+1,idx) = w(:,k,idx);
        end
    end
    figure;
    i = find(e,1);
    plot(10*log10(abs(e(i:end)).^2));
    ylim([-60 10]);
    
    [corrTx_Rx,lags] = xcorr(yy(:,idx),reference);
    [~,peaksOrdered] = sort(corrTx_Rx,'descend');
    peak = lags(peaksOrdered(1));
    error = yy(peak+1:end,idx) ~= reference(1:end - peak);
    symbErrorVolterra(idx) = sum(error(stopTrainingIt:end));
end

[~,receivedSignalSymb] = pamHardThreshold2(receivedSignal);

symbErroNoEq = sum(receivedSignalSymb ~= reference);

%plot Correlation
[corrTx_Rx,lags] = xcorr(receivedSignal,reference);
figure;
plot(lags,corrTx_Rx)
title('corr')
[~,peaksOrdered] = sort(corrTx_Rx,'descend');

% peak = lags(peaksOrdered(1));
% 
% RxSignalAligned = receivedSignal_2(peak:length(transmittedSignal) + peak-1);
% 
% figure;
% 
% stem(RxSignalAligned);
% title('Aligned')
% 
% varTx = var(transmittedSignal);
% 
% % RxSignalAligned = RxSignalAligned - mean(RxSignalAligned);
% 
% % rescaledReceivedSignal = RxSignalAligned*sqrt(varTx/var(RxSignalAligned ));
% 
% rescaledReceivedSignal = RxSignalAligned*3/max(RxSignalAligned);
% 
% figure;
% stem(rescaledReceivedSignal)
% title('Received Signal')
% 
% transmittedSignal2 = transmittedSignal*3/max(transmittedSignal);
% 
% decidedReceivedSignal = pamHardThreshold(rescaledReceivedSignal);
% decidedTransmittedSignal = pamHardThreshold(transmittedSignal2);
% 
% decDemodSignal = pamdemod(decidedReceivedSignal,4,0,'gray');
% binaryOutputData = de2bi(decDemodSignal,2);
% 
% 
% error = abs(transmittedSignal - rescaledReceivedSignal);
% 
% ser = error.'*error/length(error);
% 
% figure;
% stem(error)
% title('Error');
% 
% 
% rescaledReceivedSignal2(rescaledReceivedSignal > 0.5) = 1;
% rescaledReceivedSignal2(rescaledReceivedSignal <= 0.5) = 0;
% 
% error2 = abs(transmittedSignal - rescaledReceivedSignal2.');
% ser2 = error2.'*error2/length(error2);

if isunix
    rmpath('/home/felipe/Dropbox/VLC_exp/data/2018-08-16');
else
    rmpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\2018-08-16');
end
rmpath(['.' filesep 'Algorithms' filesep]);

