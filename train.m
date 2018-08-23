clear;
clc;
% close all;

if isunix
    addpath('/home/felipe/Dropbox/VLC_exp/data/2018-08-16');
else
    addpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\2018-08-16');
end

addpath(['.' filesep 'Algorithms' filesep]);

% load B2BO_IM_0.005_VDC_3.4V_Amp_340mVpp_B_2MHz_Beff_2.8MHz_Fs_50MSa.mat;
% load B2BO_IM_0.005_VDC_3.5V_Amp_510mVpp_B_2MHz_Beff_2.8MHz_Fs_50MSa.mat;
% load B2BO_IM_0.005_VDC_3.6V_Amp_680mVpp_B_2MHz_Beff_2.8MHz_Fs_50MSa.mat;
load B2BO_IM_0.05_VDC_3.4V_Amp_340mVpp_B_2MHz_Beff_2.8MHz_Fs_200MSa.mat;  %MI = 0.05
% load B2BO_IM_0.075_VDC_3.4V_Amp_350mVpp_B_2MHz_Beff_2.8MHz_Fs_200MSa.mat; %MI = 0.075
% load B2BO_IM_0.1_VDC_3.4V_Amp_360mVpp_B_2MHz_Beff_2.8MHz_Fs_200MSa.mat; %MI = 0.1

% transmittedSignal = waveform;
reference = Data(:,2);
receivedSignal = Data(:,4);

N = 4;
mu = 0.3;
lambda = 0.98;
stopTrainingIt = 600;
gamma = 1e-8;
delayVector = round(N/2);
w = zeros(N,size(reference{1,1},1),size(reference,1));

for idx = 1:size(reference,1)
    e = zeros(size(reference{1,1},1),1);
    Sd = zeros(N,N,size(reference{1,1},1));
    Sd(:,:,N + 10) = eye(N)*(1-lambda)/1e-1; 
    for k = N + 10:size(reference{1,1},1)
        x = receivedSignal{idx,1}(k:-1:k-N+1);
        d(k) = (reference{idx,1}(-delayVector(1) + k + 1));
        y(k,idx) = w(:,k,idx)'*x;
        [yy(k,idx),~] = pamHardThreshold2(y(k,idx));
        e(k) = conj(d(k) - y(k,idx));
        if k < stopTrainingIt
            [w(:,k+1,idx),Sd(:,:,k+1),~] = RLS_fun(w(:,k),Sd(:,:,k),x,lambda,e(k));
%             [w(:,k+1,idx),~] = affine_projection_fun(w(:,k),x,mu,e(k));
        else
            w(:,k+1,idx) = w(:,k,idx);
        end
    end
    e2(:,idx) = abs(e).^2;
%     figure;
%     i = find(e,1);
%     plot(10*log10(abs(e(i:end)).^2));
%     ylim([-60 10]);
%     
%     [corrTx_Rx,lags] = xcorr(yy(:,idx),reference{1,1});
%     [~,peaksOrdered] = sort(corrTx_Rx,'descend');
%     peak = lags(peaksOrdered(1));
%     error = yy(peak+1:end,idx) ~= reference{1,1}(1:end - peak);
%     symbError(idx) = sum(error(stopTrainingIt:end));
%     symbErrorRate(idx) = symbError(idx)*100 / length(error(stopTrainingIt:end));
end
w2 = mean(w,3);
mse = mean(e2,2);
i = find(mse,1);
figure;
plot(10*log10(mse(i:stopTrainingIt)));
hold on


N = 4;
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
w = zeros(adapFiltLength,size(reference{1,1},1),size(reference,1));

for idx = 1:size(reference,1)
    e = zeros(size(reference{1,1},1),1);
    Sd = zeros(adapFiltLength,adapFiltLength,size(reference{1,1},1));
    Sd(:,:,N + 20) = eye(adapFiltLength)*(1-lambda)/1e-1; 
% delayinSamples = 1;
    for k = N + 20:size(reference{1,1},1)
        xLin = receivedSignal{idx,1}(k:-1:k-N+1);
        xNonLin = zeros(length(l1{1}),1);
        for lIndex = 1:length(l1{1})
            xNonLin(lIndex,1) = xLin(l1{1}(lIndex))*(xLin(l2{1}(lIndex)));
        end
        
        x = [xLin;xNonLin];
        d(k) = (reference{idx,1}(-delayVector(1) + k + 1));
        y(k,idx) = w(:,k,idx)'*x;
        [yy(k,idx),~] = pamHardThreshold2(y(k,idx));
        e(k) = d(k) - y(k,idx);
        if k < stopTrainingIt
            [w(:,k+1,idx),Sd(:,:,k+1),~] = RLS_fun(w(:,k),Sd(:,:,k),x,lambda,e(k));
%             [w(:,k+1,idx),~] = affine_projection_fun(w(:,k),x,mu,e(k));
        else
            w(:,k+1,idx) = w(:,k,idx);
        end
    end
    e2(:,idx) = abs(e).^2;
%     figure;
%     i = find(e,1);
%     plot(10*log10(abs(e(i:end)).^2));
%     ylim([-60 10]);
%     
%     [corrTx_Rx,lags] = xcorr(yy(:,idx),reference);
%     [~,peaksOrdered] = sort(corrTx_Rx,'descend');
%     peak = lags(peaksOrdered(1));
%     error = yy(peak+1:end,idx) ~= reference(1:end - peak);
%     symbErrorVolterra(idx) = sum(error(stopTrainingIt:end));
end

w2Volterra = mean(w,3);
mse = mean(e2,2);
i = find(mse,1);
plot(10*log10(mse(i:stopTrainingIt)));




mu = 0.4;
stopTrainingIt = 600;
gamma = 1e-8;
feedforwardLength = 8;
feedbackLength = 5;
delayVector = round(feedforwardLength/2);

adapFiltLength = feedforwardLength + feedbackLength;
w = zeros(adapFiltLength,size(reference{1,1},1),size(reference,1));

for idx = 1:size(reference,1)
    e = zeros(size(reference{1,1},1),1);
    Sd = zeros(adapFiltLength,adapFiltLength,size(reference{1,1},1));
    Sd(:,:,feedforwardLength + 20) = eye(adapFiltLength)*(1-lambda)/1e-1; 
% delayinSamples = 1;
    for k = feedforwardLength + 20:size(reference{1,1},1)
        x = receivedSignal{idx,1}(k:-1:k-feedforwardLength+1);
        y = (d(-delayVector(1) + k:-1:-delayVector(1) + k + 1 - feedbackLength)).';
        z = [x;y];
        d(k) = (reference{idx,1}(-delayVector(1) + k + 1));
        y(k,idx) = w(:,k,idx)'*z;
        [yy(k,idx),~] = pamHardThreshold2(y(k,idx));
        e(k) = conj(d(k) - y(k,idx));
        
        if k < stopTrainingIt
            [w(:,k+1,idx),Sd(:,:,k+1),~] = RLS_fun(w(:,k),Sd(:,:,k),z,lambda,e(k));
%             [w(:,k+1,idx),~] = affine_projection_fun(w(:,k),x,mu,e(k));
        else
            w(:,k+1,idx) = w(:,k,idx);
        end
    end
    e2(:,idx) = abs(e).^2;
%     figure;
%     i = find(e,1);
%     plot(10*log10(abs(e(i:end)).^2));
%     ylim([-60 10]);
%     
%     [corrTx_Rx,lags] = xcorr(yy(:,idx),reference);
%     [~,peaksOrdered] = sort(corrTx_Rx,'descend');
%     peak = lags(peaksOrdered(1));
%     error = yy(peak+1:end,idx) ~= reference(1:end - peak);
%     symbErrorVolterra(idx) = sum(error(stopTrainingIt:end));
end

w2DFE = mean(w,3);
mse = mean(e2,2);
i = find(mse,1);
plot(10*log10(mse(i:stopTrainingIt)));


mu = 0.4;
stopTrainingIt = 600;
gamma = 1e-8;
feedforwardLength = 8;
feedbackLength = 5;
delayVector = round(feedforwardLength/2);

l1 = cell(length(feedforwardLength),1);
l2 = cell(length(feedforwardLength),1);
adapFiltLength = zeros(length(feedforwardLength),1);
for i = 1:length(feedforwardLength)
    auxMatrix = triu(ones(feedforwardLength(i)));
    [l1{i},l2{i}] = find(auxMatrix);
    adapFiltLength(i) = (feedforwardLength(i)^2+feedforwardLength(i))/2 + feedforwardLength(i);
end
adapFiltLength = adapFiltLength + feedbackLength;

w = zeros(adapFiltLength,size(reference{1,1},1),size(reference,1));


for idx = 1:size(reference,1)
    e = zeros(size(reference{1,1},1),1);
    Sd = zeros(adapFiltLength,adapFiltLength,size(reference{1,1},1));
    Sd(:,:,feedforwardLength + 20) = eye(adapFiltLength)*(1-lambda)/1e-1; 
% delayinSamples = 1;
    for k = feedforwardLength + 20:size(reference{1,1},1)
        xLin = receivedSignal{idx,1}(k:-1:k-feedforwardLength+1);
        xNonLin = zeros(length(l1{1}),1);
        for lIndex = 1:length(l1{1})
            xNonLin(lIndex,1) = xLin(l1{1}(lIndex))*(xLin(l2{1}(lIndex)));
        end
        
        x = [xLin;xNonLin];
        y = (d(-delayVector(1) + k:-1:-delayVector(1) + k + 1 - feedbackLength)).';
        z = [x;y];
        d(k) = (reference{idx,1}(-delayVector(1) + k + 1));
        y(k,idx) = w(:,k,idx)'*z;
        [yy(k,idx),~] = pamHardThreshold2(y(k,idx));
        e(k) = conj(d(k) - y(k,idx));
        
        if k < stopTrainingIt
            [w(:,k+1,idx),Sd(:,:,k+1),~] = RLS_fun(w(:,k),Sd(:,:,k),z,lambda,e(k));
%             [w(:,k+1,idx),~] = affine_projection_fun(w(:,k),x,mu,e(k));
        else
            w(:,k+1,idx) = w(:,k,idx);
        end
    end
    e2(:,idx) = abs(e).^2;
%     figure;
%     i = find(e,1);
%     plot(10*log10(abs(e(i:end)).^2));
%     ylim([-60 10]);
%     
%     [corrTx_Rx,lags] = xcorr(yy(:,idx),reference);
%     [~,peaksOrdered] = sort(corrTx_Rx,'descend');
%     peak = lags(peaksOrdered(1));
%     error = yy(peak+1:end,idx) ~= reference(1:end - peak);
%     symbErrorVolterra(idx) = sum(error(stopTrainingIt:end));
end

w2VolterraDFE = mean(w,3);
mse = mean(e2,2);
i = find(mse,1);
plot(10*log10(mse(i:stopTrainingIt)));

% figure
% stem(w2(:,end))
% figure
% stem(w2Volterra(:,end))
% figure
% stem(w2DFE(:,end))
% figure
% stem(w2VolterraDFE(:,end))

if isunix
    rmpath('/home/felipe/Dropbox/VLC_exp/data/2018-08-16');
else
    rmpath('C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\2018-08-16');
end
rmpath(['.' filesep 'Algorithms' filesep]);

