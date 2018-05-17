close all
[DAC] = connect_DAC('USB0::0x0699::0x0345::C022370::0::INSTR');
% [ADC] = connect_ADC('USB0::0x0957::0x1799::MY52163581::0::INSTR');
IP='192.168.0.85';
ADC=connect_ADC(['TCPIP0::' IP '::inst0::INSTR']);

k=100;
f=1e7;
fs=k*f; %DAC;
Fs_ADC=fs*2; % ADC

t = 0:1/fs:(1/f)*1;

signal = sin(2*pi*f*t);

Vpp = max(signal)*2;
waveform=rectpulse(signal,k);
time_window=1/f*10;
% points=k*length(signal);
points = time_window*fs;
send_to_AWG(DAC,waveform,f,Vpp);
sinal_recebido = get_from_scope_test(ADC,points,time_window);

figure;
subplot(2,1,1)
plot(waveform)
title('Sinal Transmitido com oversampling')
subplot(2,1,2)
plot(sinal_recebido)

fileName = 'sine_B2BO';

figure
c = fft(signal,1000*2^nextpow2(length(signal)));
freq  = fs*(0:(length(c)/2))/length(c);
P = c / length(c);
plot(freq,(abs(P(1:end/2+1))).');     

save([fileName '.mat'], 'sinal_recebido', 'signal');
% save(['C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\2017-05-09' fileName '.mat'], 'sinal_recebido');


