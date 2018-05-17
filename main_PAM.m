clearvars
% f2=exist('x2','var')
[DAC] = connect_DAC('USB0::0x0699::0x0345::C022370::0::INSTR');
% [ADC] = connect_ADC('USB0::0x0957::0x1799::MY52163581::0::INSTR');


IP='192.168.0.85';
ADC=connect_ADC(['TCPIP0::' IP '::inst0::INSTR']);

simbolos=100;
sf=10;
flag=zeros(sf,1);
k=100;
f=10e6;
fs=k*f; %DAC;
Fs_ADC=fs*2; % ADC

signal = randi([0 1], simbolos, 1);
% signal= real(pammod(signal, 4, 0, 'gray'))/3;


stem(signal)


Vpp = max(signal)*2;
% signal = randi([0 3], simbolos, 1);
% signal_PAM= real(pammod(signal, 4, 0, 'gray'))/2;
% sinal_flag=([flag;signal_PAM;flag]);
waveform=rectpulse(signal,k);
time_window=1/f*20;
% points=k*length(signal);
points = time_window*fs;
send_to_AWG(DAC,waveform,f,Vpp);
sinal_recebido = get_from_scope_test(ADC,points,time_window);
% sinal_recebido(sinal_recebido>=-2.5 & sinal_recebido<=2.5)=0;

rescaledReceivedSignal = sinal_recebido*sqrt(var(signal)/var(sinal_recebido));


figure;
subplot(2,1,1)
stem(waveform)
title('Sinal Transmitido com oversampling')
subplot(2,1,2)
stem(rescaledReceivedSignal)


figure

c = fft(signal,10*2^nextpow2(length(signal)));

freq  = fs*(0:(length(c)/2))/length(c);

P = c / length(c);

plot(freq,20*log10(abs(P(1:end/2+1))).');     
% fileName = 'sine_B2BO';
save(['PAM' '.mat'], 'sinal_recebido', 'signal');
% save(['C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\2017-05-09' fileName '.mat'], 'sinal_recebido');


