clearvars
close all

[DPS]=connect_USB('USB0::0x1AB1::0x0E11::DP8B160800050::0::INSTR');
[DAC]= connect_USB('USB0::0x0699::0x0345::C022370::0::INSTR');
% [ADC]= connect_IP('TCPIP0::192.168.0.85::inst0::INSTR');
[ADC]= connect_USB('USB0::0x0957::0x1799::MY52163579::0::INSTR');


Vdc=2.8;
fprintf(DPS,':OUTP CH1,ON');
fprintf(DPS,[':APPL CH1,' num2str(Vdc)]);

simbolos=10;
sf=10;
flag=zeros(sf,1);
k=50;
f=20e3;
Fs=k*f; %DAC;
Fs_ADC=Fs*2; % ADC

% waveform = randi([0 1], simbolos, 1)*0.14;
% signal= real(pammod(signal, 4, 0, 'gray'))/3;


% stem(waveform)

% t = 0:1/Fs:(1/f)*1-1/Fs;

% signal = sin(2*pi*f*t)*0.14;


signal = randi([0 3], simbolos, 1);
signal= real(pammod(signal, 4, 0, 'gray'))*0.01;
Vpp = max(signal)*2;

% sinal_flag=([flag;signal_PAM;flag]);
% waveform=rectpulse(signal,k);
time_window=1/f*2;
points=2*k;
% points = time_window*Fs;
send_to_AWG(DAC,signal,f,Vpp);
% sinal_recebido = get_from_scope_test2(ADC,Fs,time_window);
sinal_recebido= get_from_scope_Agilent(ADC,points,time_window);
% sinal_recebido(sinal_recebido>=-2.5 & sinal_recebido<=2.5)=0;

rescaledReceivedSignal = sinal_recebido*sqrt(var(signal)/var(sinal_recebido));


figure;
subplot(2,1,1)
stem(signal)
title('Sinal Transmitido com oversampling')
subplot(2,1,2)
stem(rescaledReceivedSignal)


figure

c = fft(rescaledReceivedSignal,100*2^nextpow2(length(signal)));

freq  = Fs*(0:(length(c)/2))/length(c);

P = c / length(c);

plot(freq,(abs(P(1:end/2+1))).');     
% fileName = 'sine_B2BO';
save(['..' filesep 'PAM' '.mat'], 'sinal_recebido', 'signal');
% save(['C:\Users\Felipe Barboza\Dropbox\VLC_exp\data\2017-05-09' fileName '.mat'], 'sinal_recebido');


% fprintf(DPS,':OUTP CH1,OFF');
