clearvars
[DAC] = connect_DAC('USB0::0x0699::0x0345::C022370::0::INSTR');
% [ADC] = connect_ADC('USB0::0x0957::0x1799::MY52163581::0::INSTR');
IP='192.168.0.85'
ADC=connect_ADC(['TCPIP0::' IP '::inst0::INSTR'])
simbolos=100;
sf=10;
flag=zeros(sf,1);
k=2;
Vpp=3.4;
f=1e3;
Fs_DAC=k*f/1e6; %DAC;
Fs_ADC=Fs_DAC*2; % ADC
signal = randi([0 3], simbolos, 1);
signal_PAM= real(pammod(signal, 4, 0, 'gray'))/2;
sinal_flag=([flag;signal_PAM;flag]);
waveform=rectpulse(sinal_flag,k);
time_window=1/f*k;
points=k*length(waveform);
send_to_AWG(DAC,waveform,f,Vpp);
sinal_recebido = get_from_scope(ADC,points,time_window);
sinal_recebido(sinal_recebido>=-2.5 & sinal_recebido<=2.5)=0;

figure;
subplot(2,1,1)
plot(waveform)
title('Sinal Transmitido com oversampling')
subplot(2,1,2)
plot(sinal_recebido)
