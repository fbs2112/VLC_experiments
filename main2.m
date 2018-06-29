teste = exist('m','var');
if teste ==0
    myFgen = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x0699::0x0342::C020390::0::INSTR', 'Tag', '');    %  Create the VISA-TCPIP object if it does not exist
    myFgen = visa('NI', 'USB0::0x0699::0x0342::C020390::0::INSTR');
    myFgen.OutputBufferSize = 512000; % Aumentar o buffer de saída de acordo com o sinal
    fopen(myFgen); % Ativar conexão com AFG
%     IP = '192.168.0.85';
%     port = 80;
%     OSC = instrfind('Type', 'tcpip', 'RemoteHost', '12', 'RemotePort', port, 'Tag', ''); % localizar instrumento
%     visaAddress = (['TCPIP0::' IP '::inst0::INSTR']); % endereço para VISA usado para o instrumento
%     OSC = visa('ni', visaAddress);
%     OSC.OutputBufferSize = 99999999;
%     OSC.InputBufferSize = 99999999;
%     OSC.timeout = 100;
%     fopen(OSC); % abrir conexão
end;
signal = randi([0 3], 500, 1);
waveform = real(pammod(signal, 4, 0, 'gray'))/2;
m=1;
%% Iniciar aquisição
f=1e6;
ovs=250;
Fs=ovs*f;
ncyc=1;
t=0:1/Fs:1/f*ncyc;
Vnc=150e-3; % vppmax pra 650nm 2.9V 30mA
% waveform=sin(2*pi*f*t*ncyc)*150e-3;
Vpp=3.3;
send_to_AWG(waveform,Fs,f,Vpp);
o=500; % Fator utilizado para aumentar a janela de aqs
time_window=o*max(t);
Fs2=Fs;
n2=Fs2*time_window;
Sinal_recebido= get_from_scope(Fs2,time_window);
t2=0:1/Fs2:time_window;
t2=t2(1:end-1);
% t2=t2(1:end-1);
figure;
plot(t2,Sinal_recebido);

save('Signal_Tx_4PAM_B2BO_450nm.csv','waveform','-ascii')
save('Signal_Rx_4PAM_B2BO_450nm.csv','Sinal_recebido','-ascii')