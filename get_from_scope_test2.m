function sinal = get_from_scope_test2(OSC,Fs,time_window)
% time_window=2.6e-3;
fwrite(OSC,'AUToscale'); %Autoset
pause(6);
% z1 = max(csvread('controle_scalini.csv'));
% z2 = max(csvread('controle_scalini2.csv'));
% if z1==z2
%     AGC(time_window,Fs)
% end
% fprintf(OSC, 'CHANnel1:WAVeform1:ARIThmetics AVERage');

fwrite(OSC, 'CHANnel 1 :STATe 1'); %ativar canal1
% fprintf(OSC, 'CHANnel1:COUPling AC');
fwrite(OSC, 'ACQuire:MODE RTIMe'); %modo de interpolação de sinal
fwrite(OSC, 'ACQuire:INTerpolate SINX'); %Tipo de intepolação
% fwrite(OSC, 'ACQuire:POInts:AUTO RECLength'); %Tipo de intepolação
fwrite(OSC, 'ACQuire:POInts:AUTO RES'); %Tipo de intepolação

fwrite(OSC, ['TIMebase:RANGe ' num2str(time_window)])  % Definir a janela temporal do osc, a duração do sinal
% fwrite(OSC,['ACQuire:POInts:VALue '  num2str(points)]) 
fwrite(OSC, ['ACQuire:SRATe ' num2str(Fs)]); % Definir a taxa de amostragem
% fprintf(OSC, 'CHANnel1:COUPling AC');
fprintf(OSC, 'CHANnel1:COUPling DC');
fwrite(OSC, 'EXPort:WAVeform:SOURce C1W1'); % Ativar painel para sinal simples, sem fft apenas o sinal.
fwrite(OSC, 'EXPort:WAVeform:SCOPe WFM');
fwrite(OSC, 'EXPort:WAVeform:FASTexport 1')
fwrite(OSC, 'EXPort:WAVeform:RAW OFF'); % Desativar comunicação RAW
fwrite(OSC, 'EXPort:WAVeform:INCXvalues OFF'); % Ativar no save o registro temporal do sinal
fwrite(OSC, 'EXPort:WAVeform:DLOGging OFF'); % Desativar datalogging
fwrite(OSC, 'RUNContinuous'); % Trigger continuo
pause(1);
fprintf(OSC, 'REFCurve1:STATe ON'); %REFCurve<m>:STATe <State> Habilitar reference waveform m =1 state = on
fprintf(OSC, 'REFCurve1:SOURce 1'); %REFCurve<m>:SOURce source 1 ch1wfm1
fprintf(OSC, 'REFCurve1:UPDate'); % Jogar para a memória do osciloscópio s
fprintf(OSC, 'FORMat:DATA ASCii'); %Definir formato do sinal para Ascii
fprintf(OSC, 'CHAN1:WAV1:DATA?'); % escrever comando para OSC
% pause(30);
tic
sinal = fscanf(OSC); % ler o comando
toc
% pause(20);
tic
sinal = str2num(sinal);
toc

flushinput(OSC);
flushoutput(OSC);
% pause(60);
% fwrite(OSC,'*RST');
fprintf(OSC, 'CHANnel1:COUPling AC');
fprintf(OSC, 'REFCurve1:CLEar'); %REFCurve<m>:STATe <State> Habilitar reference waveform m =1 state = on
fprintf(OSC, 'REFCurve1:STATe OFF'); %REFCurve<m>:STATe <State> Habilitar reference waveform m =1 state = on
fwrite(OSC,'*RST');
fprintf(OSC, 'CHANnel1:COUPling AC');
% pause(20);
% while length(sinal)<500000
%     get_from_scope('192.168.0.85',80,Fs,time_window)
    