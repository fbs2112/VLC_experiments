function send_to_AWG(myFgen,waveform,f,Vpp)

waveform=-waveform;
% Vppwfm= abs(max(waveform))+abs(min(waveform)); %Vpp do sinal criado no MATLAB
kgen=2^14-1;
% QF = Vppwfm/Vpp;
lwfm= length(waveform);
waveform = -waveform; % teoricamente não será mais necessário
fprintf(myFgen, '*RST');
fprintf(myFgen, '*CLS;'); 
myFgen.ByteOrder = 'littleEndian';
% f = Fs/length(waveform); % Frequência de sinal de saída no AWG  [x]
% timeVec = 0:1/Fs:1/f; % Vetor temporal formado conforme a taxa de amostragem
% timeVec = timeVec(1:end-1); % Eliminar ultimo bit do sinal para o comprimento dos dados no afg ficar correto
timeVec=length(waveform);
waveformLength = length(waveform); % obter comprimento do sinal
Vppwfm=max(waveform)+abs(min(waveform));
% VppDAC=10;
CF=Vppwfm/Vpp;
ymax=8191;
ymin=-8191;
quant_factor=(2^14-1)/Vpp;
waveform=round(waveform*quant_factor);
waveform(waveform>ymax) = ymax;
waveform(waveform<ymin) = ymin;
waveform=waveform+ymax;
waveform=uint16(waveform);
binblock = zeros(2 * waveformLength, 1); % Matriz de zeros que irá ter o dobro do tamanho do sinal,
binblock(2:2:end) = bitand(waveform, max(waveform)); % No vetor binblock, valores pares da matriz recebem valores, entre valores do sinal e outro valor
binblock(1:2:end) = bitshift(waveform, -8); % Valores ímpares do vetor são multiplicados por uma potência de 2, para ser mantido o sinal com segurança multiplicar por 2e-8
binblock = binblock'; % Vetor transposto
% Build binary block header
bytes = num2str(length(binblock));
header = ['#' num2str(length(bytes)) bytes]; % Converter o comprimento do binblock pra um valor de possível decod pro afg
% Resets the contents of edit memory and define the length of signal
fprintf(myFgen, ['DATA:DEF EMEM, ' num2str(length(timeVec)) ';']); %comprimento do sinal jogado pra memoria do AFG
% conversao do header e do binblock para a linguagem do estado islamico pro
% AFG
% Transfer the custom waveform from MATLAB to edit memory of instrument
fwrite(myFgen, [':TRACE EMEM, ' header binblock ';'], 'uint8'); % fazer download do sinal pra memoria de edicao do AFG do header binclock
% Associate the waveform in edit memory to channel 1
fprintf(myFgen, 'SOUR1:FUNC EMEM'); %  Ativar arb (da memoria de edicao no AFG)
fprintf(myFgen, ':OUTP1 ON');
f = num2str(f);
fremand = (['SOUR1:FREQ:FIXed ', f ,'Hz']);
fprintf(myFgen,fremand);
% Vpp=800e-3;
fprintf(myFgen,['SOURce1:VOLTage:LEVel:IMMediate:AMPlitude ' num2str(Vpp) 'Vpp']);
