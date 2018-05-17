function Sinal = get_from_scope_Agilent(ADC,points,time_window)
% Reset the instrument and autoscale and stop
fprintf(ADC,':CHANnel2:PROBe 1');
fprintf(ADC,':CHANnel1:PROBe 1');

fprintf(ADC,'*RST; :AUTOSCALE'); 
% fprintf(ADC,':STOP');
% Specify data from Channel 1

fprintf(ADC,':WAVEFORM:SOURCE CHAN1'); 
fprintf(ADC,':TIMEBASE:MODE MAIN');
fprintf(ADC,':ACQuire:TYPE NORMal');
% fprintf(ADC,':ACQUIRE:COUNT 1'); % Tomar wfm uma vez
fprintf(ADC,':WAV:POINTS:MODE RAW');
fprintf(ADC, [':TIMebase:SCALe ' num2str(time_window/10)]);
fprintf(ADC,[':WAV:POINTS ' num2str(points)]);
fprintf(ADC,':DIGITIZE CHAN1');
fprintf(ADC,':WAV:DATA?');
fprintf(ADC,':CHANnel1:COUPling AC')
waveform.RawData = binblockread(ADC,'uint16');
fread(ADC,1);

operationComplete = str2double(query(ADC,'*OPC?'));
while ~operationComplete
    operationComplete = str2double(query(ADC,'*OPC?'));
end
fprintf(ADC,':WAVeform:FORMat WORD');
fprintf(ADC,':WAVEFORM:BYTEORDER LSBFirst');
preambleBlock = query(ADC,':WAVEFORM:PREAMBLE?');
fprintf(ADC,':WAV:DATA?');
waveform.RawData = binblockread(ADC,'uint16'); fread(ADC,1);

% Maximum value storable in a INT16
maxVal = 2^16; 

%  split the preambleBlock into individual pieces of info
preambleBlock = regexp(preambleBlock,',','split');

% store all this information into a waveform structure for later use
waveform.Format = str2double(preambleBlock{1});     % This should be 1, since we're specifying INT16 output
waveform.Type = str2double(preambleBlock{2});
waveform.Points = str2double(preambleBlock{3});
waveform.Count = str2double(preambleBlock{4});      % This is always 1
waveform.XIncrement = str2double(preambleBlock{5}); % in seconds
waveform.XOrigin = str2double(preambleBlock{6});    % in seconds
waveform.XReference = str2double(preambleBlock{7});
waveform.YIncrement = str2double(preambleBlock{8}); % V
waveform.YOrigin = str2double(preambleBlock{9});
waveform.YReference = str2double(preambleBlock{10});
waveform.VoltsPerDiv = (maxVal * waveform.YIncrement / 8);      % V
waveform.Offset = ((maxVal/2 - waveform.YReference) * waveform.YIncrement + waveform.YOrigin);         % V
waveform.SecPerDiv = waveform.Points * waveform.XIncrement/10 ; % seconds
waveform.Delay = ((waveform.Points/2 - waveform.XReference) * waveform.XIncrement + waveform.XOrigin); % seconds

% Generate X & Y Data
waveform.XData = (waveform.XIncrement.*(1:length(waveform.RawData))) - waveform.XIncrement;
waveform.YData = (waveform.YIncrement.*(waveform.RawData - waveform.YReference)) + waveform.YOrigin; 
Sinal=waveform.YData;
T_Signal_Rx=waveform.XData;
fprintf(ADC,':CHANnel1:PROBe 1');
fprintf(ADC,':CHANnel2:PROBe 1');

    