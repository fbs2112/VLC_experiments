function [DAC] = connect_DAC(addr)
DAC = instrfind('Type', 'visa-usb', 'RsrcName', addr , 'Tag', '');
if isempty(DAC)
    DAC = visa('NI', addr);
    DAC.OutputBufferSize = 99999999;
    DAC.InputBufferSize = 99999999;
    fopen(DAC);
end
