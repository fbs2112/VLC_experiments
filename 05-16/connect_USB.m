function [DAC] = connect_USB(addr)
DAC = instrfind('Type', 'visa-usb', 'RsrcName', addr , 'Tag', '');
if isempty(DAC)
    DAC = visa('NI', addr);
    DAC.OutputBufferSize = 99999999;
    DAC.InputBufferSize = 9999999;
    fopen(DAC);
end
