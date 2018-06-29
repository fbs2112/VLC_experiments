function [ADC] = connect_ADC(addr)
ADC = instrfind('Type', 'visa-tcpip', 'RsrcName', addr, 'Tag', '');
if isempty(ADC)
    ADC = visa('NI',addr);
    ADC.OutputBufferSize = 999999999;
    ADC.InputBufferSize = 99999999;
    ADC.timeout=10;
    fopen(ADC);
end

