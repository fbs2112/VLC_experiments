function [ADC] = connect_IP(addr)
    ADC = instrfind('Type', 'visa-tcpip', 'RsrcName', addr, 'Tag', '');; % localizar instrumento
%     obj1 = instrfind('Type', 'visa-tcpip', 'RsrcName', 'TCPIP0::192.168.0.85::inst0::INSTR', 'Tag', '');

if isempty(ADC)
    ADC = visa('NI', addr);
    ADC.OutputBufferSize = 99999999;
    ADC.InputBufferSize = 9999999;
    fopen(ADC);
end
