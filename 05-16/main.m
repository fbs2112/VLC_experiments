[DPS]=connect_USB('USB0::0x1AB1::0x0E11::DP8B160800050::0::INSTR');
[DAC]= connect_USB('USB0::0x0699::0x0345::C022370::0::INSTR');
[ADC]= connect_IP('TCPIP0::192.168.0.85::inst0::INSTR');
Vdc=2.8;
fprintf(DPS,':OUTP CH1,ON');
fprintf(DPS,[':APPL CH1,' num2str(Vdc)]);
fprintf(DPS,':OUTP CH1,OFF');
