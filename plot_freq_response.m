clear;
clc;
close all

%MI_1 = 0.05
%MI_2 = 0.075
%MI_3 = 0.1

%VDC_1 = 2.8
%VDC_2 = 3


data_MI_1_VDC_1 = [159
145
130
115
100
86
74
64
56
47
42
37
33
29
26
23
21
19
18
16
15]*1e-3;

data_MI_1_VDC_2 = [165
159
145
130
116
102
89
79
70
62
55
49
45
40
37
34
31
29
26
25
22]*1e-3;

data_MI_2_VDC_1 = [221
203
185
165
145
127
110
95
82
72
64
57
50
45
40
36
33
30
28
25
23]*1e-3;

data_MI_2_VDC_2 = [264
241
218
195
171
150
132
116
104
93
83
74
67
61
56
50
46
43
40
37
34]*1e-3;

data_MI_3_VDC_1 = [274
252
230
208
185
160
140
123
107
94
83
73
65
59
53
48
44
40
37
33
31]*1e-3;

data_MI_3_VDC_2 = [347
316
286
256
227
200
176
155
137
123
110
98
89
81
74
67
61
56
52
48
45]*1e-3;

f = [0.2
0.5
1
1.5
2
2.5
3
3.5
4
4.5
5
5.5
6
6.5
7
7.5
8
8.5
9
9.5
10]*1e6;



w = linspace(200e3,10e6,1e3);

w1 = 15.5e6/(2*pi);
w2 = 3.26e6/(2*pi);
w3 = 10.86e6/(2*pi);
wc = 1e6/(2*pi);

Gb = @(wAux) exp(-wAux/w1);

figure
plot(w,(Gb(w)));
hold on;

wAux2 = w(w<wc);

Gw = exp(-wAux2/w2);

GwLength = length(Gw);

wAux3 = w(w>wc);

Gw2 = exp(-wc/w2)*exp(wc/w3)*exp(-wAux3/w3);

maxLength = max(GwLength,length(Gw2));

Gw = [Gw zeros(1,length(w) - GwLength)];

Gw(1,GwLength+1:end) = Gw2;

plot(w,(Gw),'r');

% refline(0,-3)
% 
% ylim([-10 0])
% 
% xlim([0 20e6])

grid on;

xlabel('Frequency (MHz)','interpreter','latex');
ylabel('Amplitude (dB)','interpreter','latex');

H = legend('Blue Channel','White Channel');
set(H,'interpreter','latex')



data_MI_1_VDC_1_aux = interp1(f, data_MI_1_VDC_1, w);%interpolated frequency response
data_MI_1_VDC_2_aux = interp1(f, data_MI_1_VDC_2, w);%interpolated frequency response
data_MI_2_VDC_1_aux = interp1(f, data_MI_2_VDC_1, w);%interpolated frequency response
data_MI_2_VDC_2_aux = interp1(f, data_MI_2_VDC_2, w);%interpolated frequency response
data_MI_3_VDC_1_aux = interp1(f, data_MI_3_VDC_1, w);%interpolated frequency response
data_MI_3_VDC_2_aux = interp1(f, data_MI_3_VDC_2, w);%interpolated frequency response

% figure
% plot(f,(data_MI_1_VDC_1/max(data_MI_1_VDC_1)))
% hold on
% plot(w,(data_MI_1_VDC_1_aux/max(data_MI_1_VDC_1_aux)))

figure
plot(w,(data_MI_1_VDC_1_aux/max(data_MI_1_VDC_1_aux)))
hold on
plot(w,(data_MI_1_VDC_2_aux/max(data_MI_1_VDC_2_aux)))
plot(w,(data_MI_2_VDC_1_aux/max(data_MI_2_VDC_1_aux)))
plot(w,(data_MI_2_VDC_2_aux/max(data_MI_2_VDC_2_aux)))
plot(w,(data_MI_3_VDC_1_aux/max(data_MI_3_VDC_1_aux)))
plot(w,(data_MI_3_VDC_2_aux/max(data_MI_3_VDC_2_aux)))
plot(w,(Gw));

H = legend('MI = 0.05, VDC = 2.8','MI = 0.05, VDC = 3',...
    'MI = 0.075, VDC = 2.8','MI = 0.075, VDC = 3','MI = 0.1, VDC = 2.8','MI = 0.1, VDC = 3','Model');

set(H,'interpreter','latex')
xlim([min(w) max(w)]);

figure
plot(w,10*log10(data_MI_1_VDC_1_aux/max(data_MI_1_VDC_1_aux)))
hold on
plot(w,10*log10(data_MI_1_VDC_2_aux/max(data_MI_1_VDC_2_aux)))
plot(w,10*log10(data_MI_2_VDC_1_aux/max(data_MI_2_VDC_1_aux)))
plot(w,10*log10(data_MI_2_VDC_2_aux/max(data_MI_2_VDC_2_aux)))
plot(w,10*log10(data_MI_3_VDC_1_aux/max(data_MI_3_VDC_1_aux)))
plot(w,10*log10(data_MI_3_VDC_2_aux/max(data_MI_3_VDC_2_aux)))
plot(w,10*log10(Gw));

H = legend('MI = 0.05, VDC = 2.8','MI = 0.05, VDC = 3',...
    'MI = 0.075, VDC = 2.8','MI = 0.075, VDC = 3','MI = 0.1, VDC = 2.8','MI = 0.1, VDC = 3','Model');

set(H,'interpreter','latex','Location','Best')
xlim([min(w) max(w)]);






