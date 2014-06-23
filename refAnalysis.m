%==========================================================================
% Vinay Shirhatti, 14 February 2014
% Understanding the effects of referencing schemes
%==========================================================================

clear all; close all; clc;

Fs = 2000;
dt = 1/Fs;

N = 1000;

t = dt*(1:N);

freqVals = Fs/N*(0:(N/2)-1); % Span of frequencies to be considered
sigLen = length(freqVals);

%%

st1 = sin(2*pi*6*t) + rand(1,length(t)).*sin(2*pi*26*t);
k = 0.2;
st2 = k.*st1 + 0.314.*sin(2*pi*48*t);
st3 = st1 - st2;


%%
plot(st1,'r');hold on;
plot(st2, 'k');
plot(st3, 'c');

%%

sw1 = fft(st1); sw2 = fft(st2); sw3 = fft(st3);

p1 = angle(st1); p2 = angle(st2); p3 = angle(st3);

%%
plot(freqVals,abs(sw1(1:sigLen)),'r');hold on;
plot(freqVals,abs(sw2(1:sigLen)), 'k');
plot(freqVals,abs(sw3(1:sigLen)), 'c');

%%
figure;
plot(p1,'r');hold on;
plot(p2, 'k');
plot(p3, 'c');

pdiff = p3 - p1;
figure;plot(freqVals,pdiff(1:sigLen));