%==========================================================================
% Vinay Shirhatti
% Time Frequency Analysis
%
% Comparing the methods
%==========================================================================

N = 1024;
Fs = 2000;

load ('rB9.mat');
load ('rS8.mat');

[rB9stft rB9sp] = stftTFA(rB9,16,2);
[rS8stft rS8sp] = stftTFA(rS8,16,2);

freqlowlim = 3;
freqhighlim = 90;

kLowlim = freqlowlim*(N/Fs); 
kHighlim = freqhighlim*(N/Fs);

load ('rB9_gab.mat');
rB9_gD = gaborInfo{1}{1}.gaborData;
load ('rS8_gab.mat');
rS8_gD = gaborInfo{1}{1}.gaborData;



% % % 
% % % x=gaborInfo{1}{1}.gaborData
% % % size(x)
% % % plot(x(2,:))
% % % f=1000*x/512
% % % f=1000*x(2,:)/512
% % % plot(f)
% % % close
% % % plot(f)
% % % pos=intersect(find(f>=30),find(f<80))
% % % edit reconstructSignalFromAtomsMPP.m
% % % r=reconstructSignalFromAtomsMPP(x,512,1,pos)
% % % r=reconstructSignalFromAtomsMPP(x,511,1,pos)
% % % r=reconstructSignalFromAtomsMPP(x,1024,1,pos)
% % % plot(r)
% % % close
% % % plot(r)