%==========================================================================
% Vinay Shirhatti, 31 January 2014
% Prepare data for job submission on Tyrone cluster, SERC
%==========================================================================

clear;clc;close all;

monkeyName = 'testMonk';
expDate = '310114';
protocolName = 'ANS_001';
channelNumbers = 1:10;
folderSourceString = '/media/Data/';
gridType = 'Microelectrode';

prepareDataForTyrone(monkeyName,expDate,protocolName,channelNumbers,folderSourceString,gridType);