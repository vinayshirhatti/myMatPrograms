%==========================================================================
% Vinay Shirhatti, 31 January 2014
% Prepare data for job submission on Tyrone cluster, SERC
%==========================================================================

clear;clc;close all;

monkeyName = 'alpa';
expDate = '130814';
protocolName = 'CRS_004';
channelNumbers = [65 66 69 70 73 74 77 78];
folderSourceString = '/media/store/';
gridType = 'Microelectrode';

prepareDataForTyrone(monkeyName,expDate,protocolName,channelNumbers,folderSourceString,gridType);
