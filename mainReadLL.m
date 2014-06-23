% monkeyName = 'test';
% expDate = '170214';
% protocolName = 'GRF_001';
% folderSourceString = '/media/Data/';
% gridType = 'Microelectrode';
% type = 1;

% monkeyName = 'eyeData';
% expDate = '280514';
% protocolName = 'GRF_001';
% folderSourceString = '/media/Data/';
% gridType = 'Microelectrode';
% type = 1;

% monkeyName = 'monk';
% expDate = '290514';
% protocolName = 'CRS_002';
% folderSourceString = '/media/Data/';
% gridType = 'Microelectrode';
% type = 1;

monkeyName = 'test';
expDate = '100614';
protocolName = 'GRF_002';
folderSourceString = '/media/Data/';
gridType = 'Microelectrode';
type = 1;


% saveLLData(monkeyName,expDate,protocolName,folderSourceString,gridType,type);


timeStartFromBaseLine= -0.55;
deltaT = 1.024;
electrodesToStore = [147,148,149];

extractAllData(monkeyName,expDate,protocolName,folderSourceString,timeStartFromBaseLine,deltaT,electrodesToStore,gridType)
