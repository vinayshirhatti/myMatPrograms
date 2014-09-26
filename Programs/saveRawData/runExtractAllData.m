% runExtractAllData
extractTheseIndices = 60;

monkeyName = 'rat'; gridType = 'EEG';

% [Vinay] - check if it is a Windows PC or an Linux machine and set the
% folderSourceString accordingly
if ispc
    folderSourceString = 'J:\'; % Vinay - changed source directory
else
    folderSourceString = '/media/store/';
end

[expDates,protocolNames,stimTypes] = allProtocolsTestMicroelectrode;

% stimEndTime = 3.0; interStim = 1.0; Fs = 2000; % Vinay: not required now

timeStartFromBaseLineList(1) = -0.55; deltaTList(1) = 1.024;
timeStartFromBaseLineList(2) = -1.148; deltaTList(2) = 2.048;
timeStartFromBaseLineList(3) = -1.096; deltaTList(3) = 4.096; % Vinay - for a longer stim time
timeStartFromBaseLineList(4) = -0.848; deltaTList(4) = 2.048; % 1000ms stim, 500ms interstim

electrodesToStore = [1:32 65:96];

for i=1:length(extractTheseIndices)
    expDate = expDates{extractTheseIndices(i)};
    protocolName = protocolNames{extractTheseIndices(i)};
    type = stimTypes{extractTheseIndices(i)};

    % Work on NEV data
    deltaT = deltaTList(type);
    timeStartFromBaseLine = timeStartFromBaseLineList(type);
    
    extractAllData(monkeyName,expDate,protocolName,folderSourceString,gridType,timeStartFromBaseLine,deltaT,electrodesToStore);
    %prepareDataForOrchestra(monkeyName,expDate,protocolName,electrodesToStore,folderSourceString,gridType);
    
    % extract LL Data if it exists
%     datFileName = [folderSourceString 'data\rawData\' monkeyName expDate '\' monkeyName expDate protocolName '.dat'];
%     if exist(datFileName,'file')
%         saveLLData(monkeyName,expDate,protocolName,folderSourceString,gridType,type);
%     end
end