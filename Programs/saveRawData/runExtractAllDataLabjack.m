% runExtractAllData
extractTheseIndices = 8;

monkeyName = 'test'; gridType = 'Microelectrode'; folderSourceString = 'F:\';
[expDates,protocolNames,stimTypes] = allProtocolsTestMicroelectrode;

timeStartFromBaseLineList(1) = -0.55; deltaTList(1) = 1.024;
timeStartFromBaseLineList(2) = -1.148; deltaTList(2) = 2.048;

electrodesToStore = 1:2;

for i=1:length(extractTheseIndices)
    expDate = expDates{extractTheseIndices(i)};
    protocolName = protocolNames{extractTheseIndices(i)};
    type = stimTypes{extractTheseIndices(i)};

    deltaT = deltaTList(type);
    timeStartFromBaseLine = timeStartFromBaseLineList(type);
    
    extractAllDataLabjack(monkeyName,expDate,protocolName,folderSourceString,gridType,timeStartFromBaseLine,deltaT,electrodesToStore);
    %prepareDataForOrchestra(monkeyName,expDate,protocolName,electrodesToStore,folderSourceString,gridType);
    
    % extract LL Data if it exists
    datFileName = [folderSourceString 'data\rawData\' monkeyName expDate '\' monkeyName expDate protocolName '.dat'];
    if exist(datFileName,'file')
        saveLLData(monkeyName,expDate,protocolName,folderSourceString,gridType,type);
    end
end