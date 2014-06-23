% This function is used to extract the following
% 1. Digital data
% 2. LFP data
% 3. Spike data
% 4. Segments

function extractAllDataCRS(monkeyName,expDate,protocolName,folderSourceString,timeStartFromBaseLine,deltaT,electrodesToStore,gridType)

% if ~exist('folderSourceString','var')   folderSourceString ='E:\';
% end
if ~exist('folderSourceString','var')   folderSourceString ='/media/Data/';      end % Vinay - for linux
if ~exist('timeStartFromBaseLine','var') timeStartFromBaseLine= -0.55;  end
if ~exist('deltaT','var')                deltaT = 1.024;                end

folderSourceString = appendIfNotPresent(folderSourceString,'\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% isWin=1; % variable for makeDirectory
isWin=0; % variable for makeDirectory. Because this is not Windows

% fileName = [monkeyName expDate protocolName '.nev'];
% folderName0 = [folderSourceString 'data\' monkeyName '\'];
% makeDirectory(folderName0);
% folderName0 = [folderName0 gridType '\'];
% makeDirectory(folderName0);
% folderName1 = [folderName0 expDate '\'];
% makeDirectory(folderName1);
% folderName = [folderName1 protocolName '\'];
% makeDirectory(folderName);
% 
% folderIn = [folderSourceString 'data\rawData\' monkeyName expDate '\'];
% folderExtract = [folderName 'extractedData\'];
% folderSegment = [folderName 'segmentedData'];

% Vinay - changed all the above steps for linux
fileName = [monkeyName expDate protocolName '.nev'];
folderName0 = [folderSourceString 'data/' monkeyName '/'];
makeDirectory(folderName0);
folderName0 = [folderName0 gridType '/'];
makeDirectory(folderName0);
folderName1 = [folderName0 expDate '/'];
makeDirectory(folderName1);
folderName = [folderName1 protocolName '/'];
makeDirectory(folderName);

folderIn = [folderSourceString 'data/rawData/' monkeyName expDate '/'];
folderExtract = [folderName 'extractedData/'];
folderSegment = [folderName 'segmentedData'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the NEV file

% Load the appropriate DLL
% dllName = 'C:\Supratim\programs\Matlab-Import-Filter_2_5\RequiredResources\nsNEVLibrary.dll';
dllName = '/media/Data/RequiredResources/nsNEVLibrary.dll'; % Vinay - for linux
[nsresult] = ns_SetLibrary(dllName);
if (nsresult ~= 0)      error('DLL was not found!');                    end

% Load data file and display some info about the file open data file
[nsresult, hFile] = ns_OpenFile([folderIn fileName]);
if (nsresult ~= 0)      error('Data file did not open!');               end

% Get file information
[nsresult, fileInfo] = ns_GetFileInfo(hFile);
% Gives you entityCount, timeStampResolution and timeSpan
if (nsresult ~= 0)      error('Data file information did not load!');   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Digital Codes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, entityInfo] = ns_GetEntityInfo(hFile,1:fileInfo.EntityCount);
eventList = find([entityInfo.EntityType] == 1);

digitalEventName = 'digin';
for i =1:length(eventList)
    if strcmp(digitalEventName,[entityInfo(eventList(i)).EntityLabel]);
        digID=i;
        break;
    end
end

if ~exist('digID','var')    
    error(['No event named ' digitalEventName]);
else
    numDigitalEvents=entityInfo(eventList(digID)).ItemCount;
end
disp(['Number of digital events: ' num2str(numDigitalEvents)]);

% Get the digital events
[~,digitalTimeStamps,digitalEvents] = ns_GetEventData(hFile,eventList(digID),1:numDigitalEvents);

% the codes all start with a leading 1, which means that they are greater than hex2dec(8000) = 32768.
modifiedDigitalEvents = digitalEvents(digitalEvents>32768) - 32768;
allCodesInDec = unique(modifiedDigitalEvents);

disp(['Number of distinct codes: ' num2str(length(allCodesInDec))]);

count=1;
for i=1:length(allCodesInDec)
    if isempty(digitalCodeDictionary(hex2str(dec2hex(allCodesInDec(i)))))
        error('Identified digital code');
    else
        identifiedDigitalCodes(count) = allCodesInDec(i);
        count=count+1;
    end
end

numDigitalCodes = length(identifiedDigitalCodes);
disp(['Number of distinct codes identified: ' num2str(numDigitalCodes)]);

for i=1:numDigitalCodes
    digitalCodeInfo(i).codeNumber = identifiedDigitalCodes(i); %#ok<*AGROW>
    digitalCodeInfo(i).codeName = hex2str(dec2hex(identifiedDigitalCodes(i)));
    clear codePos
    codePos = find(identifiedDigitalCodes(i) == digitalEvents-32768);
    digitalCodeInfo(i).time = digitalTimeStamps(codePos);
    digitalCodeInfo(i).value = digitalEvents(codePos+1);
end

% Write the digitalCodes
makeDirectory(folderExtract,isWin);

save([folderExtract 'NEVFileInfo.mat'], 'fileInfo', 'entityInfo');
save([folderExtract 'digitalEvents.mat'],'digitalCodeInfo','digitalTimeStamps','digitalEvents');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Get Stimulus results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strncmp(protocolName,'SRC',3)
    readDigitalCodesSRC(folderExtract); % writes stimResults and trialResults
    [goodStimNums,goodStimTimes] = getGoodStimNumsSRC(folderExtract); % Good stimuli
    getDisplayCombinationsSRC(folderExtract,goodStimNums); % writes parameterCombinations
    save([folderExtract 'goodStimNums.mat'],'goodStimNums');
    
elseif strncmp(protocolName,'GRF',3) || strncmp(protocolName,'DRF',3)
    readDigitalCodesGRF(folderExtract); % writes stimResults and trialResults
    [goodStimNums,goodStimTimes] = getGoodStimNumsGRF(folderExtract); % Good stimuli
    getDisplayCombinationsGRF(folderExtract,goodStimNums);
    save([folderExtract 'goodStimNums.mat'],'goodStimNums');
    
elseif strncmp(protocolName,'GRA',3) || strncmp(protocolName,'ANS',3)
    readDigitalCodesGRF(folderExtract); % writes stimResults and trialResults
    [goodStimNums,goodStimTimes] = getGoodStimNumsGRF(folderExtract); % Good stimuli
    getDisplayCombinationsGRA(folderExtract,goodStimNums);
    save([folderExtract 'goodStimNums.mat'],'goodStimNums');

elseif strncmp(protocolName,'CRS',3)
    readDigitalCodesCRS(folderExtract); % writes stimResults and trialResults
    [goodStimNums,goodStimTimes] = getGoodStimNumsGRF(folderExtract); % Good stimuli
    getDisplayCombinationsGRA(folderExtract,goodStimNums);
    save([folderExtract 'goodStimNums.mat'],'goodStimNums'); 
end

Fs=2000;
analogChannelsToStore = electrodesToStore;
neuralChannelsToStore = analogChannelsToStore;
getLFP=1;getSpikes=1;
getLFPandSpikes(fileName,analogChannelsToStore,folderIn,folderSegment, ...
    goodStimTimes,timeStartFromBaseLine,deltaT,Fs,hFile,neuralChannelsToStore,getLFP,getSpikes);
end