% prepareDataForTyrone is a modified version of prepareDataForOrchestra 
% which is a generic program that converts the data in a
% format that can be run on the Tyrone cluster, SERC

function prepareDataForTyrone(monkeyName,expDate,protocolName,channelNumbers,folderSourceString,gridType)

% 
% folderSourceString = appendIfNotPresent(folderSourceString,'\');
% folderNameMain = [folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName '\'];
% 
% % Input folder
% inputFolder  = [folderNameMain 'segmentedData\LFP\'];
% 
% % Output folder
% outputFolder = [folderNameMain 'mpAnalysis\'];
% makeDirectory(outputFolder);
% 
% % This is where the data will be visible on orchestra
% linuxFolder  = ['/groups/maunsell/data/' monkeyName '/' gridType '/' expDate '/' protocolName '/mpAnalysis/'];
% 

% Vinay - modified lines ahead

folderSourceString = appendIfNotPresent(folderSourceString,'/');
folderNameMain = [folderSourceString 'data/' monkeyName '/' gridType '/' expDate '/' protocolName '/'];

% Input folder
inputFolder  = [folderNameMain 'segmentedData/LFP/'];

% Output folder
outputFolder = [folderNameMain 'mpAnalysis/'];
makeDirectory(outputFolder);

% This is where the data will be visible on Tyrone
linuxFolder  = ['/localscratch/apmsvd/MP/data/' monkeyName '/' gridType '/' expDate '/' protocolName '/mpAnalysis/']; % Replace 'apmsvd' by the relevant 'username' 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Segmenting data...');
for i=1:length(channelNumbers)
    disp(channelNumbers(i));
    clear analogData analogInfo analogMatFile
    
    analogMatFile = [inputFolder 'elec' num2str(channelNumbers(i)) '.mat'];
    
    if ~exist(analogMatFile,'file')
        disp([analogMatFile ' does not exist']);
    else
        load(analogMatFile);
        X = analogData';

        if ~exist('Fs','var')
            Fs = analogInfo.SampleRate; % Initialize
        else
            if (Fs ~= analogInfo.SampleRate)
                error('Sampling rates not the same!!');
            end
        end
        
        if ~exist('L','var')
            L = size(X,1); % Initialize
        else
            if (size(X,1) ~= L)
                error('signal lengths not the same!!');
            end
        end

%         tag = ['elec' num2str(channelNumbers(i)) '\'];
        tag = ['elec' num2str(channelNumbers(i)) '/']; % Vinay - for linux
        range = [1 L];
        importData(X,outputFolder,tag,range,Fs);

        Numb_points = L;
        Max_iterations = 100; % [Vinay] - changed from 500 to 100
        prepareMPForOrchestra(outputFolder,tag,Numb_points,Max_iterations,linuxFolder);
    end
end

% write script file
writeScriptFileTyrone(monkeyName,expDate,protocolName,channelNumbers,folderSourceString,gridType);
end