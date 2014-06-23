%==========================
% temp
%=========================

% for i=1:length(ori)
% for j = 1:length(cont)
% for k = 1:length(rad)
% ocr{i}{j}{k} = find((LL.orientationDeg1 == ori(i)) & (LL.contrastPC1 == cont(j)) & (LL.radiusDeg2 == rad(k)));
% end
% end
% end

clear all; clc; close all;

% 
% monkeyName = 'test';
% expDate = '180214';
% protocolName = 'Fixate_002';
% folderSourceString = '/media/Data/';
% gridType = 'Microelectrode';
% type = '1';

monkeyName = 'eyeData';
expDate = '100414';
protocolName = 'ED_003';
folderSourceString = '/media/Data/';
gridType = 'Microelectrode';
type = '1';


% monkeyName = 'eyeData';
% expDate = '290514';
% protocolName = 'GRF_002';
% folderSourceString = '/media/Data/';
% gridType = 'Microelectrode';
% type = '1';

% monkeyName = 'test';
% expDate = '100614';
% protocolName = 'GRF_002';
% folderSourceString = '/media/Data/';
% gridType = 'Microelectrode';
% type = 1;

FsEye=200;
%FsEye=120;
eyeRangeMS{1} = [-320 320]; timePeriodMS{1} = [-500   500]; maxStimPos{1} = 20;
eyeRangeMS{2} = [-480 800]; timePeriodMS{2} = [-1000 1000]; maxStimPos{2} = 10;

% folderName    = [folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName '\'];
folderName    = [folderSourceString 'data/' monkeyName '/' gridType '/' expDate '/' protocolName '/']; % [Vinay] for linux
folderExtract = [folderName 'extractedData'];


%LL = getStimResultsLLGRF(monkeyName,expDate,protocolName,folderSourceString);
%     save([folderExtract '\LL.mat'],'LL');
%save([folderExtract '/LL.mat'],'LL'); % Vinay - for linux
    
eyeRangeMS = [-250 250];
Fs = 200;
%Fs = 120;

eyeRangePos = eyeRangeMS*Fs/1000;

% if ~exist('folderSourceString','var')   folderSourceString ='F:\';       end
if ~exist('folderSourceString','var')   folderSourceString ='/media/';       end % Vinay - for linux

% monkeyName = removeIfPresent(monkeyName,'\');
% expDate    = removeIfPresent(expDate,'\');
% protocolName = removeIfPresent(protocolName,'\');
% folderSourceString = appendIfNotPresent(folderSourceString,'\');
% Vinay - for linux
monkeyName = removeIfPresent(monkeyName,'/');
expDate    = removeIfPresent(expDate,'/');
protocolName = removeIfPresent(protocolName,'/');
folderSourceString = appendIfNotPresent(folderSourceString,'/');

% datFileName = [folderSourceString 'data\rawData\' monkeyName expDate '\' monkeyName expDate protocolName '.dat'];
datFileName = [folderSourceString 'data/rawData/' monkeyName expDate '/' monkeyName expDate protocolName '.dat']; % Vinay - for linux

% Get Lablib data
header = readLLFile('i',datFileName);

% Stimulus properties
numTrials = header.numberOfTrials;
stimNumber=1;
correctIndex=1;
trialEndIndex=1;

if (strcmp(protocolName,'GRF_001') || strcmp(protocolName,'GRF_002'))
    for i=1:numTrials
    %disp(i);
    clear trials
    trials = readLLFile('t',i);
    
    if isfield(trials,'trialEnd')
        allTrials.trialEnded(i) = 1;
        allTrials.catchTrials(trialEndIndex) = trials.trial.data.catchTrial;
        allTrials.instructTrials(trialEndIndex) = trials.trial.data.instructTrial;
        allTrials.trialCertify(trialEndIndex) = trials.trialCertify.data;
        allTrials.targetPosAllTrials(trialEndIndex) = trials.trial.data.targetIndex+1;
        allTrials.eotCodes(trialEndIndex) = trials.trialEnd.data;
        
        allTrials.fixWindowSize(trialEndIndex) = trials.fixWindowData.data.windowDeg.size.width;
        allTrials.respWindowSize(trialEndIndex) = trials.responseWindowData.data.windowDeg.size.width;
        allTrials.certifiedNonInstruction(trialEndIndex) = (allTrials.instructTrials(trialEndIndex)==0)*(allTrials.trialCertify(trialEndIndex)==0);

        if (allTrials.eotCodes(trialEndIndex)==0) &&  (allTrials.certifiedNonInstruction(trialEndIndex)==1) ...
                && (allTrials.catchTrials(trialEndIndex)==1) % Work on only Correct Trials, which are not instruction or uncertified trials. [Vinay] - work even if it is a catch trial-changed this from 0 to 1
            
            % Get Eye Data
            eyeX = trials.eyeXData.data;
            eyeY = trials.eyeYData.data;
            % eyeStartTime = trials.eyeXData.timeMS(1);  % This is wrong.
            % The eye data is synchronized with trialStartTime.
            eyeStartTime = trials.trialStart.timeMS;
            eyeAllTimes = eyeStartTime + (0:(length(eyeX)-1))*(1000/Fs);
            
            stimOnTimes  = [trials.stimulusOnTime.timeMS];
            numStimuli = allTrials.targetPosAllTrials(trialEndIndex); %=length(stimOnTimes)/3;
            
            goodTrials.targetPos(correctIndex) = numStimuli;
            goodTrials.targetTime(correctIndex) = stimOnTimes(end);
            goodTrials.fixateMS(correctIndex) = trials.fixate.timeMS;
            goodTrials.fixonMS(correctIndex) = trials.fixOn.timeMS;
            goodTrials.stimOnTimes{correctIndex} = stimOnTimes;
            
            % Find position of Gabor1
            gaborPos = find([trials.stimDesc.data.gaborIndex]==0); % could be 4 gabors for GRF protocol
            if numStimuli>1  % The first one is not the target
                for j=1:numStimuli-1

                    stimTime = stimOnTimes(gaborPos(j));
                    stp=find(eyeAllTimes>=stimTime,1);
                    
                    stimData.stimOnsetTimeFromFixate(stimNumber) = stimTime-trials.fixate.timeMS;
                    stimData.stimPos(stimNumber) = j;
                    if (j==1) % First stimulus may not have sufficient baseline
                        eyeData(stimNumber).eyePosDataX = eyeX(stp:stp+eyeRangePos(2)-1);
                        eyeData(stimNumber).eyePosDataY = eyeY(stp:stp+eyeRangePos(2)-1);
                    else
                        eyeData(stimNumber).eyePosDataX = eyeX(stp+eyeRangePos(1):stp+eyeRangePos(2)-1);
                        eyeData(stimNumber).eyePosDataY = eyeY(stp+eyeRangePos(1):stp+eyeRangePos(2)-1);  
                    end
                    eyeData(stimNumber).eyeCal = trials.eyeCalibrationData.data.cal;
                    stimNumber=stimNumber+1;
                end
            end
            correctIndex=correctIndex+1;
        end
        trialEndIndex=trialEndIndex+1;
    end
    end
    save([folderExtract '/BehaviorData.mat'],'allTrials','goodTrials','stimData'); % Vinay - for linux
%     save([folderExtract '\EyeData.mat'],'eyeData','eyeRangeMS');
    save([folderExtract '/EyeData.mat'],'eyeData','eyeRangeMS'); % Vinay - for linux
else
    for i=1:numTrials
    %disp(i);
    clear trials
    trials = readLLFile('t',i);
    
    if isfield(trials,'trialEnd')
        allTrials.trialEnded(i) = 1;
        %allTrials.catchTrials(trialEndIndex) = trials.trial.data.catchTrial;
        %allTrials.instructTrials(trialEndIndex) = trials.trial.data.instructTrial;
        allTrials.trialCertify(trialEndIndex) = trials.trialCertify.data;
        %allTrials.targetPosAllTrials(trialEndIndex) = trials.trial.data.targetIndex+1;
        allTrials.eotCodes(trialEndIndex) = trials.trialEnd.data;
        
        allTrials.fixWindowSize(trialEndIndex) = trials.fixWindowData.data.windowDeg.size.width;
        %allTrials.respWindowSize(trialEndIndex) = trials.responseWindowData.data.windowDeg.size.width;
        %allTrials.certifiedNonInstruction(trialEndIndex) = (allTrials.instructTrials(trialEndIndex)==0)*(allTrials.trialCertify(trialEndIndex)==0);

        if (allTrials.eotCodes(trialEndIndex)==0) % Work on only Correct Trials, which are not instruction or uncertified trials. [Vinay] - work even if it is a catch trial-changed this from 0 to 1
            
            % Get Eye Data
            eyeX = trials.eyeLXData.data;
            eyeY = trials.eyeLYData.data;
            % eyeStartTime = trials.eyeXData.timeMS(1);  % This is wrong.
            % The eye data is synchronized with trialStartTime.
            eyeStartTime = trials.trialStart.timeMS;
            eyeAllTimes = eyeStartTime + (0:(length(eyeX)-1))*(1000/Fs);
            
            stimOnTimes  = [trials.fixOn.timeMS];
            %numStimuli = allTrials.targetPosAllTrials(trialEndIndex); %=length(stimOnTimes)/3;
            
            %goodTrials.targetPos(correctIndex) = numStimuli;
            %goodTrials.targetTime(correctIndex) = stimOnTimes(end);
            goodTrials.fixateMS(correctIndex) = trials.fixate.timeMS(1);
            goodTrials.fixonMS(correctIndex) = trials.fixOn.timeMS;
            %goodTrials.stimOnTimes{correctIndex} = stimOnTimes;
            
            % Find position of Gabor1
            %gaborPos = find([trials.stimDesc.data.gaborIndex]==1); % could be 4 gabors for GRF protocol
%             if numStimuli>1  % The first one is not the target
%                 for j=1:numStimuli-1

                    stimTime = trials.fixate.timeMS(1);
                    stp=find(eyeAllTimes>=stimTime,1);
                    
                    %stimData.stimOnsetTimeFromFixate(stimNumber) = stimTime-trials.fixate.timeMS;
                    %stimData.stimPos(stimNumber) = j;
%                     if (j==1) % First stimulus may not have sufficient baseline
                        eyeData(stimNumber).eyePosDataX = eyeX(stp:stp+eyeRangePos(2)-1);
                        eyeData(stimNumber).eyePosDataY = eyeY(stp:stp+eyeRangePos(2)-1);
%                     else
%                         eyeData(stimNumber).eyePosDataX = eyeX(stp+eyeRangePos(1):stp+eyeRangePos(2)-1);
%                         eyeData(stimNumber).eyePosDataY = eyeY(stp+eyeRangePos(1):stp+eyeRangePos(2)-1);  
%                     end
                    eyeData(stimNumber).eyeCal = trials.eyeLeftCalibrationData.data.cal;
                    stimNumber=stimNumber+1;
%                 end
            end
            correctIndex=correctIndex+1;
        end
        trialEndIndex=trialEndIndex+1;
    end
    save([folderExtract '/BehaviorData.mat'],'allTrials','goodTrials'); % Vinay - for linux
%     save([folderExtract '\EyeData.mat'],'eyeData','eyeRangeMS');
    save([folderExtract '/EyeData.mat'],'eyeData','eyeRangeMS'); % Vinay - for linux
end


% save([folderExtract '/BehaviorData.mat'],'allTrials','goodTrials','stimData'); % Vinay - for linux
% %     save([folderExtract '\EyeData.mat'],'eyeData','eyeRangeMS');
% save([folderExtract '/EyeData.mat'],'eyeData','eyeRangeMS'); % Vinay - for linux


folderName    = [folderSourceString 'data/' monkeyName '/' gridType '/' expDate '/' protocolName '/'];
folderExtract = [folderName 'extractedData/'];

clear eyeData 
load([folderExtract 'EyeData.mat']);

[eyeDataDegX,eyeDataDegY] = convertEyeDataToDeg(eyeData,1);

for i=1:length(eyeDataDegX)
plot(eyeDataDegX{i},eyeDataDegY{i}); axis([-2 2 -2 2]);
hold on;
disp(i);
pause;
end

for i=1:stimNumber-1
    m11Val(i) = eyeData(i).eyeCal.m11;
    m12Val(i) = eyeData(i).eyeCal.m12;
    m21Val(i) = eyeData(i).eyeCal.m21;
    m22Val(i) = eyeData(i).eyeCal.m22;
    tXVal(i) = eyeData(i).eyeCal.tX;
    tYVal(i) = eyeData(i).eyeCal.tY;
end