% This function is used to extract the digital data for the CRS Protocol.
% Vinay - modified from extractDigitalDataGRF
% 11 July 2014

function goodStimTimes = extractDigitalDataCRS(digitalEvents,digitalTimeStamps,folderExtract,ignoreTargetStimFlag,frameRate)

if ~exist('ignoreTargetStimFlag','var')   ignoreTargetStimFlag=0;       end
if ~exist('frameRate','var')              frameRate=60;                 end

useSingleITC18Flag=1;

ignoreCompErrors=0; % Vinay - set this to 1 if you want to ignore computer errors
hideSomeDigitalCode=1; % Vinay - set this to 1 if you have sent partial code based on the protocol

% Special cases in case a singleITC is used.
if useSingleITC18Flag
    % First, find the reward signals
    rewardOnPos = find(rem(digitalEvents,2)==0);
    rewardOffPos = find(digitalEvents==2^16-1);
    
    if length(rewardOnPos)~=length(rewardOffPos)
        disp('Unequal number of reward on and reward off!!');
    else
        rewardPos = [rewardOnPos(:) ; rewardOffPos(:)];
        disp([num2str(length(rewardPos)) ' are reward signals and will be discarded' ]);
        digitalEvents(rewardPos)=[];
        digitalTimeStamps(rewardPos)=[];
    end
    digitalEvents=digitalEvents-1; % [Vinay] - The last bit corresponds to 
    % the reward bit, it is 1 by default/OFF and goes to 0 when reward is 
    % ON. Subtract this bit to get the actual digital code
end

% All digital codes all start with a leading 1, which means that they are greater than hex2dec(8000) = 32768.
modifiedDigitalEvents = digitalEvents(digitalEvents>32768) - 32768;
allCodesInDec = unique(modifiedDigitalEvents);
disp(['Number of distinct codes: ' num2str(length(allCodesInDec))]);
allCodesInStr = convertDecCodeToStr(allCodesInDec,useSingleITC18Flag);

clear identifiedDigitalCodes badDigitalCodes
count=1; badCount=1;
for i=1:length(allCodesInDec)
    if ~digitalCodeDictionaryCRS(allCodesInStr(i,:)) % Vinay - for CRS
        disp(['Unidentified digital code: ' allCodesInStr(i,:) ', bin: ' dec2bin(allCodesInDec(i),16) ', dec: ' num2str(allCodesInDec(i)) ', occured ' num2str(length(find(modifiedDigitalEvents==allCodesInDec(i))))]);
        badDigitalCodes(badCount) = allCodesInDec(i);
        badCount=badCount+1;
    else
        identifiedDigitalCodes(count) = allCodesInDec(i);
        count=count+1;
    end
end

if badCount>1
    error(['The following Digital Codes are bad: ' num2str(badDigitalCodes)]);
end

numDigitalCodes = length(identifiedDigitalCodes);
disp(['Number of distinct codes identified: ' num2str(numDigitalCodes)]);

for i=1:numDigitalCodes
    digitalCodeInfo(i).codeNumber = identifiedDigitalCodes(i); %#ok<*AGROW>
    digitalCodeInfo(i).codeName = convertDecCodeToStr(identifiedDigitalCodes(i));
    clear codePos
    codePos = find(identifiedDigitalCodes(i) == digitalEvents-32768);
    digitalCodeInfo(i).time = digitalTimeStamps(codePos);
    digitalCodeInfo(i).value = digitalEvents(codePos+1);
end

% Write the digitalCodes
makeDirectory(folderExtract);
save([folderExtract 'digitalEvents.mat'],'digitalCodeInfo','digitalTimeStamps','digitalEvents');

%%%%%%%%%%%%%%%%%%%%%%% Get Stimulus results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
readDigitalCodesCRS(folderExtract,frameRate,hideSomeDigitalCode); % writes stimResults and trialResults
[goodStimNums,goodStimTimes] = getGoodStimNumsCRS(folderExtract,ignoreTargetStimFlag,ignoreCompErrors); % Good stimuli
getDisplayCombinationsCRS(folderExtract,goodStimNums);
save([folderExtract 'goodStimNums.mat'],'goodStimNums');

end

% GRF Specific protocols
function [stimResults,trialResults,trialEvents] = readDigitalCodesCRS(folderOut,frameRate,hideSomeDigitalCode)

if ~exist('frameRate','var')              frameRate=60;                 end
kForceQuit=7;

% Get the values of the following trial events for comparison with the dat
% file from lablib
trialEvents{1} = 'TS'; % Trial start
trialEvents{2} = 'TE'; % Trial End

folderOut = appendIfNotPresent(folderOut,'\');
load([folderOut 'digitalEvents.mat']);

allDigitalCodesInDec = [digitalCodeInfo.codeNumber];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the times and values of the events in trialEvents

for i=1:length(trialEvents)
    pos = find(convertStrCodeToDec(trialEvents{i})==allDigitalCodesInDec);
    if isempty(pos)
        disp(['Code ' trialEvents{i} ' not found!!']);
    else
        trialResults(i).times = [digitalCodeInfo(pos).time]; %#ok<*AGROW>
        trialResults(i).value = [digitalCodeInfo(pos).value];
    end
end

% Vinay - it seems that when the task is stopped it gives an extra trialEnd
% code. Have to snub this, else it gives an error.
trialEndError = 0; % initialize
if ((length(trialResults(2).times) - length(trialResults(1).times)) == 1)
    trialResults(2).times(length(trialResults(2).times)) = [];
    trialResults(2).value(length(trialResults(2).times)) = [];
    trialEndError = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stimulus properties
azimuth          = [digitalCodeInfo(find(convertStrCodeToDec('AZ')==allDigitalCodesInDec)).value];
elevation        = [digitalCodeInfo(find(convertStrCodeToDec('EL')==allDigitalCodesInDec)).value];
contrast         = [digitalCodeInfo(find(convertStrCodeToDec('CO')==allDigitalCodesInDec)).value];
temporalFrequency= [digitalCodeInfo(find(convertStrCodeToDec('TF')==allDigitalCodesInDec)).value];
radius           = [digitalCodeInfo(find(convertStrCodeToDec('RA')==allDigitalCodesInDec)).value];
sigma            = [digitalCodeInfo(find(convertStrCodeToDec('SI')==allDigitalCodesInDec)).value];
spatialFrequency = [digitalCodeInfo(find(convertStrCodeToDec('SF')==allDigitalCodesInDec)).value];
orientation      = [digitalCodeInfo(find(convertStrCodeToDec('OR')==allDigitalCodesInDec)).value];
spatialPhase     = [digitalCodeInfo(find(convertStrCodeToDec('SP')==allDigitalCodesInDec)).value]; % Vinay - added spatial phase becuase it is being used in CRS

%-------------------------------
% Vinay - adjust for the additional gain of 10 for the sent contrast data for the ring gabor
% Remove this part after making the corrections in CRSMap in XCode
contrast(contrast == 10000) = 1000;
contrast(contrast == 20000) = 2000;
temporalFrequency(temporalFrequency == 800) = 80;
%--------------------------------

% Get timing
trialStartTimes = [digitalCodeInfo(find(convertStrCodeToDec('TS')==allDigitalCodesInDec)).time];
taskGaborTimes  = [digitalCodeInfo(find(convertStrCodeToDec('TG')==allDigitalCodesInDec)).time];
mapping0Times   = [digitalCodeInfo(find(convertStrCodeToDec('M0')==allDigitalCodesInDec)).time];
mapping1Times   = [digitalCodeInfo(find(convertStrCodeToDec('M1')==allDigitalCodesInDec)).time];
mapping2Times   = [digitalCodeInfo(find(convertStrCodeToDec('M2')==allDigitalCodesInDec)).time]; % Vinay - added for CRS: the centre gabor
% Vinay, CRS - S = gabor0, R = gabor1, C = gabor2
numTrials = length(trialStartTimes);

% Vinay - Read the protocolNumber
stimResults.protocolNumber = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('PN')==allDigitalCodesInDec)).value]);

% Adjust the stimResults.parameter based on the protocolNumber if partial
% code was sent
if hideSomeDigitalCode
    switch stimResults.protocolNumber
        case 3
            azi = [azimuth,azimuth]; azimuth = reshape(azi',[],1);
            ele = [elevation,elevation]; elevation = reshape(ele',[],1);
            ori = [orientation,orientation]; orientation = reshape(ori',[],1);
            tf = [temporalFrequency,temporalFrequency]; temporalFrequency = reshape(tf',[],1);
            sig = [sigma,sigma]; sigma = reshape(sig',[],1);
            sf = [spatialFrequency,spatialFrequency]; spatialFrequency = reshape(sf',[],1);
            sp = [spatialPhase,spatialPhase]; spatialPhase = reshape(sp',[],1);
        case 4 % Dual Orientation Protocol
            azi = [azimuth,azimuth]; azimuth = reshape(azi',[],1);
            ele = [elevation,elevation]; elevation = reshape(ele',[],1);
            con = [contrast,contrast]; contrast = reshape(con',[],1);
            tf = [temporalFrequency,temporalFrequency]; temporalFrequency = reshape(tf',[],1);
            sig = [sigma,sigma]; sigma = reshape(sig',[],1);
            sf = [spatialFrequency,spatialFrequency]; spatialFrequency = reshape(sf',[],1);
            sp = [spatialPhase,spatialPhase]; spatialPhase = reshape(sp',[],1);
    end
end
            

if (max(diff([length(azimuth) length(elevation) length(contrast) length(temporalFrequency) ...
    length(radius) length(sigma) length(spatialFrequency) length(orientation) length(spatialPhase)])) > 0 )

    error('Length of stimulus properties are not even');
else
    if((length(azimuth) == length(mapping0Times)) && isempty(mapping1Times) && isempty(mapping2Times))
        disp('Only Mapping 0 is used i.e. S ON, C & R OFF'); % Vinay, CRS: S ON, C & R OFF
        stimResults.azimuth = convertUnits(azimuth',100);
        stimResults.elevation = convertUnits(elevation',100);
        stimResults.contrast = convertUnits(contrast',10);
        stimResults.temporalFrequency = convertUnits(temporalFrequency',100);
        stimResults.radius = convertUnits(radius',100);
        stimResults.sigma = convertUnits(sigma',100);
        stimResults.orientation = convertUnits(orientation');
        stimResults.spatialFrequency = convertUnits(spatialFrequency',100);
        stimResults.spatialPhase = convertUnits(spatialPhase'); % Vinay - for CRS
        stimResults.side=0;
        
    elseif((length(azimuth) == length(mapping1Times)) && isempty(mapping0Times) && isempty(mapping2Times))
        disp('Only Mapping 1 is used i.e. R ON, C & S OFF'); % Vinay, CRS: R ON, C & S OFF
        stimResults.azimuth = convertUnits(azimuth',100);
        stimResults.elevation = convertUnits(elevation',100);
        stimResults.contrast = convertUnits(contrast',10);
        stimResults.temporalFrequency = convertUnits(temporalFrequency',10);
        stimResults.radius = convertUnits(radius',100);
        stimResults.sigma = convertUnits(sigma',100);
        stimResults.orientation = convertUnits(orientation');
        stimResults.spatialFrequency = convertUnits(spatialFrequency',100);
        stimResults.spatialPhase = convertUnits(spatialPhase'); % Vinay - for CRS
        stimResults.side=1;
    
    elseif((length(azimuth) == length(mapping2Times)) && isempty(mapping0Times) && isempty(mapping2Times))
        disp('Only Mapping 2 is used i.e. C ON, R & S OFF'); % Vinay, CRS: C ON, R & S OFF
        stimResults.azimuth = convertUnits(azimuth',100);
        stimResults.elevation = convertUnits(elevation',100);
        stimResults.contrast = convertUnits(contrast',10);
        stimResults.temporalFrequency = convertUnits(temporalFrequency',10);
        stimResults.radius = convertUnits(radius',100);
        stimResults.sigma = convertUnits(sigma',100);
        stimResults.orientation = convertUnits(orientation');
        stimResults.spatialFrequency = convertUnits(spatialFrequency',100);
        stimResults.spatialPhase = convertUnits(spatialPhase'); % Vinay - for CRS
        stimResults.side=2;
    
    elseif((length(azimuth) == (length(mapping0Times) + length(mapping1Times))) && isempty(mapping2Times))
        disp('Only Mapping 0 and 1 are used i.e S & R ON, C OFF'); % Vinay, CRS: S & R ON, C OFF 
        stimResults.azimuth = convertUnits(azimuth',100);
        stimResults.elevation = convertUnits(elevation',100);
        stimResults.contrast = convertUnits(contrast',10);
        stimResults.temporalFrequency = convertUnits(temporalFrequency',10);
        stimResults.radius = convertUnits(radius',100);
        stimResults.sigma = convertUnits(sigma',100);
        stimResults.orientation = convertUnits(orientation');
        stimResults.spatialFrequency = convertUnits(spatialFrequency',100);
        stimResults.spatialPhase = convertUnits(spatialPhase'); % Vinay - for CRS
        stimResults.side=[0 1];
        
    elseif((length(azimuth) == (length(mapping0Times) + length(mapping2Times))) && isempty(mapping1Times))
        disp('Only Mapping 0 and 2 are used i.e. S & C ON, R OFF'); % Vinay, CRS: S & C ON, R OFF
        stimResults.azimuth = convertUnits(azimuth',100);
        stimResults.elevation = convertUnits(elevation',100);
        stimResults.contrast = convertUnits(contrast',10);
        stimResults.temporalFrequency = convertUnits(temporalFrequency',10);
        stimResults.radius = convertUnits(radius',100);
        stimResults.sigma = convertUnits(sigma',100);
        stimResults.orientation = convertUnits(orientation');
        stimResults.spatialFrequency = convertUnits(spatialFrequency',100);
        stimResults.spatialPhase = convertUnits(spatialPhase'); % Vinay - for CRS
        stimResults.side=[0 2];
        
    elseif((length(azimuth) == (length(mapping1Times) + length(mapping2Times))) && isempty(mapping0Times))
        disp('Only Mapping 1 and 2 are used i.e. R & C ON, S OFF'); % Vinay, CRS: R & C ON, S OFF
        stimResults.azimuth = convertUnits(azimuth',100);
        stimResults.elevation = convertUnits(elevation',100);
        stimResults.contrast = convertUnits(contrast',10);
        stimResults.temporalFrequency = convertUnits(temporalFrequency',10);
        stimResults.radius = convertUnits(radius',100);
        stimResults.sigma = convertUnits(sigma',100);
        stimResults.orientation = convertUnits(orientation');
        stimResults.spatialFrequency = convertUnits(spatialFrequency',100);
        stimResults.spatialPhase = convertUnits(spatialPhase'); % Vinay - for CRS
        stimResults.side=[1 2]; 
        
    else
        disp('Digital codes from all gabors!!!');
        stimResults.azimuth = convertUnits(azimuth',100);
        stimResults.elevation = convertUnits(elevation',100);
        stimResults.contrast = convertUnits(contrast',10);
        stimResults.temporalFrequency = convertUnits(temporalFrequency',10);
        stimResults.radius = convertUnits(radius',100);
        stimResults.sigma = convertUnits(sigma',100);
        stimResults.orientation = convertUnits(orientation');
        stimResults.spatialFrequency = convertUnits(spatialFrequency',100);
        stimResults.spatialPhase = convertUnits(spatialPhase'); % Vinay - for CRS
        stimResults.side=[0 1 2];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Vinay - Read the protocolNumber
% stimResults.protocolNumber = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('PN')==allDigitalCodesInDec)).value]);

% Instruction trials
instructionTrials = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('IT')==allDigitalCodesInDec)).value])';
if length(instructionTrials) ~= numTrials 
    error('Number of instruction trial entries different from numTrials');
end

% Catch trials
catchTrials = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('CT')==allDigitalCodesInDec)).value])';
if length(catchTrials) ~= numTrials 
    error('Number of catch trial entries different from numTrials');
end

% TrialCertify & TrialEnd (eotCode)
% These two entries may be repeated twice during force quit
trialCertify = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('TC')==allDigitalCodesInDec)).value])';
eotCodes = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('TE')==allDigitalCodesInDec)).value])';

forceQuits = find(eotCodes==kForceQuit);
numForceQuits = length(forceQuits);

% Vinay - if the final code of 'eotCodes' is 5 (False Alarm/quit) then it
% might be becuase of the stipulated blocks getting exhausted and the
% experiment stopping as a result of that. Check for this condition and if
% it is found to be so then truncate the last eotCode and trialCertify.
% Also, it seems that when the task is stopped it gives an extra trialEnd
% code if the blocks are not exhausted. So the last eotCodes ans
% trialCertify have to be dropped. This has been taken care of in the
% trialResults reading code earlier and the trialEndError indicates if this
% error occured or not
if ((eotCodes(length(eotCodes))==5) || (trialEndError == 1))
    eotCodes(length(eotCodes)) = [];
    trialCertify(length(eotCodes)) = [];
end

if length(eotCodes)-numForceQuits == numTrials
    disp(['numTrials: ' num2str(numTrials) ' numEotCodes: '  ...
        num2str(length(eotCodes)) ', ForceQuits: ' num2str(numForceQuits)]);
    goodEOTPos = find(eotCodes ~=kForceQuit);
    eotCodes = eotCodes(goodEOTPos);
    trialCertify = trialCertify(goodEOTPos);
else
     disp(['numTrials: ' num2str(numTrials) ' numEotCodes: '  ...
        num2str(length(eotCodes)) ', forcequits: ' num2str(numForceQuits)]);
    error('ForceQuit pressed after trial started'); % TODO - deal with this case
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numStimTask  = getStimPosPerTrial(trialStartTimes, taskGaborTimes);
numStimMap0  = getStimPosPerTrial(trialStartTimes, mapping0Times);
numStimMap1  = getStimPosPerTrial(trialStartTimes, mapping1Times);
numStimMap2  = getStimPosPerTrial(trialStartTimes, mapping2Times); % Vinay - for CRS

% % Check if Task and Mapping stim nums are the same for non-intruction trials
% nonInstructionTrials = find(instructionTrials==0);
% if (max(abs(numStimTask(nonInstructionTrials) - numStimMap0(nonInstructionTrials)))==0)
%     disp('Mapping0 and Task times are the same');
%     numStims = numStimMap0;
%     stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M0')==allDigitalCodesInDec)).time]';
% elseif (max(abs(numStimTask(nonInstructionTrials) - numStimMap1(nonInstructionTrials)))==0)
%     disp('Mapping1 and Task times are the same');
%     numStims = numStimMap1;
%     stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M1')==allDigitalCodesInDec)).time]';
% elseif (max(abs(numStimTask(nonInstructionTrials) - numStimMap2(nonInstructionTrials)))==0) % Vinay - for CRS
%     disp('Mapping2 and Task times are the same');
%     numStims = numStimMap2;
%     stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M2')==allDigitalCodesInDec)).time]';
% else
%     error('Mapping0/1/2 and Task times not the same');
% end

% Vinay - undo the checking ahead if the task gabor is hidden
taskGabor = 0; % Vinay - implies that the task gabor is visible or not in the task
if taskGabor
% Check if Task and Mapping stim nums are the same for non-intruction trials
    nonInstructionTrials = find(instructionTrials==0);
    if (max(abs(numStimTask(nonInstructionTrials) - numStimMap0(nonInstructionTrials)))==0)
        disp('Mapping0 and Task times are the same');
        numStims = numStimMap0;
        stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M0')==allDigitalCodesInDec)).time]';
    elseif (max(abs(numStimTask(nonInstructionTrials) - numStimMap1(nonInstructionTrials)))==0)
        disp('Mapping1 and Task times are the same');
        numStims = numStimMap1;
        stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M1')==allDigitalCodesInDec)).time]';
    elseif (max(abs(numStimTask(nonInstructionTrials) - numStimMap2(nonInstructionTrials)))==0) % Vinay - for centre gabor, CRS
        disp('Mapping2 and Task times are the same');
        numStims = numStimMap2;
        stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M2')==allDigitalCodesInDec)).time]';
    else
        error('Mapping0/1/2 and Task times not the same');
    end

    % elseif (sum(numStimMap0) == 0)
    %        numStims = numStimMap1;
    %        stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M1')==allDigitalCodesInDec)).time]';
    % else
    %     numStims = numStimMap0;
    %     stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M0')==allDigitalCodesInDec)).time]';    
    % end
else
    % Vinay - store the stim times of the C,R,S gabors
    stimResults.time0 = [digitalCodeInfo(find(convertStrCodeToDec('M0')==allDigitalCodesInDec)).time]';
    stimResults.time1 = [digitalCodeInfo(find(convertStrCodeToDec('M1')==allDigitalCodesInDec)).time]';
    stimResults.time2 = [digitalCodeInfo(find(convertStrCodeToDec('M2')==allDigitalCodesInDec)).time]';
    % Vinay - assign numStims as any of numStimMap0/1/2 based on whichever
    % is used. The numbers for the gabors which are drawn should be equal,
    % so it doesn't matter which one is assigned to numStims. numStims says
    % how many stimuli were shown on each of the trials
    if (sum(numStimMap0) ~= 0) % Vinay - this would be zero only if gabor0 (surround) were hidden. Similarly for other gabors.
           numStims = numStimMap0;
           % [Vinay] - define another variable 'time' in stimResults to
           % store the onset time of stimulus display. This variable is
           % assigned the times corresponding to gabor0(surround) if it
           % is drawn (since it is the fist one to be drawn), else the
           % times corresponding to gabor1 (ring) and if even gabor1 is
           % undrawn then the times of gabor2
           stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M0')==allDigitalCodesInDec)).time]';
    elseif (sum(numStimMap1) ~= 0)
           numStims = numStimMap1;
           stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M1')==allDigitalCodesInDec)).time]';
    elseif (sum(numStimMap2) ~= 0)
           numStims = numStimMap2;
           stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M2')==allDigitalCodesInDec)).time]';
    else
        error('No CRS gabor stimulus shown');
    end
end

taskType = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('TG')==allDigitalCodesInDec)).value])';
posTask = 0;
pos=0;
for i=1:numTrials
    taskTypeThisTrial = taskType(posTask+1:posTask+numStimTask(i));
    if (numStims(i)>0)
        stimResults.type(pos+1:pos+numStims(i)) = taskTypeThisTrial;
        stimResults.trialNumber(pos+1:pos+numStims(i)) = i;
        stimResults.stimPosition(pos+1:pos+numStims(i)) = 1:numStims(i);
        
%         if stimResults.side==0
%             stimResults.stimOnFrame(pos+1:pos+numStims(i)) = ...
%                 (mapping0Times(pos+1:pos+numStims(i)) - mapping0Times(pos+1))*frameRate;
%         elseif stimResults.side==1
%             stimResults.stimOnFrame(pos+1:pos+numStims(i)) = ...
%                 (mapping1Times(pos+1:pos+numStims(i)) - mapping1Times(pos+1))*frameRate;
%         elseif stimResults.side==2 % Vinay - for CRS
%             stimResults.stimOnFrame(pos+1:pos+numStims(i)) = ...
%                 (mapping2Times(pos+1:pos+numStims(i)) - mapping2Times(pos+1))*frameRate;
%         end
        
        if ~isempty(stimResults.side(find(stimResults.side)==0)) % Vinay - checks if the corresponding gabor (surround) is present or not
            stimResults.stimOnFrame0(pos+1:pos+numStims(i)) = ...
                (mapping0Times(pos+1:pos+numStims(i)) - mapping0Times(pos+1))*frameRate; % Vinay - gives the frame where the corresponding stimulus ON occured
        end
        if ~isempty(stimResults.side(find(stimResults.side)==1)) % checks if the corresponding gabor (ring) is present or not
            stimResults.stimOnFrame1(pos+1:pos+numStims(i)) = ...
                (mapping1Times(pos+1:pos+numStims(i)) - mapping1Times(pos+1))*frameRate;
        end
        if ~isempty(stimResults.side(find(stimResults.side)==2)) % checks if the corresponding gabor (centre) is present or not
            stimResults.stimOnFrame2(pos+1:pos+numStims(i)) = ...
                (mapping2Times(pos+1:pos+numStims(i)) - mapping2Times(pos+1))*frameRate;
        end
        
        stimResults.instructionTrials(pos+1:pos+numStims(i)) = instructionTrials(i); %always zero
        stimResults.catch(pos+1:pos+numStims(i)) = catchTrials(i);
        stimResults.eotCodes(pos+1:pos+numStims(i)) = eotCodes(i);
        stimResults.trialCertify(pos+1:pos+numStims(i)) = trialCertify(i);
        pos = pos+numStims(i);
    end
    posTask = posTask+numStimTask(i);
end

% Save in folderOut
save([folderOut 'stimResults.mat'],'stimResults');
save([folderOut 'trialResults.mat'],'trialEvents','trialResults');

end
function [goodStimNums,goodStimTimes] = getGoodStimNumsCRS(folderOut,ignoreTargetStimFlag,ignoreCompErrors)

if ~exist('ignoreTargetStimFlag','var')       ignoreTargetStimFlag=0;   end

folderOut = appendIfNotPresent(folderOut,'\');
load([folderOut 'stimResults.mat']);

totalStims = length(stimResults.eotCodes);
disp(['Number of trials: ' num2str(max(stimResults.trialNumber))]);
disp(['Number of stimuli: ' num2str(totalStims)]);

% exclude uncertified trials, catch trials and instruction trials
tc = find(stimResults.trialCertify==1);
if ignoreCompErrors
    tc = []; % [Vinay] - not checking for computer errors
end
it = find(stimResults.instructionTrials==1);
ct = find(stimResults.catch==1);

if ignoreTargetStimFlag
    badStimNums = [it tc]; % catch trials are now considered good
else
    badStimNums = [it tc ct];
end

%eottypes
% 0 - correct, 1 - wrong, 2-failed, 3-broke, 4-ignored, 5-False
% Alarm/quit, 6 - distracted, 7 - force quit
%disp('Analysing correct, wrong and failed trials');
%badEOTs = find(stimResults.eotCodes>2); 
%disp('Analysing correct and wrong trials')
%badEOTs = find(stimResults.eotCodes>1); 
disp('Analysing only correct trials')
badEOTs = find(stimResults.eotCodes>0); 
badStimNums = [badStimNums badEOTs];

goodStimNums = setdiff(1:totalStims,unique(badStimNums));

% stim types
% 0 - Null, 1 - valid, 2 - target, 3 - frontpadding, 4 - backpadding
if ~ignoreTargetStimFlag
    disp('Only taking valid stims ');
    validStims = find(stimResults.type==1);
    goodStimNums = intersect(goodStimNums,validStims);
    
    %%%%%%%%%%%%%% Remove bad stimuli after target %%%%%%%%%%%%%%%%%%%%
    
    clear trialNums stimPos
    trialNums = stimResults.trialNumber(goodStimNums);
    stimPos   = stimResults.stimPosition(goodStimNums);
    
    % Get the target positions of the trialNums
    clear goodTrials
    goodTrials = unique(trialNums);
    
    clear targetPos
    for i=1:length(goodTrials)
        allStimWithThisTrialNum = find(stimResults.trialNumber==goodTrials(i));
        
        if sum(stimResults.catch(allStimWithThisTrialNum))>0        % catch trials
            targetPos(trialNums==goodTrials(i)) = inf; %#ok<*AGROW>
        else
            targetPos(trialNums==goodTrials(i)) = find(stimResults.type(allStimWithThisTrialNum)==2);
        end
    end
    
    validStimuliAfterTarget = find(stimPos>targetPos);
    if ~isempty(validStimuliAfterTarget)
        disp([num2str(length(validStimuliAfterTarget)) ' out of ' num2str(length(goodStimNums)) ' stimuli after target']);
        save([folderOut 'validStimAfterTarget.mat'],'validStimuliAfterTarget');
    end
    
    goodStimNums(validStimuliAfterTarget)=[];
end
disp(['Number of good stimuli: ' num2str(length(goodStimNums))]);
goodStimTimes = stimResults.time(goodStimNums);
end
function parameterCombinations = getDisplayCombinationsCRS(folderOut,goodStimNums)

folderOut = appendIfNotPresent(folderOut,'\');
load([folderOut 'stimResults.mat']);

% Nine parameters are chosen: for CRS [Vinay]
% 1. Azimuth
% 2. Elevation
% 3. Sigma
% 4. Spatial Frequency
% 5. Orientation
% 6. Contrast
% 7. Temporal Frequency
% 8. Spatial Phase
% 9. Radius

% Parameters index
parameters{1} = 'azimuth';
parameters{2} = 'elevation';
parameters{3} = 'sigma';
parameters{4} = 'spatialFrequency';
parameters{5} = 'orientation';
parameters{6} = 'contrast';
parameters{7} = 'temporalFrequency'; %#ok<NASGU>
parameters{8} = 'spatialPhase'; % Vinay - for CRS
parameters{9} = 'radius'; % Vinay - for CRS
parameters{10} = 'kGaborNumber';

% [Vinay] - read the parameter values for all the gabors together
aValsAllGabor  = stimResults.azimuth;
eValsAllGabor  = stimResults.elevation;
sValsAllGabor  = stimResults.sigma;
fValsAllGabor  = stimResults.spatialFrequency;
oValsAllGabor  = stimResults.orientation;
cValsAllGabor  = stimResults.contrast;
tValsAllGabor  = stimResults.temporalFrequency;
pValsAllGabor  = stimResults.spatialPhase; % Vinay - for CRS
rValsAllGabor  = stimResults.radius; % Vinay - for CRS

% [Vinay] - The parameters stored above are for all the gabors together.
% They appear in repeated ordered groups (gabor0,1,2). So the relevant
% parameter values have to be picked corresponding to individual gabors.

% [Vinay] - If there are any null gabors then assign '0' values to its
% parameters and all lengths as 1

numOfGabors = length(stimResults.side);

if numOfGabors ~= 3
    
    % Vinay - If the number of gabors is not 3 there are null gabors.
    % Reshape the ValsAllGabor matrices so that each row corresponds to the
    % values for a particular gabor.
    aTemp = reshape(aValsAllGabor,numOfGabors,length(aValsAllGabor)/numOfGabors);
    eTemp = reshape(eValsAllGabor,numOfGabors,length(eValsAllGabor)/numOfGabors);
    sTemp = reshape(sValsAllGabor,numOfGabors,length(sValsAllGabor)/numOfGabors);
    fTemp = reshape(fValsAllGabor,numOfGabors,length(fValsAllGabor)/numOfGabors);
    oTemp = reshape(oValsAllGabor,numOfGabors,length(oValsAllGabor)/numOfGabors);
    cTemp = reshape(cValsAllGabor,numOfGabors,length(cValsAllGabor)/numOfGabors);
    tTemp = reshape(tValsAllGabor,numOfGabors,length(tValsAllGabor)/numOfGabors);
    pTemp = reshape(pValsAllGabor,numOfGabors,length(pValsAllGabor)/numOfGabors);
    rTemp = reshape(rValsAllGabor,numOfGabors,length(rValsAllGabor)/numOfGabors);
    
    countGabor = 0;
    for gi = 1:3
        if isempty(stimResults.side(find(stimResults.side)==(gi-1))) % check if this gabor is null or not
            % If it is null then assign '0's to the corresponding row in a
            % temp matrix below
            aValsAllGaborTmp(gi,1:size(aTemp,2))  = 0;
            eValsAllGaborTmp(gi,1:size(eTemp,2))  = 0;
            sValsAllGaborTmp(gi,1:size(sTemp,2))  = 0;
            fValsAllGaborTmp(gi,1:size(fTemp,2))  = 0;
            oValsAllGaborTmp(gi,1:size(oTemp,2))  = 0;
            cValsAllGaborTmp(gi,1:size(cTemp,2))  = 0;
            tValsAllGaborTmp(gi,1:size(tTemp,2))  = 0;
            pValsAllGaborTmp(gi,1:size(pTemp,2))  = 0;
            rValsAllGaborTmp(gi,1:size(rTemp,2))  = 0;
        else
            % If it is not null then extract the corresponding row from the
            % aTemp (or the related one) matrix and assign the values to
            % the corresponding row of the Temp matrix below
            countGabor=countGabor+1;
            aValsAllGaborTmp(gi,1:size(aTemp,2))  = aTemp(countGabor,1:size(aTemp,2));
            eValsAllGaborTmp(gi,1:size(eTemp,2))  = eTemp(countGabor,1:size(eTemp,2));
            sValsAllGaborTmp(gi,1:size(sTemp,2))  = sTemp(countGabor,1:size(sTemp,2));
            fValsAllGaborTmp(gi,1:size(fTemp,2))  = fTemp(countGabor,1:size(fTemp,2));
            oValsAllGaborTmp(gi,1:size(oTemp,2))  = oTemp(countGabor,1:size(oTemp,2));
            cValsAllGaborTmp(gi,1:size(cTemp,2))  = cTemp(countGabor,1:size(cTemp,2));
            tValsAllGaborTmp(gi,1:size(tTemp,2))  = tTemp(countGabor,1:size(tTemp,2));
            pValsAllGaborTmp(gi,1:size(pTemp,2))  = pTemp(countGabor,1:size(pTemp,2));
            rValsAllGaborTmp(gi,1:size(rTemp,2))  = rTemp(countGabor,1:size(rTemp,2));
        end
    end
    
    % Redefine the ValsAllGabor matrices by reshaping the Temp matrix 
    % obtained above 
    aValsAllGabor = reshape(aValsAllGaborTmp,1,[]);
    eValsAllGabor = reshape(eValsAllGaborTmp,1,[]);
    sValsAllGabor = reshape(sValsAllGaborTmp,1,[]);
    fValsAllGabor = reshape(fValsAllGaborTmp,1,[]);
    oValsAllGabor = reshape(oValsAllGaborTmp,1,[]);
    cValsAllGabor = reshape(cValsAllGaborTmp,1,[]);
    tValsAllGabor = reshape(tValsAllGaborTmp,1,[]);
    pValsAllGabor = reshape(pValsAllGaborTmp,1,[]);
    rValsAllGabor = reshape(rValsAllGaborTmp,1,[]);
    
end
    
for gaborNum = 1:3
    
    %     gaborIndices = gaborNum:numOfGabors:length(aValsAllGabor); % [Vinay] - index values corresponding to 'gaborNum' gabor
    % Vinay - since the ValsAllGood matrices have been redefined to contain the
    % information for all the gabors we increment by 3 instead of
    % numOfGabors
    
    gaborIndices = gaborNum:3:length(aValsAllGabor); % [Vinay] - index values corresponding to 'gaborNum' gabor
    
    aValsAll = aValsAllGabor(gaborIndices);
    eValsAll = eValsAllGabor(gaborIndices);
    sValsAll = sValsAllGabor(gaborIndices);
    fValsAll = fValsAllGabor(gaborIndices);
    oValsAll = oValsAllGabor(gaborIndices);
    cValsAll = cValsAllGabor(gaborIndices);
    tValsAll = tValsAllGabor(gaborIndices);
    pValsAll = pValsAllGabor(gaborIndices); % Vinay - for CRS
    rValsAll = rValsAllGabor(gaborIndices); % Vinay - for CRS
    
    if ~isempty(aValsAll)
        % Get good stim
        if ~exist('goodStimNums','var')
            goodStimNums = getGoodStimNumsCRS(folderOut);
        end

        aValsGood = aValsAll(goodStimNums);
        eValsGood = eValsAll(goodStimNums);
        sValsGood = sValsAll(goodStimNums);
        fValsGood = fValsAll(goodStimNums);
        oValsGood = oValsAll(goodStimNums);
        cValsGood = cValsAll(goodStimNums);
        tValsGood = tValsAll(goodStimNums);
        pValsGood = pValsAll(goodStimNums); % Vinay - for CRS
        rValsGood = rValsAll(goodStimNums); % Vinay - for CRS

        aValsUnique{gaborNum} = unique(aValsGood); aLen = length(aValsUnique{gaborNum});
        eValsUnique{gaborNum} = unique(eValsGood); eLen = length(eValsUnique{gaborNum});
        sValsUnique{gaborNum} = unique(sValsGood); sLen = length(sValsUnique{gaborNum});
        fValsUnique{gaborNum} = unique(fValsGood); fLen = length(fValsUnique{gaborNum});
        oValsUnique{gaborNum} = unique(oValsGood); oLen = length(oValsUnique{gaborNum});
        cValsUnique{gaborNum} = unique(cValsGood); cLen = length(cValsUnique{gaborNum});
        tValsUnique{gaborNum} = unique(tValsGood); tLen = length(tValsUnique{gaborNum});
        pValsUnique{gaborNum} = unique(pValsGood); pLen = length(pValsUnique{gaborNum}); % Vinay - for CRS
        rValsUnique{gaborNum} = unique(rValsGood); rLen = length(rValsUnique{gaborNum}); % Vinay - for CRS

        % display
        disp(['Gabor ' num2str(gaborNum - 1)]); % Vinay - display kGabor number
        disp(['Number of unique azimuths: ' num2str(aLen)]);
        disp(['Number of unique elevations: ' num2str(eLen)]);
        disp(['Number of unique sigmas: ' num2str(sLen)]);
        disp(['Number of unique Spatial freqs: ' num2str(fLen)]);
        disp(['Number of unique orientations: ' num2str(oLen)]);
        disp(['Number of unique contrasts: ' num2str(cLen)]);
        disp(['Number of unique temporal freqs: ' num2str(tLen)]);
        disp(['Number of unique spatial phases: ' num2str(pLen)]); % Vinay - for CRS
        disp(['Number of unique radii: ' num2str(rLen)]); % Vinay - for CRS

        % If more than one value, make another entry with all values
        if (aLen > 1)           aLen=aLen+1;                    end
        if (eLen > 1)           eLen=eLen+1;                    end
        if (sLen > 1)           sLen=sLen+1;                    end
        if (fLen > 1)           fLen=fLen+1;                    end
        if (oLen > 1)           oLen=oLen+1;                    end
        if (cLen > 1)           cLen=cLen+1;                    end
        if (tLen > 1)           tLen=tLen+1;                    end
        if (pLen > 1)           pLen=pLen+1;                    end % Vinay - for CRS
        if (rLen > 1)           rLen=rLen+1;                    end % Vinay - for CRS

        allPos = 1:length(goodStimNums);
        disp(['total combinations: ' num2str((aLen)*(eLen)*(sLen)*(fLen)*(oLen)*(cLen)*(tLen)*(pLen)*(rLen))]);

        for a=1:aLen
            if a==aLen
                aPos = allPos;
            else
                aPos = find(aValsGood == aValsUnique{gaborNum}(a));
            end

            for e=1:eLen
                if e==eLen
                    ePos = allPos;
                else
                    ePos = find(eValsGood == eValsUnique{gaborNum}(e));
                end

                for s=1:sLen
                    if s==sLen
                        sPos = allPos;
                    else
                        sPos = find(sValsGood == sValsUnique{gaborNum}(s));
                    end

                    for f=1:fLen
                        if f==fLen
                            fPos = allPos;
                        else
                            fPos = find(fValsGood == fValsUnique{gaborNum}(f));
                        end

                        for o=1:oLen
                            if o==oLen
                                oPos = allPos;
                            else
                                oPos = find(oValsGood == oValsUnique{gaborNum}(o));
                            end

                            for c=1:cLen
                                if c==cLen
                                    cPos = allPos;
                                else
                                    cPos = find(cValsGood == cValsUnique{gaborNum}(c));
                                end

                                for t=1:tLen
                                    if t==tLen
                                        tPos = allPos;
                                    else
                                        tPos = find(tValsGood == tValsUnique{gaborNum}(t));
                                    end

                                    for p=1:pLen % Vinay - for CRS
                                        if p==pLen
                                            pPos = allPos;
                                        else
                                            pPos = find(pValsGood == pValsUnique{gaborNum}(p));
                                        end

                                        for r=1:rLen % Vinay - for CRS
                                            if r==rLen
                                                rPos = allPos;
                                            else
                                                rPos = find(rValsGood == rValsUnique{gaborNum}(r));
                                            end


                                    aePos = intersect(aPos,ePos);
                                    aesPos = intersect(aePos,sPos);
                                    aesfPos = intersect(aesPos,fPos);
                                    aesfoPos = intersect(aesfPos,oPos);
                                    aesfocPos = intersect(aesfoPos,cPos);
                                    aesfoctPos = intersect(aesfocPos,tPos);
                                    aesfoctpPos = intersect(aesfoctPos,pPos); % Vinay - for CRS
                                    aesfoctprPos = intersect(aesfoctpPos,rPos); % Vinay - for CRS

                                    parameterCombinations{a,e,s,f,o,c,t,p,r,gaborNum} = aesfoctprPos; %#ok<AGROW> % [Vinay] - added additional parameters p and r for CRS
                                    % [Vinay] - the last dimension has been
                                    % added to save the gabor number
                                    % corresponding to its parameter
                                    % combinations. gabor0 (surround) has index
                                    % value 1, gabor1 (ring) -> 2,
                                    % gabor2(centre) -> 3
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
%     % save
%     save([folderOut 'parameterCombinations.mat'],'parameters','parameterCombinations', ...
%         'aValsUnique','eValsUnique','sValsUnique','fValsUnique','oValsUnique','cValsUnique','tValsUnique','pValsUnique','rValsUnique'); % [Vinay] - added additional parameters p and r for CRS

end

    % save
    save([folderOut 'parameterCombinations.mat'],'parameters','parameterCombinations', ...
        'aValsUnique','eValsUnique','sValsUnique','fValsUnique','oValsUnique','cValsUnique','tValsUnique','pValsUnique','rValsUnique'); % [Vinay] - added additional parameters p and r for CRS

end
    
function outNum = convertUnits(num,f,useSingleITC18Flag)

if ~exist('f','var')                        f=1;                        end
if ~exist('useSingleITC18Flag','var')       useSingleITC18Flag=1;       end

for i=1:length(num)
    if num(i) > 16384
        num(i)=num(i)-32768;
    end
end
outNum = num/f;

if useSingleITC18Flag
    outNum=outNum/2;
end
end
function [numStim,stimOnPos] = getStimPosPerTrial(trialStartTimes, stimStartTimes)

numTrials = length(trialStartTimes);

stimOnPos = cell(1,numTrials);
numStim   = zeros(1,numTrials);

for i=1:numTrials-1
    stimOnPos{i} = intersect(find(stimStartTimes>=trialStartTimes(i)),find(stimStartTimes<trialStartTimes(i+1)));
    numStim(i) = length(stimOnPos{i});
end
stimOnPos{numTrials} = find(stimStartTimes>=trialStartTimes(numTrials));
numStim(numTrials) = length(stimOnPos{numTrials});
end
