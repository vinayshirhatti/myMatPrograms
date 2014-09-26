% This function is used to extract the digital data for the CRS Protocol.
% Vinay - modified from extractDigitalDataGRF
% 11 July 2014

function goodStimTimes = extractDigitalDataCRS(digitalEvents,digitalTimeStamps,folderExtract,ignoreTargetStimFlag,frameRate)

if ~exist('ignoreTargetStimFlag','var')   ignoreTargetStimFlag=0;       end
if ~exist('frameRate','var')              frameRate=60;                 end

useSingleITC18Flag=1;

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
    digitalEvents=digitalEvents-1;
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
readDigitalCodesCRS(folderExtract,frameRate); % writes stimResults and trialResults
[goodStimNums,goodStimTimes] = getGoodStimNumsCRS(folderExtract,ignoreTargetStimFlag); % Good stimuli
getDisplayCombinationsCRS(folderExtract,goodStimNums);
save([folderExtract 'goodStimNums.mat'],'goodStimNums');

end

% GRF Specific protocols
function [stimResults,trialResults,trialEvents] = readDigitalCodesCRS(folderOut,frameRate)

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


% Get timing
trialStartTimes = [digitalCodeInfo(find(convertStrCodeToDec('TS')==allDigitalCodesInDec)).time];
taskGaborTimes  = [digitalCodeInfo(find(convertStrCodeToDec('TG')==allDigitalCodesInDec)).time];
mapping0Times   = [digitalCodeInfo(find(convertStrCodeToDec('M0')==allDigitalCodesInDec)).time];
mapping1Times   = [digitalCodeInfo(find(convertStrCodeToDec('M1')==allDigitalCodesInDec)).time];
mapping2Times   = [digitalCodeInfo(find(convertStrCodeToDec('M2')==allDigitalCodesInDec)).time]; % Vinay - added for CRS: the centre gabor
numTrials = length(trialStartTimes);

% Check the default case - only mapping0/1 is on, and only its stimulus properties are put out.

if (max(diff([length(azimuth) length(elevation) length(contrast) length(temporalFrequency) ...
    length(radius) length(sigma) length(spatialFrequency) length(orientation) length(spatialPhase)])) > 0 )

    error('Length of stimulus properties are not even');
else
    if((length(azimuth) == length(mapping0Times)) && isempty(mapping1Times) && isempty(mapping2Times))
        disp('Only Mapping 0 is used');
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
        disp('Only Mapping 1 is used');
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
        disp('Only Mapping 2 is used');
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
    
    elseif((length(azimuth) == length(mapping0Times)) && (length(azimuth)== length(mapping1Times)) && isempty(mapping2Times))
        disp('Only Mapping 0 and 1 are used');
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
        
    elseif((length(azimuth) == length(mapping0Times)) && (length(azimuth)== length(mapping2Times)) && isempty(mapping1Times))
        disp('Only Mapping 0 and 2 are used');
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
        
    elseif((length(azimuth) == length(mapping1Times)) && (length(azimuth)== length(mapping2Times)) && isempty(mapping0Times))
        disp('Only Mapping 1 and 2 are used');
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
elseif (max(abs(numStimTask(nonInstructionTrials) - numStimMap2(nonInstructionTrials)))==0) % Vinay - for CRS
    disp('Mapping2 and Task times are the same');
    numStims = numStimMap2;
    stimResults.time = [digitalCodeInfo(find(convertStrCodeToDec('M2')==allDigitalCodesInDec)).time]';
else
    error('Mapping0/1/2 and Task times not the same');
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
        
        if stimResults.side==0
            stimResults.stimOnFrame(pos+1:pos+numStims(i)) = ...
                (mapping0Times(pos+1:pos+numStims(i)) - mapping0Times(pos+1))*frameRate;
        elseif stimResults.side==1
            stimResults.stimOnFrame(pos+1:pos+numStims(i)) = ...
                (mapping1Times(pos+1:pos+numStims(i)) - mapping1Times(pos+1))*frameRate;
        elseif stimResults.side==2 % Vinay - for CRS
            stimResults.stimOnFrame(pos+1:pos+numStims(i)) = ...
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
function [goodStimNums,goodStimTimes] = getGoodStimNumsCRS(folderOut,ignoreTargetStimFlag)

if ~exist('ignoreTargetStimFlag','var')       ignoreTargetStimFlag=0;   end

folderOut = appendIfNotPresent(folderOut,'\');
load([folderOut 'stimResults.mat']);

totalStims = length(stimResults.eotCodes);
disp(['Number of trials: ' num2str(max(stimResults.trialNumber))]);
disp(['Number of stimuli: ' num2str(totalStims)]);

% exclude uncertified trials, catch trials and instruction trials
tc = find(stimResults.trialCertify==1);
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

% get Contrast
aValsAll  = stimResults.azimuth;
eValsAll  = stimResults.elevation;
sValsAll  = stimResults.sigma;
fValsAll  = stimResults.spatialFrequency;
oValsAll  = stimResults.orientation;
cValsAll  = stimResults.contrast;
tValsAll  = stimResults.temporalFrequency;
pValsAll  = stimResults.spatialPhase; % Vinay - for CRS
rValsAll  = stimResults.radius; % Vinay - for CRS

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

    aValsUnique = unique(aValsGood); aLen = length(aValsUnique);
    eValsUnique = unique(eValsGood); eLen = length(eValsUnique);
    sValsUnique = unique(sValsGood); sLen = length(sValsUnique);
    fValsUnique = unique(fValsGood); fLen = length(fValsUnique);
    oValsUnique = unique(oValsGood); oLen = length(oValsUnique);
    cValsUnique = unique(cValsGood); cLen = length(cValsUnique);
    tValsUnique = unique(tValsGood); tLen = length(tValsUnique);
    pValsUnique = unique(pValsGood); pLen = length(pValsUnique); % Vinay - for CRS
    rValsUnique = unique(rValsGood); rLen = length(rValsUnique); % Vinay - for CRS

    % display
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
            aPos = find(aValsGood == aValsUnique(a));
        end

        for e=1:eLen
            if e==eLen
                ePos = allPos;
            else
                ePos = find(eValsGood == eValsUnique(e));
            end

            for s=1:sLen
                if s==sLen
                    sPos = allPos;
                else
                    sPos = find(sValsGood == sValsUnique(s));
                end

                for f=1:fLen
                    if f==fLen
                        fPos = allPos;
                    else
                        fPos = find(fValsGood == fValsUnique(f));
                    end

                    for o=1:oLen
                        if o==oLen
                            oPos = allPos;
                        else
                            oPos = find(oValsGood == oValsUnique(o));
                        end
                        
                        for c=1:cLen
                            if c==cLen
                                cPos = allPos;
                            else
                                cPos = find(cValsGood == cValsUnique(c));
                            end
                            
                            for t=1:tLen
                                if t==tLen
                                    tPos = allPos;
                                else
                                    tPos = find(tValsGood == tValsUnique(t));
                                end
                                
                                for p=1:pLen % Vinay - for CRS
                                    if p==pLen
                                        pPos = allPos;
                                    else
                                        pPos = find(pValsGood == pValsUnique(p));
                                    end
                                    
                                    for r=1:rLen % Vinay - for CRS
                                        if r==rLen
                                            rPos = allPos;
                                        else
                                            rPos = find(rValsGood == rValsUnique(r));
                                        end


                                aePos = intersect(aPos,ePos);
                                aesPos = intersect(aePos,sPos);
                                aesfPos = intersect(aesPos,fPos);
                                aesfoPos = intersect(aesfPos,oPos);
                                aesfocPos = intersect(aesfoPos,cPos);
                                aesfoctPos = intersect(aesfocPos,tPos);
                                aesfoctpPos = intersect(aesfoctPos,pPos); % Vinay - for CRS
                                aesfoctprPos = intersect(aesfoctpPos,rPos); % Vinay - for CRS
                                
                                parameterCombinations{a,e,s,f,o,c,t,p,r} = aesfoctprPos; %#ok<AGROW> % [Vinay] - added additional parameters p and r for CRS
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    % save
    save([folderOut 'parameterCombinations.mat'],'parameters','parameterCombinations', ...
        'aValsUnique','eValsUnique','sValsUnique','fValsUnique','oValsUnique','cValsUnique','tValsUnique','pValsUnique','rValsUnique'); % [Vinay] - added additional parameters p and r for CRS

end
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
