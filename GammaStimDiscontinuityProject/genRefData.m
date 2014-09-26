% Generate reference data - bipolar, average, CSD
% 10-20 EEG system
% Vinay Shirhatti, 19 September 2014
%==========================================================================

function genRefData(monkeyName,expDate,protocolName,folderSourceString,gridType,refType)


if isunix
    folderSourceString = appendIfNotPresent(folderSourceString,'/');
    folderNameMain = [folderSourceString 'data/' monkeyName '/' gridType '/' expDate '/' protocolName '/'];

    % Input folder
    lfpFolder  = [folderNameMain 'segmentedData/LFP/'];

    % Output folder
    if refType == 1 % bipolar
        outputFolder = [lfpFolder 'bipolar/'];
        makeDirectory(outputFolder);
    elseif refType == 2 % average
        outputFolder = lfpFolder;
        makeDirectory(outputFolder);
    elseif refType == 3 % csd
        outputFolder = [lfpFolder 'csd/'];
        makeDirectory(outputFolder);
    end
elseif ispc
    folderSourceString = appendIfNotPresent(folderSourceString,'\');
    folderNameMain = [folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName '\'];

    % Input folder
    lfpFolder  = [folderNameMain 'segmentedData\LFP\'];

    % Output folder
    if refType == 1 % bipolar
        outputFolder = [lfpFolder 'bipolar\'];
        makeDirectory(outputFolder);
    elseif refType == 2 % average
        outputFolder = lfpFolder;
        makeDirectory(outputFolder);
    elseif refType == 3 % csd
        outputFolder = [lfpFolder 'csd\'];
        makeDirectory(outputFolder);
    end
end

load([lfpFolder 'lfpInfo.mat']);

numElecs = length(electrodesStored);


if refType == 1
    disp('Generating bipolar reference data....');
    for i = 1:numElecs
        neighbourElec = findPair(electrodesStored(i),gridType); 
        for j = 1:length(neighbourElec)
            if electrodesStored(i)<neighbourElec(j)
                clear analogData1 analogData2 analogDataNotched1 analogDataNotched2
                load([lfpFolder 'elec' num2str(electrodesStored(i)) '.mat']);
                analogData1 = analogData;
                analogDataNotched1 = analogDataNotched;
                clear analogData analogDataNotched

                load([lfpFolder 'elec' num2str(neighbourElec(j)) '.mat']);
                analogData2 = analogData;
                analogDataNotched2 = analogDataNotched;
                clear analogData analogDataNotched

                analogData = analogData1 - analogData2;
                analogDataNotched  = analogDataNotched1 - analogDataNotched2;

                electrodesStoredString = [num2str(electrodesStored(i)) num2str(neighbourElec(j))];
                
                allStoredElectrodes(1,i) = str2num(electrodesStoredString);
                
                disp(['saving elec:' electrodesStoredString]);
                save([outputFolder 'elec' electrodesStoredString '.mat'],'analogData','analogDataNotched');
            else
                disp(['Already stored for this pair:' num2str(i) ',' num2str(j)]);
            end
        end
        
        % store w.r.t Cz
        neighbourElec = 18; % Cz
        clear analogData1 analogData2 analogDataNotched1 analogDataNotched2
        load([lfpFolder 'elec' num2str(electrodesStored(i)) '.mat']);
        analogData1 = analogData;
        analogDataNotched1 = analogDataNotched;
        clear analogData analogDataNotched

        load([lfpFolder 'elec' num2str(neighbourElec) '.mat']);
        analogData2 = analogData;
        analogDataNotched2 = analogDataNotched;
        clear analogData analogDataNotched

        analogData = analogData1 - analogData2;
        analogDataNotched  = analogDataNotched1 - analogDataNotched2;

        electrodesStoredString = [num2str(electrodesStored(i)) num2str(neighbourElec)];
        
        allStoredElectrodes(2,i) = str2num(electrodesStoredString);
        
        disp(['saving elec:' electrodesStoredString]);
        save([outputFolder 'elec' electrodesStoredString '.mat'],'analogData','analogDataNotched');
        
        electrodesStoredPair = allStoredElectrodes(1,:);
        electrodesStoredCz = allStoredElectrodes(2,:);
        
        electrodesStoredPair(electrodesStoredPair==0) = []; % to remove the instances where the calculation was skipped
        
        save([outputFolder 'lfpInfo.mat'],'electrodesStoredPair','electrodesStoredCz');
        
    end
    
elseif refType == 2 % average
    disp('Generating average reference data....');
    clear analogData analogNotchedData analogDataSum analogDataNotchedSum
    analogDataSum = [];
    analogDataNotchedSum = [];
    for i = 1:numElecs
        load([lfpFolder 'elec' num2str(electrodesStored(i)) '.mat']);
        if i == 1
            analogDataSum = analogData;
            analogDataNotchedSum = analogDataNotched;
        else
            analogDataSum= analogDataSum + analogData;
            analogDataNotchedSum= analogDataNotchedSum + analogDataNotched;
        end
    end
    analogData = analogDataSum/numElecs;
    analogDataNotched = analogDataNotchedSum/numElecs;
    
    disp('saving average reference data....');
    save([outputFolder 'avgRef.mat'],'analogData','analogDataNotched');
    
elseif refType == 3 % csd
    disp('Generating CSD data....');
    for i = 1:numElecs
        neighbourElec = findNeighbours(electrodesStored(i),gridType); 
        
        clear analogData1 analogData2 analogDataNotched1 analogDataNotched2
        load([lfpFolder 'elec' num2str(electrodesStored(i)) '.mat']);
        analogData1 = analogData;
        analogDataNotched1 = analogDataNotched;
        clear analogData analogDataNotched
        
        for j = 1:length(neighbourElec)
            
            load([lfpFolder 'elec' num2str(neighbourElec(j)) '.mat']);
            analogData2 = analogData1 - analogData/length(neighbourElec);
            analogDataNotched2 = analogDataNotched1 - analogData/length(neighbourElec);
            clear analogData analogDataNotched
            
            analogData = analogData2;
            analogDataNotched  = analogDataNotched2;
            
        end
        disp(['saving elec:' num2str(electrodesStored(i))]);
        save([outputFolder 'elec' num2str(electrodesStored(i)) '.mat'],'analogData','analogDataNotched');
    end
    
else
    disp('Give one of these options for refType: 1 (bipolar), 2 (average), 3 (csd)');
end

end

%--------------------------------------------------------------------------

function neighbourElec = findNeighbours(electrodeNum,gridType)

if strcmp(gridType,'EEG') % [Vinay] - adding the EEG 10-20 system grid
    
%     electrodeArray = ...
%        [00 01 00 02 00;
%         03 15 17 16 04;
%         05 13 18 14 06;
%         07 11 19 12 08;
%         00 09 00 10 00];
      
    if electrodeNum<1 || electrodeNum>21
        disp('Electrode Number out of range');
    else
        switch electrodeNum
            case 1
                neighbourElec = [2,3,15,17];
            case 2
                neighbourElec = [1,4,16,17];
            case 3
                neighbourElec = [1,5,13,15];
            case 4
                neighbourElec = [2,6,14,16];
            case 5
                neighbourElec = [3,7,13];
            case 6
                neighbourElec = [4,8,14];
            case 7
                neighbourElec = [5,9,11];
            case 8
                neighbourElec = [6,10,12];
            case 9
                neighbourElec = [7,10,11,19];
            case 10
                neighbourElec = [8,9,12,19];
            case 11
                neighbourElec = [7,9,13,19];
            case 12
                neighbourElec = [8,10,14,19];
            case 13
                neighbourElec = [5,11,15,18];
            case 14
                neighbourElec = [6,12,16,18];
            case 15
                neighbourElec = [1,3,13,17];
            case 16
                neighbourElec = [2,4,14,17];
            case 17
                neighbourElec = [1,2,15,16,18];
            case 18
                neighbourElec = [13,14,17,19];
            case 19
                neighbourElec = [9,10,11,12,18];
            case 20
                neighbourElec = [8,10,12,19];
            case 21
                neighbourElec = [7,9,11,19];
            otherwise
                neighbourElec = [];
        end
        
        
    end
else
    disp('Implemented only for EEG at the moment');
end    


end

%--------------------------------------------------------------------------

function neighbourElec = findPair(electrodeNum,gridType)

if strcmp(gridType,'EEG') % [Vinay] - adding the EEG 10-20 system grid
    
%     electrodeArray = ...
%        [00 01 00 02 00;
%         03 15 17 16 04;
%         05 13 18 14 06;
%         07 11 19 12 08;
%         00 09 00 10 00];
      
    if electrodeNum<1 || electrodeNum>21
        disp('Electrode Number out of range');
    else
        switch electrodeNum
            case 1
                neighbourElec = 2;
            case 2
                neighbourElec = 1;
            case 3
                neighbourElec = 4;
            case 4
                neighbourElec = 3;
            case 5
                neighbourElec = 6;
            case 6
                neighbourElec = 5;
            case 7
                neighbourElec = 8;
            case 8
                neighbourElec = 7;
            case 9
                neighbourElec = 10;
            case 10
                neighbourElec = 9;
            case 11
                neighbourElec = 12;
            case 12
                neighbourElec = 11;
            case 13
                neighbourElec = 14;
            case 14
                neighbourElec = 13;
            case 15
                neighbourElec = 16;
            case 16
                neighbourElec = 15;
            case 17
                neighbourElec = 19;
            case 18
                neighbourElec = 18;
            case 19
                neighbourElec = 17;
            case 20
                neighbourElec = 21;
            case 21
                neighbourElec = 20;
            otherwise
                neighbourElec = [];
        end
        
        
    end
else
    disp('Implemented only for EEG at the moment');
end    


end
            
         
