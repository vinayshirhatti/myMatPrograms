%==========================================================================
% Vinay Shirhatti, 03 September 2014
% Prepare data for job submission on the local Host machine
% Also run MP decomposition and save gaborInfo in the respective mpAnalysis
% folder for ready reference during analysis
%==========================================================================

clear;clc;close all;

monkeyName = 'murty';
expDate = '180914';
protocolName = 'GRF_002';
channelNumbers = 1:19;
folderSourceString = '/media/store/';
gridType = 'EEG';

%%
disp(['Number of electrodes/channels: ' num2str(length(channelNumbers))]);
disp('Preparing data....')
prepareDataForHost(monkeyName,expDate,protocolName,channelNumbers,folderSourceString,gridType);


%% Run MP decomposition on the local machine

disp('Running MP decomposition....');
for i = 1:length(channelNumbers)
    if isunix
        localCtlFile = ([folderSourceString 'data/' monkeyName '/' gridType '/' expDate '/' protocolName ...
            '/mpAnalysis/elec' num2str(channelNumbers(i)) '/ImportData_SIG/GaborMP/local.ctl']);
    else
        localCtlFile = ([folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName ...
            '\mpAnalysis\elec' num2str(channelNumbers(i)) '\ImportData_SIG\GaborMP\local.ctl']);
    end
    if ~exist(localCtlFile,'file')
        disp([localCtlFile ' does not exist']);
    else
        disp(['electrode/channel number:' num2str(channelNumbers(i))]);
        runMPDecomp(localCtlFile); % simply runs gabord on the above local.ctl
    end
end

%% Save gaborInfo for each electrode

disp('Saving gaborInfo....');
saveGaborInfo(monkeyName,expDate,protocolName,folderSourceString,gridType,channelNumbers);

%% Reconstruct energy spectrum for each trial and store the energy matrix
% Vinay - not a good startegy this! Creates huge files!
% 
% for i = 1:length(channelNumbers)
%     if ispc
%         mpFolder = ([folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName ...
%                 '\mpAnalysis\elec' num2str(channelNumbers(i)) '\']);
%         load([mpFolder '\gaborInfo.mat']);
%     else
%         mpFolder = ([folderSourceString 'data/' monkeyName '/' gridType '/' expDate '/' protocolName ...
%                 '/mpAnalysis/elec' num2str(channelNumbers(i)) '/']);
%         load([mpFolder '/gaborInfo.mat']);
%     end
%     
%     goodPos = 1:length(gaborInfo);
%     
% %     gaborInfoGoodPos = gaborInfo(goodPos);
% 
%     L = 4096;
%     wrap = [];
%     numAtomsMP = 100;
%     atomList = (1:numAtomsMP);
%     
%     mpEnergy = [];
%     mpEfileName = ([mpFolder 'mpEnergy.mat']);
%     save(mpEfileName,'mpEnergy');
%     mpE = matfile(mpEfileName,'Writable',true);
%     
%     disp(['Reconstructing Energy from:' num2str(numAtomsMP) 'atoms, and'  num2str(length(goodPos)) 'trials']);
%     disp(['Saving mpEnergy for electrode ' num2str(channelNumbers(i))]);
%     for m=1:length(goodPos)
%         disp(['trial number: ' num2str(m) '(actual trial - )' num2str(goodPos(m))]);
%         mpEnergy = reconstructEnergyFromAtomsMPP(gaborInfo{m}.gaborData,L,wrap,atomList);
%         mpE.mpEnergy(m,1:size(mpEnergy,1),1:size(mpEnergy,2)) = mpEnergy;
%         clear mpEnergy;
%     end
%     
%     
% %     save([mpFolder 'mpEnergy.mat'], 'mpEnergy');
%     clear gaborInfo
% end