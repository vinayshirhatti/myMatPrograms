% writes the script file that must be run on the tyrone.
% modified from writeScriptFile for orchestra

function writeScriptFileTyrone(monkeyName,expDate,protocolName,channelNumbers,folderSourceString,gridType)

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
% MPSourceFolder = '/home/sr122/MP/source/';

% Vinay - the modified lines are ahead

folderSourceString = appendIfNotPresent(folderSourceString,'/');
folderNameMain = [folderSourceString 'data/' monkeyName '/' gridType '/' expDate '/' protocolName '/'];

% Input folder
inputFolder  = [folderNameMain 'segmentedData/LFP/'];

% Output folder
outputFolder = [folderNameMain 'mpAnalysis/'];
makeDirectory(outputFolder);

% This is where the data will be visible on Tyrone
linuxFolder  = ['/localscratch/apmsvd/MP/data/' monkeyName '/' gridType '/' expDate '/' protocolName '/mpAnalysis/'];
MPSourceFolder = '/localscratch/apmsvd/MP/source/';

% 
% %%%%%%%% Now generate the script file %%%%%%%%
% disp('Writing script file...');
% scriptFile = [outputFolder 'commandFile.sh'];
% fid = fopen(scriptFile,'wt');
% fprintf(fid, '#!/bin/bash\n');
% 
% for i=1:length(channelNumbers)
%     analogMatFile = [inputFolder 'elec' num2str(channelNumbers(i)) '.mat'];
% 
%     if ~exist(analogMatFile,'file')
%         disp([analogMatFile ' does not exist']);
%     else
%         tag = ['elec' num2str(channelNumbers(i)) '/'];
%         outputFile = [linuxFolder tag 'GaborMP/orchestraOutput'];
%         ctlFile = [linuxFolder tag 'ImportData_SIG/GaborMP/local.ctl'];
%         fprintf(fid,['bsub -o ' outputFile ' -r ' MPSourceFolder 'gabord ' ctlFile '\n']);
%     end
% end
% 
% fprintf(fid,'\n');
% fprintf(fid,'exit 0;');
% fclose(fid);

% Vinay - modified lines ahead

%%%%%%%% Now generate the script file %%%%%%%%
disp('Writing script file...');
scriptFile = [outputFolder 'commandFile.sh'];
fid = fopen(scriptFile,'wt');
fprintf(fid, '#!/bin/sh\n');
fprintf(fid, '#PBS -N mpjob\n');
fprintf(fid, '#PBS -l nodes=4:ppn=32:regular\n');
fprintf(fid, '#PBS -l walltime=24:00:00\n');
fprintf(fid, '\n');
fprintf(fid, 'cd /localscratch/apmsvd/MP/source/\n');
fprintf(fid, '\n');
fprintf(fid, 'NPROCS=`wc -l < $PBS_NODEFILE`\n');
fprintf(fid, 'HOSTS=`cat $PBS_NODEFILE | uniq | tr ''\\n'' "," | sed ''s|,$||''`\n'); % To print an apostrophe here we need to specify it as two consecutive apostrophes i.e. ''

for i=1:length(channelNumbers)
    analogMatFile = [inputFolder 'elec' num2str(channelNumbers(i)) '.mat'];

    if ~exist(analogMatFile,'file')
        disp([analogMatFile ' does not exist']);
    else
        tag = ['elec' num2str(channelNumbers(i)) '/'];
%         outputFile = [linuxFolder tag 'GaborMP/orchestraOutput'];
        outputFile = [linuxFolder tag 'GaborMP/tyroneOutput']; % Vinay - changed orchestra to tyrone
        ctlFile = [linuxFolder tag 'ImportData_SIG/GaborMP/local.ctl'];
%         fprintf(fid,['bsub -o ' outputFile ' -r ' MPSourceFolder 'gabord ' ctlFile '\n']);
        fprintf(fid,['mpirun -np $NPROCS --host $HOSTS gabord ' ctlFile '>' outputFile '\n']);
    end
end

fprintf(fid,'\n');
% fprintf(fid,'exit 0;'); % Vinay not required here
fclose(fid);

end