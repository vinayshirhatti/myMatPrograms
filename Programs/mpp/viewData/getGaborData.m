% gaborInfo = getGaborData(folderName,tag,channelNum)
% This program reads the information about gabor atoms from the files
% generated by MP.

% Inputs
% folderName and tag: Standard inputs that specify the data
% channelNum: channel number of interest

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray, 2008 
% Distributed under the General Public License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add option for reading in windows framework

function gaborInfo = getGaborData(folderName,tag,channelNum,isWin)

if ~exist('isWin','var')            isWin=0;                        end

if isWin
    folderName=appendIfNotPresent(folderName,'\');
    tag=appendIfNotPresent(tag,'\');
    fn = [folderName tag 'GaborMP\mp0.bok.' conv2Str(channelNum-1)];
else
    folderName=appendIfNotPresent(folderName,'/');
    tag=appendIfNotPresent(tag,'/');
    fn = [folderName tag 'GaborMP/mp0.bok.' conv2Str(channelNum-1)];
end

if isWin
    fn = convertPlatform(fn,'win');
end
fd = fopen(fn,'r');

numAtoms=0;
trialNum=1;

while (~isempty(numAtoms))
    
    junk=fread(fd,1,'int'); % nvars
    junk=fread(fd,1,'int'); % window numb
    numAtoms=fread(fd,1,'int'); %no of atoms
    sigE=fread(fd,1,'double');
    sigE=sigE^2;
    
    if (~isempty(numAtoms))       
        gaborInfo{trialNum}.numAtoms    = numAtoms;
        gaborInfo{trialNum}.sigE        = sigE;
        
        x=zeros(5,numAtoms);
        for i=1:numAtoms
            clear k;
            k=fread(fd,4,'short');
            x(1:3,i)=k(2:4);
            x(4:5,i)=fread(fd,2,'double');
        end
        gaborInfo{trialNum}.gaborData = x;
        clear x;
        trialNum=trialNum+1;
    end
end

fclose(fd);
end