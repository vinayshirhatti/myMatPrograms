% importData(X,folderName,tag,range,Fs,channelLabels,IDnum,firstName)

% Inputs
% X: A 3-D array of signals of size LxMxN, where L = signal length, M = Number of trials, N = number of Channels
% folderName: (default: pwd) is the folder where is directories are created 
% tag: An additional string to distinguish the dataset. Folders are created in foldername/tag. Default: 'test'
% signalRange: the portion of signal to be analysed. If not specified, the entire signal L is analyzed.
% Fs: Sampling Rate (default 1000)
% channelLabels: The Channel numbers of the N channels. Default: 1:N
% IDnum: The identity number of the set. Default: 1
% firstName: Used for patient data. Default: the input tag

function importData(X,folderName,tag,signalRange,Fs,channelLabels,IDnum,firstName)

sizeX = size(X);
if length(sizeX) == 2
    % just one Channel;
    L = sizeX(1); M = sizeX(2); N = 1;
else
    L = sizeX(1); M = sizeX(2); N = sizeX(3);
end

if ~exist('folderName','var')
	folderName = pwd;
end

if ~exist('tag','var')                tag = 'test/';                  end
if ~exist('signalRange','var')        signalRange = [1 L];            end
if ~exist('Fs','var')                 Fs = 1000;                      end
if ~exist('channelLabels','var')      channelLabels = 1:N;            end
if ~exist('IDnum','var')          	  IDnum = 1;                      end
if ~exist('firstName','var')      	  firstName = tag;                end 

writeMPfiles(X,folderName,tag,signalRange);
[EDF, goodChannels, numTrials] = getEDF(X,Fs,channelLabels);

signalLength = signalRange(2)- signalRange(1)+1;
writeheaderfile(folderName,tag,EDF,goodChannels, IDnum, firstName, numTrials,signalLength);

end