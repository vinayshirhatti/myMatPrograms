% Generate gaborInfo and save for each electrode
% Vinay Shirhatti, 24 August 2014
%--------------------------------------------------------------------------

function saveGaborInfo(monkeyName,expDate,protocolName,folderSourceString,gridType,channelNumbers)

if isunix
    folderName = [folderSourceString 'data/' monkeyName '/' gridType '/' ...
        expDate '/' protocolName '/mpAnalysis/'];
    slash = '/';
else
    folderName = [folderSourceString 'data\' monkeyName '\' gridType '\' ...
        expDate '\' protocolName '\mpAnalysis\'];
    slash = '\';
end

for i = 1:length(channelNumbers)
    tag{i} = ['elec' num2str(channelNumbers(i)) slash];
end


for i = 1:length(tag)
    elecName = [folderName tag{i}];
    if ~exist(elecName,'file')
        disp(['Electrode/Channel ' channelNumbers(i) ' does not exist']);
    else
        gaborInfo = getGaborData(folderName,tag{i},1);
        disp(['Saving gaborInfo for electrode ' num2str(channelNumbers(i))]);
        save([folderName tag{i} 'gaborInfo.mat'], 'gaborInfo');
    end
end

end
