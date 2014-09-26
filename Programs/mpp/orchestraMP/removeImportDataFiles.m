% This program removes the ImportData files

function removeImportDataFiles(monkeyName,expDate,protocolName,folderSourceString,electrodeNumbers,gridType)

numElectrodes = length(electrodeNumbers);
success = zeros(1,numElectrodes);

for i=1:numElectrodes
    
    fileName = [folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName ...
        '\mpAnalysis\elec' num2str(electrodeNumbers(i)) '\ImportData_SIG\sig.dat.000'];
    if exist(fileName,'file');
        success(i)=1;
        delete(fileName);
    else
        %disp(['No ImportData in ' expDate protocolName 'elec' num2str(electrodeNumbers(i))]);
    end
end

disp([monkeyName expDate protocolName ': ' num2str(sum(success)) ' importData file(s) deleted']);
end