% check whether MP has been done properly
function checkGaborComputation(monkeyName,expDate,protocolName,electrodeNumbers,gridType)

numElectrodes = length(electrodeNumbers);
success = zeros(1,numElectrodes);

for i=1:numElectrodes
    eNum = electrodeNumbers(i);
    
    success(i)=0;
    orchestraOutputFile = ['N:\maunsell\data\' monkeyName '\' gridType '\' expDate '\' ...
            protocolName '\mpAnalysis\elec' num2str(eNum) '\GaborMP\orchestraOutput']; 
    if exist(orchestraOutputFile,'file')
        output = textread(orchestraOutputFile,'%s'); 
        
        if isempty(output)
            disp([monkeyName expDate protocolName ', elec ' num2str(eNum) ' is empty']);
        elseif strcmpi(output{end},'Done')
            success(i)=1;
        else
            disp([monkeyName expDate protocolName ', elec ' num2str(eNum) ' not done']);
        end
    else
        disp([monkeyName expDate protocolName ', elec ' num2str(eNum) ' in progress']);
    end
end

if min(success)==1
    disp([monkeyName expDate protocolName ' done']);
end