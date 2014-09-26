% Bad stimList. This is specific for each monkey and grid type because of
% differences in the good electrodes and thresholds.

% monkeyName='abu'; gridType='ECoG'; folderSourceString='H:\';
% [expDates,protocolNames] = allProtocolsAbuECoG;

% monkeyName='sid'; gridType='Microelectrode'; folderSourceString='J:\';
% [expDates,protocolNames] = allProtocolsTestMicroelectrode;

monkeyName='murty'; gridType='EEG'; 
if ispc
    folderSourceString='J:\';
else
    folderSourceString='/media/store/';
end
[expDates,protocolNames] = allProtocolsTestMicroelectrode;


%findBadTrialsForTheseIndices=4:7; checkTheseElectrodes=[1 8 32 33 72 96]; showElectrodes=33; threshold=6; maxLimit=600; %022211

%findBadTrialsForTheseIndices=8:11; checkTheseElectrodes=[1 8 32 33]; showElectrodes=33; threshold=6; maxLimit=400; %022811
%findBadTrialsForTheseIndices=12:17; checkTheseElectrodes=[1 8 32 33]; showElectrodes=33; threshold=6; maxLimit=400; %030111
%findBadTrialsForTheseIndices=18:23; checkTheseElectrodes=[1 8 32 33]; showElectrodes=33; threshold=6; maxLimit=400; %030211
%findBadTrialsForTheseIndices=24:27; checkTheseElectrodes=[1 8 32 33]; showElectrodes=33; threshold=6; maxLimit=400; %030211
%findBadTrialsForTheseIndices=28:33; checkTheseElectrodes=[1 8 32]; showElectrodes=33; threshold=6; maxLimit=400; minLimit=-800; %030211
%findBadTrialsForTheseIndices=34:36; checkTheseElectrodes=[1 8 32]; showElectrodes=33; threshold=6; maxLimit=400; minLimit=-600; %030211
%findBadTrialsForTheseIndices=37:39; checkTheseElectrodes=[1 8 32]; showElectrodes=33; threshold=6; maxLimit=400; minLimit=-600; %030211
%findBadTrialsForTheseIndices=40:42; checkTheseElectrodes=[1 8 32]; showElectrodes=33; threshold=6; maxLimit=400; minLimit=-500; %030211
%findBadTrialsForTheseIndices=43:45; checkTheseElectrodes=[1 8 32]; showElectrodes=33; threshold=6; maxLimit=400; minLimit=-500; %030211
%findBadTrialsForTheseIndices=46:49; checkTheseElectrodes=[1 8 32]; showElectrodes=33; threshold=6; maxLimit=400; minLimit=-450; %030211
%findBadTrialsForTheseIndices=50:53; checkTheseElectrodes=[1 8 67]; showElectrodes=67; threshold=6; maxLimit=400; minLimit=-400; %030211
%findBadTrialsForTheseIndices=54:56; checkTheseElectrodes=[1 8 29]; showElectrodes=29; threshold=6; maxLimit=400; minLimit=-400; %030211
%findBadTrialsForTheseIndices=57:60; checkTheseElectrodes=[1 8 16]; showElectrodes=16; threshold=6; maxLimit=300; minLimit=-400; %030211
%findBadTrialsForTheseIndices=61:63; checkTheseElectrodes=[13]; showElectrodes=13; threshold=6; maxLimit=300; minLimit=-400; %030211
%findBadTrialsForTheseIndices=64:67; checkTheseElectrodes=[27]; showElectrodes=27; threshold=6; maxLimit=400; minLimit=-400; %030211
%findBadTrialsForTheseIndices=68:71; checkTheseElectrodes=[71]; showElectrodes=71; threshold=6; maxLimit=400; minLimit=-400; %030211
%findBadTrialsForTheseIndices=72:76; checkTheseElectrodes=[28]; showElectrodes=28; threshold=6; maxLimit=400; minLimit=-400; %030211
%findBadTrialsForTheseIndices=77:81; checkTheseElectrodes=[6]; showElectrodes=6; threshold=6; maxLimit=400; minLimit=-400; %030211
% findBadTrialsForTheseIndices=82:84; checkTheseElectrodes=[9]; showElectrodes=9; threshold=6; maxLimit=400; minLimit=-400; %030211

% findBadTrialsForTheseIndices=28; checkTheseElectrodes=[65 66 69 70 89 90 93 94]; showElectrodes=[65 66 69 70 89 90 93 94]; threshold=6; maxLimit=100; minLimit=-100; %200614

% findBadTrialsForTheseIndices=40; checkTheseElectrodes=[65 66 69 70]; showElectrodes=[65 66 69 70]; threshold=6; maxLimit=100; minLimit=-100; %060814

% findBadTrialsForTheseIndices=41; checkTheseElectrodes=[65 66 69 70]; showElectrodes=[65 66 69 70]; threshold=6; maxLimit=100; minLimit=-100; %070814

% findBadTrialsForTheseIndices=42; checkTheseElectrodes=[65 66 69 70]; showElectrodes=[65 66 69 70]; threshold=6; maxLimit=100; minLimit=-100; %070814

% findBadTrialsForTheseIndices=46; checkTheseElectrodes=[65 66 69 70 73 74 78]; showElectrodes=[65 66 69 70 73 74 78]; threshold=6; maxLimit=200; minLimit=-200; %130814

% findBadTrialsForTheseIndices=47:48; checkTheseElectrodes=1:10; showElectrodes=1:10; threshold=6; maxLimit=50; minLimit=-50; %080914

% findBadTrialsForTheseIndices=49; checkTheseElectrodes=1:10; showElectrodes=1:19; threshold=6; maxLimit=50; minLimit=-50; %080914

% findBadTrialsForTheseIndices=51; checkTheseElectrodes=1:10; showElectrodes=1:10; threshold=6; maxLimit=50; minLimit=-50; %090914

findBadTrialsForTheseIndices=55:56; checkTheseElectrodes=1:19; showElectrodes=1:19; threshold=6; maxLimit=50; minLimit=-50; %180914

% findBadTrialsForTheseIndices=57; checkTheseElectrodes=1:19; showElectrodes=1:19; threshold=6; maxLimit=50; minLimit=-50; %180914

% findBadTrialsForTheseIndices=58; checkTheseElectrodes=1:19; showElectrodes=1:19; threshold=6; maxLimit=30; minLimit=-30; %180914


for i=1:length(findBadTrialsForTheseIndices)
    expDate = expDates{findBadTrialsForTheseIndices(i)};
    protocolName = protocolNames{findBadTrialsForTheseIndices(i)};
    disp([expDate protocolName]);
    [allBadTrials,badTrials] = findBadTrialsWithLFP(monkeyName,expDate,protocolName,folderSourceString,gridType,checkTheseElectrodes,threshold,maxLimit,showElectrodes,minLimit);
    pause;
    clf;
end