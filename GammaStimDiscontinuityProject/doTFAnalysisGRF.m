% Displays data from a single electrode, bipolar signals
% Pivots around any two selected parameters and shows the signal plots and 
% TFA plots for all their combinations 
% Vinay - modified from displaySingleChannelGRFv3.m for GRF protocol and
% doTFAnalysisCRS.m for CRS protocol
% 20 September 2014

function doTFAnalysisGRF(monkeyName,expDate,protocolName,folderSourceString,gridType)

if ~exist('folderSourceString','var')   folderSourceString='/media/store/';        end
if ~exist('gridType','var')             gridType='EEG';                             end

%---Vinay - set the paths based on the whether it is a Windows PC or a
%linux machine
if ispc
    folderName = [folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName '\'];

    % Get folders
    folderName = appendIfNotPresent(folderName,'\');
    folderExtract = [folderName 'extractedData\'];
    folderSegment = [folderName 'segmentedData\'];
    folderLFP = [folderSegment 'LFP\'];
    folderSpikes = [folderSegment 'Spikes\'];
end

if isunix
    folderName = [folderSourceString 'data/' monkeyName '/' gridType '/' expDate '/' protocolName '/'];

    % Get folders
    folderName = appendIfNotPresent(folderName,'/');
    folderExtract = [folderName 'extractedData/'];
    folderSegment = [folderName 'segmentedData/'];
    folderLFP = [folderSegment 'LFP/'];
    folderSpikes = [folderSegment 'Spikes/'];
end
%-----

% load LFP Information
[analogChannelsStored,timeVals,~,analogInputNums] = loadlfpInfo(folderLFP);
[neuralChannelsStored,SourceUnitIDs] = loadspikeInfo(folderSpikes);

% Get Combinations
[~,aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract);

% Get properties of the Stimulus
% stimResults = loadStimResults(folderExtract);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display main options
% fonts
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16; fontSizeTiny = 8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Panels
panelHeight = 0.34; panelStartHeight = 0.61;
staticPanelWidth = 0.25; staticStartPos = 0.025;
dynamicPanelWidth = 0.25; dynamicStartPos = 0.275;
timingPanelWidth = 0.25; timingStartPos = 0.525;
plotOptionsPanelWidth = 0.2; plotOptionsStartPos = 0.775;
backgroundColor = 'w';

% Vinay - define a flag for notching the line noise
notchData = 0;
% [Vinay] - define a flag to take bipolar (i.e. diff) signals: V1-V2
useBipolar = 0;
% [Vinay]- define a flag to plot SEM for the ERP plots
plotSEM = 1;
% [Vinay] - define a flag to save MP data if it is set
saveMPFlag = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Static Panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%staticTitle = [monkeyName '_' expDate '_' protocolName];
% if 0 % don't plot the static panel
%     hStaticPanel = uipanel('Title','Information','fontSize', fontSizeLarge, ...
%         'Unit','Normalized','Position',[staticStartPos panelStartHeight staticPanelWidth panelHeight]);
% 
%     staticText = [{ '   '};
%         {['Monkey Name: ' monkeyName]}; ...
%         {['Date: ' expDate]}; ...
%         {['Protocol Name: ' protocolName]}; ...
%         {'   '}
%         {['Orientation  (Deg): ' num2str(stimResults.orientation)]}; ...
%         {['Spatial Freq (CPD): ' num2str(stimResults.spatialFrequency)]}; ...
%         {['Eccentricity (Deg): ' num2str(stimResults.eccentricity)]}; ...
%         {['Polar angle  (Deg): ' num2str(stimResults.polarAngle)]}; ...
%         {['Sigma        (Deg): ' num2str(stimResults.sigma)]}; ...
%         {['Radius       (Deg): ' num2str(stimResults.radius)]}; ...
%         ];
% 
%     tStaticText = uicontrol('Parent',hStaticPanel,'Unit','Normalized', ...
%         'Position',[0 0 1 1], 'Style','text','String',staticText,'FontSize',fontSizeSmall);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dynamic panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynamicHeight = 0.06; dynamicGap=0.015; dynamicTextWidth = 0.6;
hDynamicPanel = uipanel('Title','Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[dynamicStartPos panelStartHeight dynamicPanelWidth panelHeight]);

% Analog channel
% Analog channel
[analogChannelStringList,analogChannelStringArray] = getAnalogStringFromValues(analogChannelsStored,analogInputNums);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Analog Channel (1:2)','FontSize',fontSizeTiny);
hAnalogChannel = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-(dynamicHeight+dynamicGap) (0.5*(1-dynamicTextWidth)) dynamicHeight], ...
    'Style','popup','String',analogChannelStringList,'FontSize',fontSizeTiny);
% [Vinay] - adding analog channel 2 for chosing bipolar signals V1-V2 
hAnalogChannel2 = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(dynamicTextWidth+(0.5*(1-dynamicTextWidth))) 1-(dynamicHeight+dynamicGap) (0.5*(1-dynamicTextWidth)) dynamicHeight], ...
    'Style','popup','String',(['none|' analogChannelStringList]),'FontSize',fontSizeTiny,'Callback',{@bipolar_Callback});

% Neural channel
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-2*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Neural Channel','FontSize',fontSizeTiny);
    
if ~isempty(neuralChannelsStored)
    neuralChannelString = getNeuralStringFromValues(neuralChannelsStored,SourceUnitIDs);
    hNeuralChannel = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-2*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',neuralChannelString,'FontSize',fontSizeTiny);
else
    hNeuralChannel = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position', [dynamicTextWidth 1-2*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','text','String','Not found','FontSize',fontSizeTiny);
end
% Sigma
sigmaString = getStringFromValues(sValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-3*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Sigma (Deg)','FontSize',fontSizeTiny);
hSigma = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-3*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',sigmaString,'FontSize',fontSizeTiny);

% Spatial Frequency
spatialFreqString = getStringFromValues(fValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-4*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Spatial Freq (CPD)','FontSize',fontSizeTiny);
hSpatialFreq = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-4*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',spatialFreqString,'FontSize',fontSizeTiny);

% Orientation
orientationString = getStringFromValues(oValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-5*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Orientation (Deg)','FontSize',fontSizeTiny);
hOrientation = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-5*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',orientationString,'FontSize',fontSizeTiny);

% Contrast
contrastString = getStringFromValues(cValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-6*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Contrast (%)','FontSize',fontSizeTiny);
hContrast = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-6*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',contrastString,'FontSize',fontSizeTiny);

% Temporal Frequency
temporalFreqString = getStringFromValues(tValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-7*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Temporal Freq (Hz)','FontSize',fontSizeTiny);
hTemporalFreq = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-7*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',temporalFreqString,'FontSize',fontSizeTiny);

% Azimuth
azimuthString = getStringFromValues(aValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-8*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Azimuth (Deg)','FontSize',fontSizeTiny);
hAzimuth = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-8*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',azimuthString,'FontSize',fontSizeTiny);

% Elevation
elevationString = getStringFromValues(eValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-9*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Elevation (Deg)','FontSize',fontSizeTiny);
hElevation = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position',...
    [dynamicTextWidth 1-9*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',elevationString,'FontSize',fontSizeTiny);


%Conjunction plot parameter and value. This parameter will be fixed to a
%value and the rest of the plots are taken by averaging across only this
%value of this parameter. Eg. All plots at contrast 50% of centre gabor

parametersString = 'NA|azimuth|elevation|sigma|spatialFreq|orientation|contrast|temporalFreq';

uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-10*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','pivot1','FontSize',fontSizeMedium);
% hGabor7 = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
%     'BackgroundColor', backgroundColor, 'Position', ...
%     [0 1-10*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
%     'Style','popup','String',gaborString,'FontSize',fontSizeTiny,'Callback',{@resetPivotParams_Callback});
hParam7 = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-10*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',parametersString,'FontSize',fontSizeTiny,'Callback',{@resetPivotParams_Callback});

%Conjunction plot parameter and value. This parameter will be fixed to a
%value and the rest of the plots are taken by averaging across only this
%value of this parameter. Eg. All plots at contrast 50% of centre gabor
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-12*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','pivot2','FontSize',fontSizeMedium);
% hGabor8 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
%     'BackgroundColor', backgroundColor, 'Position', ...
%     [textWidth 1-8*selectParamHeight selGWidth selectParamHeight], ...
%     'Style','popup','String',gaborString,'FontSize',fontSizeTiny,'Callback',{@resetPivotParams_Callback});
hParam8 = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-12*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',parametersString,'FontSize',fontSizeTiny,'Callback',{@resetPivotParams_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Timing panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timingHeight = 0.1; timingTextWidth = 0.5; timingBoxWidth = 0.25;
hTimingPanel = uipanel('Title','Timing','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[timingStartPos panelStartHeight timingPanelWidth panelHeight]);

signalRange = [-0.8 1.1];
fftRange = [0 250];
baseline = [-0.4 0];
stimPeriod = [0.4 0.8];

% Signal Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Parameter','FontSize',fontSizeMedium);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth 1-timingHeight timingBoxWidth timingHeight], ...
    'Style','text','String','Min','FontSize',fontSizeMedium);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth+timingBoxWidth 1-timingHeight timingBoxWidth timingHeight], ...
    'Style','text','String','Max','FontSize',fontSizeMedium);

% Stim Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-3*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim Range (s)','FontSize',fontSizeSmall);
hStimMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(1)),'FontSize',fontSizeSmall);
hStimMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(2)),'FontSize',fontSizeSmall);

% FFT Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-5*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','FFT Range (Hz)','FontSize',fontSizeSmall);
hFFTMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftRange(1)),'FontSize',fontSizeSmall);
hFFTMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftRange(2)),'FontSize',fontSizeSmall);

% Baseline
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-6*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Basline (s)','FontSize',fontSizeSmall);
hBaselineMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(1)),'FontSize',fontSizeSmall);
hBaselineMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(2)),'FontSize',fontSizeSmall);

% Stim Period
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-7*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim period (s)','FontSize',fontSizeSmall);
hStimPeriodMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-7*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(1)),'FontSize',fontSizeSmall);
hStimPeriodMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-7*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(2)),'FontSize',fontSizeSmall);

% Y Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-8*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Y Range','FontSize',fontSizeSmall);
hYMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-8*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','0','FontSize',fontSizeSmall);
hYMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-8*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','1','FontSize',fontSizeSmall);

% STA length
staLen = [-0.05 0.05]; 
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-9*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','STA len (s)','FontSize',fontSizeSmall);
hSTAMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-9*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(staLen(1)),'FontSize',fontSizeSmall);
hSTAMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-9*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(staLen(2)),'FontSize',fontSizeSmall);
hRemoveMeanSTA = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0 1-10*timingHeight 1 timingHeight], ...
    'Style','togglebutton','String','remove mean STA','FontSize',fontSizeMedium);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotOptionsHeight = 0.1;
hPlotOptionsPanel = uipanel('Title','Plotting Options','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[plotOptionsStartPos panelStartHeight plotOptionsPanelWidth panelHeight]);

% Button for Plotting
[colorString, colorNames] = getColorString;
uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 1-plotOptionsHeight 0.6 plotOptionsHeight], ...
    'Style','text','String','Color','FontSize',fontSizeTiny);

hChooseColor = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.6 1-plotOptionsHeight 0.4 plotOptionsHeight], ...
    'Style','popup','String',colorString,'FontSize',fontSizeSmall);

% Vinay - adding an option to adjust lineWidth of plots
linewidths = [1,1.5,2,2.5,3];
linewidthString = '1|1.5|2|2.5|3|';
uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 1-2*plotOptionsHeight 0.6 plotOptionsHeight], ...
    'Style','text','String','lineWidth','FontSize',fontSizeTiny);

hChooseWidth = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.6 1-2*plotOptionsHeight 0.4 plotOptionsHeight], ...
    'Style','popup','String',linewidthString,'FontSize',fontSizeTiny);

% Analysis Type
analysisTypeString = 'ERP|Firing Rate|Raster|FFT|delta FFT|STA';
uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 1-3*plotOptionsHeight 0.6 plotOptionsHeight], ...
    'Style','text','String','Analysis Type','FontSize',fontSizeSmall);
hAnalysisType = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.6 1-3*plotOptionsHeight 0.4 plotOptionsHeight], ...
    'Style','popup','String',analysisTypeString,'FontSize',fontSizeSmall);

% Vinay - adding a toggle button to use notched data
hNotchData = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 5*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','togglebutton','String','raw/lineNoiseRemoved','FontSize',fontSizeMedium, ...
    'Callback',{@notch_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 4*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','cla','FontSize',fontSizeMedium, ...
    'Callback',{@cla_Callback});

hHoldOn = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 3*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','togglebutton','String','hold on','FontSize',fontSizeMedium, ...
    'Callback',{@holdOn_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 2*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','rescale Y','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleY_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','rescale X','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleData_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 0 1 plotOptionsHeight], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium, ...
    'Callback',{@plotData_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get plots and message handles

% Get electrode array information
electrodeGridPos = [staticStartPos panelStartHeight staticPanelWidth panelHeight];
hElectrodes = showElectrodeLocations(electrodeGridPos,analogChannelsStored(get(hAnalogChannel,'val')), ...
    colorNames(get(hChooseColor,'val')),[],1,0,gridType);

% Set up the plot box and its dimensions
startXPos = staticStartPos; endXPos = 0.95; startYPos = 0.05; endYPos = 0.53;
centerGap = 0.05;
plotsWidthX = (endXPos-startXPos-centerGap);
plotsWidthY = (endYPos-startYPos-centerGap);
gap = 0.01; gapSmall = 0.002;
pivotGridPos=[startXPos startYPos endXPos endYPos];


uicontrol('Unit','Normalized','Position',[0 0.975 1 0.025],...
    'Style','text','String',[monkeyName expDate protocolName],'FontSize',fontSizeLarge);

param7 = get(hParam7, 'val');
param8 = get(hParam8, 'val');

hPivotParamPlot = getPlotHandles(length(getValsUnique(param7)),length(getValsUnique(param8)),pivotGridPos,gapSmall);


%**************************************************************************
% [Vinay] - TIME-FREQUENCY ANALYSIS
% Plot the time-frequency distributions on a separate figure
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hTFfig = figure(2);

% Define variables/flags required in TF analysis with their default values
plotStyle = 3; % pcolor, imagesc, line
spectrumType = 1; % raw, difference
cmin = -1; cmax = 3; % caxis limits
tfMethod = 1; % MTM, MP

% Default MTM params
mtmParams.Fs = 2000;
mtmParams.tapers=[1 1]; % [1 1] case is simply stft with dpss window
mtmParams.trialave=1;
mtmParams.err=0;
mtmParams.pad=0;

movingWin = [0.5 0.01];

% Default MP parameters
numAtomsMP = 100;


%__________________________________________________________________________
% The tf plots panel for the 6 chosen parameters
% ---------------tfPlotsPanel----------------------------------------------
%--------------------------------------------------------------------------

xGap = 0.01; yGap = 0.01;
tfPlotsHeight = 0.14; tfPlotsGap =0.02;
tfPanelWidth = 0.98; tfPanelHeight = 0.81;
hTFPlotsPanel = uipanel('Title','TF Plots','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[xGap yGap tfPanelWidth tfPanelHeight]);

startXPos = 2*xGap; startYPos = 2*yGap; % start from left bottom corner

pivotGridTFPos = [startXPos startYPos tfPanelWidth-2*xGap tfPanelHeight-5*yGap];

hPivotParamTFPlot = getPlotHandles(length(getValsUnique(param7)),length(getValsUnique(param8)),pivotGridTFPos,gapSmall);


%__________________________________________________________________________
% Parameters selection panels for TF methods
% ---------------tfParamsPanel----------------------------------------------
%--------------------------------------------------------------------------

xGap = 0.002;
tfParamsHeight = 1-(tfPanelHeight + yGap); tfParamsGap = tfPlotsGap;
tfParamsWidth = 0.4 - xGap; tfParamsPanelHeight = tfParamsHeight+tfPlotsGap;
tfParamsPanelxStart = (1-tfParamsWidth);
tfParamsPanelyStart = (1-tfParamsPanelHeight);

% hTFParamsPanel = uipanel('Title','TF Parameters','fontSize', fontSizeLarge, ...
%     'Unit','Normalized','Position',[tfParamsPanelxStart tfParamsPanelyStart tfParamsWidth tfParamsPanelHeight]);

textWidth = 0.24;
textHeight = 0.18;

% ====================TF Params Panel=====================================
hTFParamsPanel = uipanel('Title','MTM, MP, HHT Parameters','fontSize', fontSizeMedium, ...
    'Unit','Normalized','Position',[(tfParamsPanelxStart) (tfParamsPanelyStart) tfParamsWidth/2 tfParamsPanelHeight]);

% Parameters to divide the panel into rows and columns
xPanel1 = 0; xPanel2 = 0.5;

% Fs, sampling frequency
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-textHeight) textWidth textHeight], ...
    'Style','text','String','Fs','FontSize',fontSizeSmall);

hMTMFs = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-textHeight) textWidth textHeight], ...
    'Style','edit','String',mtmParams.Fs,'FontSize',fontSizeTiny,'Callback',{@resetMTMParams_Callback}); % initialize with default Fs = 2000

%------------------MTM parameters-----------------------
% Tapers TW
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-2*textHeight-tfParamsGap) textWidth textHeight], ...
    'Style','text','String','TW','FontSize',fontSizeSmall);

hMTMTapersTW = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-2*textHeight-tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',mtmParams.tapers(1),'FontSize',fontSizeTiny,'Callback',{@resetMTMParams_Callback}); % initialize TW = 1

% Tapers 'k'
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','text','String','k','FontSize',fontSizeSmall);

hMTMTapersK = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',mtmParams.tapers(2),'FontSize',fontSizeTiny,'Callback',{@resetMTMParams_Callback}); % initialize k = 1

% Window length
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-4*textHeight-3*tfParamsGap) textWidth textHeight], ...
    'Style','text','String','wLen','FontSize',fontSizeSmall);

hMTMwLen = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-4*textHeight-3*tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',movingWin(1),'FontSize',fontSizeTiny,'Callback',{@resetMTMParams_Callback}); % initialize wLen = 0.1 s

% Window translation step
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-5*textHeight-4*tfParamsGap) textWidth textHeight], ...
    'Style','text','String','wStep','FontSize',fontSizeSmall);
hMTMwStep = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-5*textHeight-4*tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',movingWin(2),'FontSize',fontSizeTiny,'Callback',{@resetMTMParams_Callback}); % initialize wStep = 0.01 s

%----MP parameters-----
% number of atoms for reconstruction
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel2 (1-textHeight) textWidth textHeight], ...
    'Style','text','String','numAtoms MP','FontSize',fontSizeSmall);
hMPnumAtoms = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel2+textWidth) (1-textHeight) textWidth textHeight], ...
    'Style','edit','String',numAtomsMP,'FontSize',fontSizeTiny,'Callback',{@resetMPParams_Callback}); % initialize numAtoms = 100


%========= Plotting settings ==============================
hTFPlotSettingsPanel = uipanel('Title','TF Plot settings','fontSize', fontSizeMedium, ...
    'Unit','Normalized','Position',[(tfParamsPanelxStart+(tfParamsWidth/2)+xGap) (tfParamsPanelyStart) ((tfParamsWidth/2)-xGap) tfParamsPanelHeight]);

% pcolor or imagesc or raw
uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-textHeight) 2*textWidth textHeight], ...
    'Style','text','String','Plot style','FontSize',fontSizeSmall);
hTFPlotStyle = uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+2*textWidth) (1-textHeight) 2*textWidth textHeight], ...
    'Style','pop','String','pcolor|imagesc|line','FontSize',fontSizeTiny,'Callback',{@resetTFSettings_Callback}); % initialize style = pcolor

% spectrum type: raw or difference (change from baseline)
uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-2*textHeight-tfParamsGap) 2*textWidth textHeight], ...
    'Style','text','String','Spectrum type','FontSize',fontSizeSmall);
hTFSpectrumType = uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+2*textWidth) (1-2*textHeight-tfParamsGap) 2*textWidth textHeight], ...
    'Style','pop','String','raw|difference','FontSize',fontSizeTiny,'Callback',{@resetTFSettings_Callback}); % initialize style = pcolor


% color axis settings
uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','text','String','cmin','FontSize',fontSizeSmall);
hTFcmin = uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',cmin,'FontSize',fontSizeTiny,'Callback',{@rescaleCaxis_Callback}); % initialize cmin = -1

uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[xPanel2 (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','text','String','cmax','FontSize',fontSizeSmall);
hTFcmax = uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel2+textWidth) (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',cmax,'FontSize',fontSizeTiny,'Callback',{@rescaleCaxis_Callback}); % initialize cmax = 3


% Method used
uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-5*textHeight-4*tfParamsGap) 2*textWidth 1.5*textHeight], ...
    'Style','text','String','Method','FontSize',fontSizeMedium);
hTFPlotMethod = uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+2*textWidth) (1-5*textHeight-4*tfParamsGap) 2*textWidth 1.5*textHeight], ...
    'Style','pop','String','MTM|MP','FontSize',fontSizeTiny,'Callback',{@resetTFSettings_Callback}); % initialize method = MTM


% --------- TF Plot callback button ---------------------------
uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[0 0 1 textHeight], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium, ...
    'Callback',{@plotTFData_Callback});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
    function plotData_Callback(~,~)
        a=get(hAzimuth,'val');
        e=get(hElevation,'val');
        s=get(hSigma,'val');
        f=get(hSpatialFreq,'val');
        o=get(hOrientation,'val');
        c=get(hContrast,'val');
        t=get(hTemporalFreq,'val');
        analysisType = get(hAnalysisType,'val');
        plotColor = colorNames(get(hChooseColor,'val'));
        BLMin = str2double(get(hBaselineMin,'String'));
        BLMax = str2double(get(hBaselineMax,'String'));
        STMin = str2double(get(hStimPeriodMin,'String'));
        STMax = str2double(get(hStimPeriodMax,'String'));
        STAMin = str2double(get(hSTAMin,'String'));
        STAMax = str2double(get(hSTAMax,'String'));
        holdOnState = get(hHoldOn,'val');
        removeMeanSTA = get(hRemoveMeanSTA,'val');
        notchData = get(hNotchData, 'val');
        
        plotLineWidth = linewidths(get(hChooseWidth,'val')); % Vinay - for plot line width
        
        hPivotParamPlot = getPlotHandles(length(getValsUnique(param7)),length(getValsUnique(param8)),pivotGridPos,gap);

        if analysisType==6 % Spike triggered average
            analogChannelPos = get(hAnalogChannel,'val');
            analogChannelString = analogChannelStringArray{analogChannelPos};
            spikeChannelPos = get(hNeuralChannel,'val');
            spikeChannelNumber = neuralChannelsStored(spikeChannelPos);
            unitID = SourceUnitIDs(spikeChannelPos);
            
            plotColors{1} = 'g';
            plotColors{2} = 'k';
            plotSTA1Channel(plotHandles,analogChannelString,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
                s,f,o,c,t,timeVals,plotColors,BLMin,BLMax,STMin,STMax,folderName,[STAMin STAMax],removeMeanSTA);
            
            % Write code for this
            %plotSTA1Parameter1Channel(hOrientationPlot,analogChannelString,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
            %    a,e,s,f,[],timeVals,plotColors,BLMin,BLMax,STMin,STMax,folderName);
            %plotSTA1Parameter1Channel(hSpatialFreqPlot,analogChannelString,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
            %    a,e,s,[],o,timeVals,plotColors,BLMin,BLMax,STMin,STMax,folderName);
            
            if analogChannelPos<=length(analogChannelsStored)
                analogChannelNumber = analogChannelsStored(analogChannelPos);
            else
                analogChannelNumber = 0;
            end
            channelNumber = [analogChannelNumber spikeChannelNumber];
            
        elseif analysisType == 2 || analysisType == 3
            channelPos = get(hNeuralChannel,'val');
            channelNumber = neuralChannelsStored(channelPos);
            unitID = SourceUnitIDs(channelPos);
            plotSpikeData1Channel(hPivotParamTFPlot,channelNumber,s,f,o,c,t,folderSpikes,...
                analysisType,timeVals,plotColor,unitID,folderName);
%             plotSpikeData1Parameter1Channel(hTemporalFreqPlot,channelNumber,a,e,s,f,o,c,[],folderSpikes,...
%                 analysisType,timeVals,plotColor,unitID,folderName);
%             plotSpikeData1Parameter1Channel(hContrastPlot,channelNumber,a,e,s,f,o,[],t,folderSpikes,...
%                 analysisType,timeVals,plotColor,unitID,folderName);
%             plotSpikeData1Parameter1Channel(hOrientationPlot,channelNumber,a,e,s,f,[],c,t,folderSpikes,...
%                 analysisType,timeVals,plotColor,unitID,folderName);
%             plotSpikeData1Parameter1Channel(hSpatialFreqPlot,channelNumber,a,e,s,[],o,c,t,folderSpikes,...
%                 analysisType,timeVals,plotColor,unitID,folderName);
%             plotSpikeData1Parameter1Channel(hSigmaPlot,channelNumber,a,e,[],f,o,c,t,folderSpikes,...
%                 analysisType,timeVals,plotColor,unitID,folderName);
        else
            analogChannelPos = get(hAnalogChannel,'val');
            analogChannelString = analogChannelStringArray{analogChannelPos};
            
            % [Vinay] - read the analog channel 2 string for bipolar case
            analogChannelPos2 = get(hAnalogChannel2,'val');
            analogChannelString2 = ('none');
            if analogChannelPos2 ~= 1
                analogChannelString2 = analogChannelStringArray{analogChannelPos2-1};
            end
            
%             rfMapVals = plotLFPData1Channel(plotHandles,analogChannelString,s,f,o,c,t,folderLFP,...
%                 analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, notchData); % Vinay - added notchData
%             plotLFPData1Parameter1Channel(hTemporalFreqPlot,analogChannelString,a,e,s,f,o,c,[],folderLFP,...
%                 analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, notchData); % Vinay - added notchData
%             plotLFPData1Parameter1Channel(hContrastPlot,analogChannelString,a,e,s,f,o,[],t,folderLFP,...
%                 analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, notchData); % Vinay - added notchData
%             plotLFPData1Parameter1Channel(hOrientationPlot,analogChannelString,a,e,s,f,[],c,t,folderLFP,...
%                 analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, notchData); % Vinay - added notchData
%             plotLFPData1Parameter1Channel(hSpatialFreqPlot,analogChannelString,a,e,s,[],o,c,t,folderLFP,...
%                 analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, notchData); % Vinay - added notchData
%             plotLFPData1Parameter1Channel(hSigmaPlot,analogChannelString,a,e,[],f,o,c,t,folderLFP,...
%                 analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, notchData); % Vinay - added notchData

            plotLFPDataPivotParameters1Channel(hPivotParamPlot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,param7, param8, folderLFP,...
                analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, notchData,useBipolar,plotSEM,holdOnState,plotLineWidth);

            if analogChannelPos<=length(analogChannelsStored)
                channelNumber = analogChannelsStored(analogChannelPos);
            else
                channelNumber = 0;
            end
            
        end

        if analysisType<=3  % ERP or spikes
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
        elseif analysisType <=5
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
        else
            xMin = str2double(get(hSTAMin,'String'));
            xMax = str2double(get(hSTAMax,'String'));
        end

%         rescaleData(plotHandles,xMin,xMax,getYLims(plotHandles));
%         rescaleData(hTemporalFreqPlot,xMin,xMax,getYLims(hTemporalFreqPlot));
%         rescaleData(hContrastPlot,xMin,xMax,getYLims(hContrastPlot));
%         rescaleData(hOrientationPlot,xMin,xMax,getYLims(hOrientationPlot));
%         rescaleData(hSpatialFreqPlot,xMin,xMax,getYLims(hSpatialFreqPlot));
%         rescaleData(hSigmaPlot,xMin,xMax,getYLims(hSigmaPlot));

        rescaleData(hPivotParamPlot,xMin,xMax,getYLims(hPivotParamPlot));
        
        showElectrodeLocations(electrodeGridPos,channelNumber,plotColor,hElectrodes,holdOnState,0,gridType);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleY_Callback(~,~)

        analysisType = get(hAnalysisType,'val');
        
        if analysisType<=3 % ERP or spikes
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
        else
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
        end

        yLims = [str2double(get(hYMin,'String')) str2double(get(hYMax,'String'))];
%         rescaleData(plotHandles,xMin,xMax,yLims);
%         rescaleData(hTemporalFreqPlot,xMin,xMax,yLims);
%         rescaleData(hContrastPlot,xMin,xMax,yLims);
%         rescaleData(hOrientationPlot,xMin,xMax,yLims);
%         rescaleData(hSpatialFreqPlot,xMin,xMax,yLims);
%         rescaleData(hSigmaPlot,xMin,xMax,yLims);
        
        rescaleData(hPivotParamPlot,xMin,xMax,yLims);
        
        % [Vinay] - rescale the TF plots
        figure(2);
        
        if plotStyle ~= 3 % pcolor, imagesc
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
            yLims = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
        else % line
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
            yLims = getYLims(hPivotParamTFPlot);
        end
        
        rescaleData(hPivotParamTFPlot,xMin,xMax,yLims);
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleData_Callback(~,~)

        analysisType = get(hAnalysisType,'val');

        if analysisType<=3 % ERP or spikes
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
        elseif analysisType==6
            xMin = str2double(get(hSTAMin,'String'));
            xMax = str2double(get(hSTAMax,'String'));
        else    
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
        end

%         rescaleData(plotHandles,xMin,xMax,getYLims(plotHandles));
%         rescaleData(hTemporalFreqPlot,xMin,xMax,getYLims(hTemporalFreqPlot));
%         rescaleData(hContrastPlot,xMin,xMax,getYLims(hContrastPlot));
%         rescaleData(hOrientationPlot,xMin,xMax,getYLims(hOrientationPlot));
%         rescaleData(hSpatialFreqPlot,xMin,xMax,getYLims(hSpatialFreqPlot));
%         rescaleData(hSigmaPlot,xMin,xMax,getYLims(hSigmaPlot));

        rescaleData(hPivotParamPlot,xMin,xMax,getYLims(hPivotParamPlot));
        
        % [Vinay] - rescale the TF plots
        figure(2);
        
        if plotStyle ~= 3 % pcolor, imagesc
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
            yLims = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
        else % line
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
            yLims = getYLims(hPivotParamTFPlot);
        end
        
        rescaleData(hPivotParamTFPlot,xMin,xMax,yLims);
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function holdOn_Callback(source,~)
        holdOnState = get(source,'Value');
        
%         holdOnGivenPlotHandle(plotHandles,holdOnState);
%         holdOnGivenPlotHandle(hTemporalFreqPlot,holdOnState);
%         holdOnGivenPlotHandle(hContrastPlot,holdOnState);
%         holdOnGivenPlotHandle(hOrientationPlot,holdOnState);
%         holdOnGivenPlotHandle(hSpatialFreqPlot,holdOnState);
%         holdOnGivenPlotHandle(hSigmaPlot,holdOnState);

        holdOnGivenPlotHandle(hPivotParamPlot,holdOnState);
        
        holdOnGivenPlotHandle(hPivotParamTFPlot,holdOnState);
        
        if holdOnState
            set(hElectrodes,'Nextplot','add');
        else
            set(hElectrodes,'Nextplot','replace');
        end

        function holdOnGivenPlotHandle(plotHandles,holdOnState)
            
            [numRows,numCols] = size(plotHandles);
            if holdOnState
                for i=1:numRows
                    for j=1:numCols
                        set(plotHandles(i,j),'Nextplot','add');

                    end
                end
            else
                for i=1:numRows
                    for j=1:numCols
                        set(plotHandles(i,j),'Nextplot','replace');
                    end
                end
            end
        end 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cla_Callback(~,~)
        
%         claGivenPlotHandle(plotHandles);
%         claGivenPlotHandle(hTemporalFreqPlot);
%         claGivenPlotHandle(hContrastPlot);
%         claGivenPlotHandle(hOrientationPlot);
%         claGivenPlotHandle(hSpatialFreqPlot);
%         claGivenPlotHandle(hSigmaPlot);
        
        claGivenPlotHandle(hPivotParamPlot);
        
%         cla(hRFMapPlot);cla(hcenterRFMapPlot);
        
        claGivenPlotHandle(hPivotParamTFPlot);
        
        function claGivenPlotHandle(plotHandles)
            [numRows,numCols] = size(plotHandles);
            for i=1:numRows
                for j=1:numCols
                    cla(plotHandles(i,j));
                end
            end
        end
    end
    
    function notch_Callback(hObject,~,~)
       notchData = get(hObject,'val');
    end

    function resetPivotParams_Callback(~,~)
%         adjFactor = 5;
%         gabor1 = adjFactor - get(hGabor1,'val'); % If sel is C, then value = 2 (NA = 1, C = 2, R = 3, S = 4) => gabor1 = 5 - 2 = 3. If R then 2, if S then 1
%         gabor2 = adjFactor - get(hGabor2,'val');
%         gabor3 = adjFactor - get(hGabor3,'val');
%         gabor4 = adjFactor - get(hGabor4,'val');
%         gabor5 = adjFactor - get(hGabor5,'val');
%         gabor6 = adjFactor - get(hGabor6,'val');
%         gabor7 = adjFactor - get(hGabor7,'val');
%         gabor8 = adjFactor - get(hGabor8,'val');
        
%         param1 = get(hParam1, 'val');
%         param2 = get(hParam2, 'val');
%         param3 = get(hParam3, 'val');
%         param4 = get(hParam4, 'val');
%         param5 = get(hParam5, 'val');
%         param6 = get(hParam6, 'val');
        param7 = get(hParam7, 'val');
        param8 = get(hParam8, 'val');
        
%         hParam1Plot = getPlotHandles(1,length(getValsUnique(gabor1, param1)),param1Grid,0.002);
%         hParam2Plot = getPlotHandles(1,length(getValsUnique(gabor2, param2)),param2Grid,0.002);
%         hParam3Plot = getPlotHandles(1,length(getValsUnique(gabor3, param3)),param3Grid,0.002);
%         hParam4Plot = getPlotHandles(1,length(getValsUnique(gabor4, param4)),param4Grid,0.002);
%         hParam5Plot = getPlotHandles(1,length(getValsUnique(gabor5, param5)),param5Grid,0.002);
%         hParam6Plot = getPlotHandles(1,length(getValsUnique(gabor6, param6)),param6Grid,0.002);
        
        hPivotParamPlot = getPlotHandles(length(getValsUnique(param7)),length(getValsUnique(param8)),pivotGridPos,gap);
        
        % [Vinay] - reset the TF grid as well
        figure(2);
%         hParam1TFPlot = getPlotHandles(1,length(getValsUnique(gabor1, param1)),param1TFGrid,0.002);
%         hParam2TFPlot = getPlotHandles(1,length(getValsUnique(gabor2, param2)),param2TFGrid,0.002);
%         hParam3TFPlot = getPlotHandles(1,length(getValsUnique(gabor3, param3)),param3TFGrid,0.002);
%         hParam4TFPlot = getPlotHandles(1,length(getValsUnique(gabor4, param4)),param4TFGrid,0.002);
%         hParam5TFPlot = getPlotHandles(1,length(getValsUnique(gabor5, param5)),param5TFGrid,0.002);
%         hParam6TFPlot = getPlotHandles(1,length(getValsUnique(gabor6, param6)),param6TFGrid,0.002);
        
        hPivotParamTFPlot = getPlotHandles(length(getValsUnique(param7)),length(getValsUnique(param8)),pivotGridTFPos,gap);
        
    end

    
    function valsUnique = getValsUnique(paramNumber)
        
        switch (paramNumber-1)
            case 1
                valsUnique = aValsUnique;
            case 2
                valsUnique = eValsUnique;
            case 3
                valsUnique = sValsUnique;
            case 4
                valsUnique = fValsUnique;
            case 5
                valsUnique = oValsUnique;
            case 6
                valsUnique = cValsUnique;
            case 7
                valsUnique = tValsUnique;
            otherwise
                valsUnique = [];
        end
    end

    
    function bipolar_Callback(hObject,~,~)
       if (get(hObject,'val') ~= 1)
           useBipolar = 1;
       else
           useBipolar = 0;
       end
    end

    %---------------TF analysis Functions----------------------------
    function resetMTMParams_Callback(~,~)
        mtmParams.Fs = str2double(get(hMTMFs,'String'));
        mtmParams.tapers(1) = str2double(get(hMTMTapersTW,'String'));
        mtmParams.tapers(2) = str2double(get(hMTMTapersK,'String'));
        
        movingWin(1) = str2double(get(hMTMwLen,'String'));
        movingWin(2) = str2double(get(hMTMwStep,'String'));
    end

    function resetMPParams_Callback(~,~)
        numAtomsMP = str2double(get(hMPnumAtoms,'String'));
    end

    function resetTFSettings_Callback(~,~)
        plotStyle = get(hTFPlotStyle,'val');
        spectrumType = get(hTFSpectrumType,'val');
        
        cmin = str2double(get(hTFcmin,'String'));
        cmax = str2double(get(hTFcmax,'String'));
        
        tfMethod = get(hTFPlotMethod,'val');
    end

    function rescaleCaxis_Callback(~,~)
        cmin = str2double(get(hTFcmin,'String'));
        cmax = str2double(get(hTFcmax,'String'));
        
%         set(hParam1TFPlot,'cLim',[cmin cmax]);
%         set(hParam2TFPlot,'cLim',[cmin cmax]);
%         set(hParam3TFPlot,'cLim',[cmin cmax]);
%         set(hParam4TFPlot,'cLim',[cmin cmax]);
%         set(hParam5TFPlot,'cLim',[cmin cmax]);
%         set(hParam6TFPlot,'cLim',[cmin cmax]);
        
        set(hPivotParamTFPlot,'cLim',[cmin cmax]);
        
%         for i = 1:size(tfplotHandles,1)
%             for j = 1:size(tfplotHandles,2)
%                 set(tfplotHandles(i,j),'cLim',[cmin cmax]);
%             end
%         end
        
    end
        
        
        

    %----main TF plotting function----------------
    function plotTFData_Callback(~,~)
        a=get(hAzimuth,'val');
        e=get(hElevation,'val');
        s=get(hSigma,'val');
        f=get(hSpatialFreq,'val');
        o=get(hOrientation,'val');
        c=get(hContrast,'val');
        t=get(hTemporalFreq,'val');
        
        analysisType = get(hAnalysisType,'val');
        plotColor = colorNames(get(hChooseColor,'val'));
        BLMin = str2double(get(hBaselineMin,'String'));
        BLMax = str2double(get(hBaselineMax,'String'));
        STMin = str2double(get(hStimPeriodMin,'String'));
        STMax = str2double(get(hStimPeriodMax,'String'));
        STAMin = str2double(get(hSTAMin,'String'));
        STAMax = str2double(get(hSTAMax,'String'));
        holdOnState = get(hHoldOn,'val');
        removeMeanSTA = get(hRemoveMeanSTA,'val');
        notchData = get(hNotchData, 'val');
        
        plotLineWidth = linewidths(get(hChooseWidth,'val')); % Vinay - for plot line width
        
        % Load TF Parameters
        % MTM
        mtmParams.Fs = str2double(get(hMTMFs,'String'));
        mtmParams.tapers(1) = str2double(get(hMTMTapersTW,'String'));
        mtmParams.tapers(2) = str2double(get(hMTMTapersK,'String'));
        movingWin(1) = str2double(get(hMTMwLen,'String'));
        movingWin(2) = str2double(get(hMTMwStep,'String'));
        
        % MP
        numAtomsMP = str2double(get(hMPnumAtoms,'String'));
        
        % Plot settings
        plotStyle = get(hTFPlotStyle,'val');
        spectrumType = get(hTFSpectrumType,'val');
        
        cmin = str2double(get(hTFcmin,'String'));
        cmax = str2double(get(hTFcmax,'String'));
        
        tfMethod = get(hTFPlotMethod,'val');
        
        
        % [Vinay] - get the plot handles here again
%         hParam1TFPlot = getPlotHandles(1,length(getValsUnique(gabor1, param1)),param1TFGrid,0.002);
%         hParam2TFPlot = getPlotHandles(1,length(getValsUnique(gabor2, param2)),param2TFGrid,0.002);
%         hParam3TFPlot = getPlotHandles(1,length(getValsUnique(gabor3, param3)),param3TFGrid,0.002);
%         hParam4TFPlot = getPlotHandles(1,length(getValsUnique(gabor4, param4)),param4TFGrid,0.002);
%         hParam5TFPlot = getPlotHandles(1,length(getValsUnique(gabor5, param5)),param5TFGrid,0.002);
%         hParam6TFPlot = getPlotHandles(1,length(getValsUnique(gabor6, param6)),param6TFGrid,0.002);
        
        hPivotParamTFPlot = getPlotHandles(length(getValsUnique(param7)),length(getValsUnique(param8)),pivotGridTFPos,gap);
        
        
        analogChannelPos = get(hAnalogChannel,'val');
        analogChannelString = analogChannelStringArray{analogChannelPos};

        % [Vinay] - read the analog channel 2 string for bipolar case
        analogChannelPos2 = get(hAnalogChannel2,'val');
        analogChannelString2 = ('none');
        if analogChannelPos2 ~= 1
            analogChannelString2 = analogChannelStringArray{analogChannelPos2-1};
        end
        
        
        tfplotLFPDataPivotParameters1Channel(hPivotParamTFPlot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,param7,param8,folderLFP,...
            timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, notchData,useBipolar,...
            tfMethod,mtmParams,movingWin,numAtomsMP,plotStyle,spectrumType,cmin,cmax, holdOnState, saveMPFlag,plotLineWidth);
        
    end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function that plots the data
% function rfMapVals = plotLFPData1Channel(plotHandles,channelString,s,f,o,c,t,folderLFP,...
% analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName,notchData) % Vinay - added notchData
% 
% if ispc
%     folderExtract = [folderName 'extractedData\'];
%     folderSegment = [folderName 'segmentedData\'];
% else
%     folderExtract = [folderName 'extractedData/'];
%     folderSegment = [folderName 'segmentedData/'];
% end
% 
% titleFontSize = 10;
% 
% [parameterCombinations,aValsUnique,eValsUnique] = loadParameterCombinations(folderExtract);
% [numRows,numCols] = size(plotHandles);
% 
% % Get the data
% removeAvgRef = 0;
% if removeAvgRef
%     disp('Removing average reference');
%     load([folderLFP 'avgRef']);
%     avgRef = analogData;
% end
% clear signal analogData
% load([folderLFP channelString]);
% if removeAvgRef
%     analogData = analogData-avgRef;
% end
% 
% % Get bad trials
% badTrialFile = [folderSegment 'badTrials.mat'];
% if ~exist(badTrialFile,'file')
%     disp('Bad trial file does not exist...');
%     badTrials=[];
% else
%     badTrials = loadBadTrials(badTrialFile);
%     disp([num2str(length(badTrials)) ' bad trials']);
% end
% 
% rfMapVals = zeros(numRows,numCols);
% for i=1:numRows
%     e = numRows-i+1;
%     for j=1:numCols
%         a = j;
%         clear goodPos
%         goodPos = parameterCombinations{a,e,s,f,o,c,t};
%         goodPos = setdiff(goodPos,badTrials);
%       
%         if isempty(goodPos)
%             disp('No entries for this combination..')
%         else
%             disp(['pos=(' num2str(i) ',' num2str(j) ') ,n=' num2str(length(goodPos))]);
%     
%             Fs = round(1/(timeVals(2)-timeVals(1)));
%             BLRange = uint16((BLMax-BLMin)*Fs);
%             STRange = uint16((STMax-STMin)*Fs);
%             BLPos = find(timeVals>=BLMin,1)+ (1:BLRange);
%             STPos = find(timeVals>=STMin,1)+ (1:STRange);
% 
%             xsBL = 0:1/(BLMax-BLMin):Fs-1/(BLMax-BLMin);
%             xsST = 0:1/(STMax-STMin):Fs-1/(STMax-STMin);
%             
% %             % Vinay - added next two 'if' conditions to take care of the case
% %             % when the lengths of BLPos(STPos) and xsBL(xsST) differ by 1
% %             if ((abs(length(BLRange)-length(xsBL)))==1)
% %                 BLRange = ceil((BLMax-BLMin)*Fs); % Vinay - added ceil so that the Range isn't a non-integer (else there's a mismatch in vector lengths ahead)
% %             end
% %         
% %             if ((abs(length(STRange)-length(xsST)))==1)
% %                 STRange = ceil((STMax-STMin)*Fs); % Vinay - added ceil so that the Range isn't a non-integer (else there's a mismatch in vector lengths ahead)
% %             end
%             
%             if notchData
%                 analogData = analogDataNotched;
%             end
% 
%             if analysisType == 1        % compute ERP
%                 clear erp
%                 erp = mean(analogData(goodPos,:),1); %#ok<*NODEF>
%                 plot(plotHandles(i,j),timeVals,erp,'color',plotColor);
%                 
%                 rfMapVals(e,a) = rms(erp(STPos));
% 
%             elseif analysisType == 2  ||   analysisType == 3 % compute Firing rates
%                 disp('Use plotSpikeData instead of plotLFPData...');
%             else
%                 
%                 fftBL = abs(fft(analogData(goodPos,BLPos),[],2));
%                 fftST = abs(fft(analogData(goodPos,STPos),[],2));
% 
%                 if analysisType == 4
%                     plot(plotHandles(i,j),xsBL,log10(mean(fftBL)),'g');
%                     set(plotHandles(i,j),'Nextplot','add');
%                     plot(plotHandles(i,j),xsST,log10(mean(fftST)),'k');
%                     set(plotHandles(i,j),'Nextplot','replace');
%                 end
% 
%                 if analysisType == 5
%                     if xsBL == xsST %#ok<BDSCI>
%                         plot(plotHandles(i,j),xsBL,log10(mean(fftST))-log10(mean(fftBL)),'color',plotColor);
%                     else
%                         disp('Choose same baseline and stimulus periods..');
%                     end
%                 end
%             end
%             
%             % Display title
%             if (i==1)
%                 if (j==1)
%                     title(plotHandles(i,j),['Azi: ' num2str(aValsUnique(a))],'FontSize',titleFontSize);
%                 else
%                     title(plotHandles(i,j),num2str(aValsUnique(a)),'FontSize',titleFontSize);
%                 end
%             end
%                 
%             if (j==numCols)
%                 if (i==1)
%                  title(plotHandles(i,j),[{'Ele'} {num2str(eValsUnique(e))}],'FontSize',titleFontSize,...
%                      'Units','Normalized','Position',[1.25 0.5]);
%                 else
%                     title(plotHandles(i,j),num2str(eValsUnique(e)),'FontSize',titleFontSize,...
%                      'Units','Normalized','Position',[1.25 0.5]);
%                 end
%             end
%         end
%     end
% end
% 
% if analysisType~=1
%     rfMapVals=[];
% end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotLFPDataPivotParameters1Channel(plotHandles,channelString,analogChannelString2,a,e,s,f,o,c,t,paramNum1,paramNum2,folderLFP,...
analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName,notchData,useBipolar,plotSEM,holdOnState,plotLineWidth) % Vinay - added notchData

if ispc
    folderExtract = [folderName 'extractedData\'];
    folderSegment = [folderName 'segmentedData\'];
else
    folderExtract = [folderName 'extractedData/'];
    folderSegment = [folderName 'segmentedData/'];
end

titleFontSize = 10;

timeForComputation = [40 100]/1000; % ms
freqForComputation = [40 60]; % Hz

[parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract);

numRows = size(plotHandles,1);
numCols = size(plotHandles,2);

% Get the data
clear signal analogData analogDataNotched
load([folderLFP channelString]);

% Get bad trials
badTrialFile = [folderSegment 'badTrials.mat'];
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end

% % Out of a,e,s,f,o,c and t only one parameter is empty
% if isempty(a)
%     aList = 1:length(aValsUnique); 
%     eList = e+zeros(1,numCols); 
%     sList = s+zeros(1,numCols);
%     fList = f+zeros(1,numCols); 
%     oList = o+zeros(1,numCols);
%     cList = c+zeros(1,numCols);
%     tList = t+zeros(1,numCols);
%     titleParam = 'Azi: ';
%     titleList = aValsUnique;
% end
% 
% if isempty(e) 
%     aList = a+zeros(1,numCols);
%     eList = 1:length(eValsUnique);
%     sList = s+zeros(1,numCols);
%     fList = f+zeros(1,numCols); 
%     oList = o+zeros(1,numCols);
%     cList = c+zeros(1,numCols);
%     tList = t+zeros(1,numCols);
%     titleParam = 'Ele: ';
%     titleList = eValsUnique;
% end
% 
% if isempty(s) 
%     aList = a+zeros(1,numCols); 
%     eList = e+zeros(1,numCols);
%     sList = 1:length(sValsUnique);
%     fList = f+zeros(1,numCols); 
%     oList = o+zeros(1,numCols);
%     cList = c+zeros(1,numCols);
%     tList = t+zeros(1,numCols);
%     titleParam = 'Sigma: ';
%     titleList = sValsUnique;
% end
% 
% if isempty(f)
%     aList = a+zeros(1,numCols); 
%     eList = e+zeros(1,numCols);
%     sList = s+zeros(1,numCols);
%     fList = 1:length(fValsUnique); 
%     oList = o+zeros(1,numCols);
%     cList = c+zeros(1,numCols);
%     tList = t+zeros(1,numCols);
%     titleParam = 'SF: ';
%     titleList = fValsUnique;
% end
% 
% if isempty(o) 
%     aList = a+zeros(1,numCols); 
%     eList = e+zeros(1,numCols);
%     sList = s+zeros(1,numCols);
%     fList = f+zeros(1,numCols); 
%     oList = 1:length(oValsUnique);
%     cList = c+zeros(1,numCols);
%     tList = t+zeros(1,numCols);
%     titleParam = 'Ori: ';
%     titleList = oValsUnique;
% end
% 
% if isempty(c) 
%     aList = a+zeros(1,numCols); 
%     eList = e+zeros(1,numCols);
%     sList = s+zeros(1,numCols);
%     fList = f+zeros(1,numCols); 
%     oList = o+zeros(1,numCols);
%     cList = 1:length(cValsUnique);
%     tList = t+zeros(1,numCols);
%     titleParam = 'Con: ';
%     titleList = cValsUnique;
% end
% 
% if isempty(t) 
%     aList = a+zeros(1,numCols); 
%     eList = e+zeros(1,numCols);
%     sList = s+zeros(1,numCols);
%     fList = f+zeros(1,numCols); 
%     oList = o+zeros(1,numCols);
%     cList = c+zeros(1,numCols);
%     tList = 1:length(tValsUnique);
%     titleParam = 'TF: ';
%     titleList = tValsUnique;
% end

% [Vinay] - repeat the set of parameters for each gabor depending on the
% number of columns to be drawn
aList = repmat(a,numRows,numCols);
eList = repmat(e,numRows,numCols);
sList = repmat(s,numRows,numCols);
fList = repmat(f,numRows,numCols);
oList = repmat(o,numRows,numCols);
cList = repmat(c,numRows,numCols);
tList = repmat(t,numRows,numCols);


    % [Vinay] decide the row parameter based on paramNum1
        for row = 1:numRows
            switch (paramNum1-1)
                case 1
                    aList(row,:) = row;                    
                    % Every new row takes a new value of param1.
                    % So assign the 'row' value in the list above
                    titleParam1 = 'Azi: ';
                    titleList1 = aValsUnique;
                case 2
                    eList(row,:) = row; 
                    titleParam1 = 'Ele: ';
                    titleList1 = eValsUnique;
                case 3
                    sList(row,:) = row; 
                    titleParam1 = 'Sigma: ';
                    titleList1 = sValsUnique;
                case 4
                    fList(row,:) = row; 
                    titleParam1 = 'SF: ';
                    titleList1 = fValsUnique;
                case 5
                    oList(row,:) = row; 
                    titleParam1 = 'Ori: ';
                    titleList1 = oValsUnique;
                case 6
                    cList(row,:) = row; 
                    titleParam1 = 'Contr: ';
                    titleList1 = cValsUnique;
                case 7
                    tList(row,:) = row; 
                    titleParam1 = 'TF: ';
                    titleList1 = tValsUnique;
                otherwise
                    disp('No particular gabor or parameter selected');
            end
        end
        
        % [Vinay] decide the column parameter based on paramNum2
        for row = 1:numRows
            switch (paramNum2-1)
                case 1
                    aList(row,:) = 1:length(aValsUnique);
                    % For every row we now have to assign indices
                    % incrementing along the columns.
                    % Basically the column entries go from 1 to 
                    % length{nValsUnique}
                    titleParam2 = 'Azi: ';
                    titleList2 = aValsUnique;
                case 2
                    eList(row,:) = 1:length(eValsUnique); 
                    titleParam2 = 'Ele: ';
                    titleList2 = eValsUnique;
                case 3
                    sList(row,:) = 1:length(sValsUnique);
                    titleParam2 = 'Sigma: ';
                    titleList2 = sValsUnique;
                case 4
                    fList(row,:) = 1:length(fValsUnique);
                    titleParam2 = 'SF: ';
                    titleList2 = fValsUnique;
                case 5
                    oList(row,:) = 1:length(oValsUnique);
                    titleParam2 = 'Ori: ';
                    titleList2 = oValsUnique;
                case 6
                    cList(row,:) = 1:length(cValsUnique);
                    titleParam2 = 'Contr: ';
                    titleList2 = cValsUnique;
                case 7
                    tList(row,:) = 1:length(tValsUnique);
                    titleParam2 = 'TF: ';
                    titleList2 = tValsUnique;
                otherwise
                    disp('No particular gabor or parameter selected');
            end
        end

        
% [Vinay] - Get the lengths of indices in parameterCombinations

    aLen = length(aValsUnique);
    eLen = length(eValsUnique);
    sLen = length(sValsUnique);
    fLen = length(fValsUnique);
    oLen = length(oValsUnique);
    cLen = length(cValsUnique);
    tLen = length(tValsUnique);
    
    % If more than one value, then length is one greater for all the values
    % together
    if (aLen> 1)           aLen=aLen+1;                    end
    if (eLen> 1)           eLen=eLen+1;                    end
    if (sLen> 1)           sLen=sLen+1;                    end
    if (fLen> 1)           fLen=fLen+1;                    end
    if (oLen> 1)           oLen=oLen+1;                    end
    if (cLen> 1)           cLen=cLen+1;                    end
    if (tLen> 1)           tLen=tLen+1;                    end


% Main loop
computationVals=zeros(1,numCols);
for k = 1:numRows
    for j=1:numCols
        clear goodPos
        goodPos = parameterCombinations{aList(k,j),eList(k,j),sList(k,j),fList(k,j),oList(k,j),cList(k,j),tList(k,j)};
        goodPos = setdiff(goodPos,badTrials);

        if isempty(goodPos)
            disp('No entries for this combination..')
        else
            disp(['pos=' num2str(k) ',' num2str(j) ',n=' num2str(length(goodPos))]);

            Fs = round(1/(timeVals(2)-timeVals(1)));
            BLRange = uint16((BLMax-BLMin)*Fs);
            STRange = uint16((STMax-STMin)*Fs);
            BLPos = find(timeVals>=BLMin,1)+ (1:BLRange);
            STPos = find(timeVals>=STMin,1)+ (1:STRange);

            xsBL = 0:1/(BLMax-BLMin):Fs-1/(BLMax-BLMin);
            xsST = 0:1/(STMax-STMin):Fs-1/(STMax-STMin);

    %         % Vinay - added next two 'if' conditions to take care of the case
    %         % when the lengths of BLPos(STPos) and xsBL(xsST) differ by 1
    %         if ((abs(length(BLRange)-length(xsBL)))==1)
    %             BLRange = ceil((BLMax-BLMin)*Fs); % Vinay - added ceil so that the Range isn't a non-integer (else there's a mismatch in vector lengths ahead)
    %         end
    %         
    %         if ((abs(length(STRange)-length(xsST)))==1)
    %             STRange = ceil((STMax-STMin)*Fs); % Vinay - added ceil so that the Range isn't a non-integer (else there's a mismatch in vector lengths ahead)
    %         end

            xsComputation = intersect(find(timeVals>=timeForComputation(1)),find(timeVals<timeForComputation(2)));
            freqComputation = intersect(find(xsST>=freqForComputation(1)),find(xsST<=freqForComputation(2)));

            % Vinay - added this notch data check
            if notchData
                analogData = analogDataNotched;
            end
            
            if useBipolar
                analogChannelString1 = channelString; % first electrode selected
                clear analogData analogDataNotched
                load([folderLFP analogChannelString1]);
                
                if notchData
                    analogData1 = analogDataNotched;
                else
                    analogData1 = analogData; % Vinay - these are the 
                    % freshly loaded values for electrode 1
                end
                % Vinay - store the analogData 
                % for electrode 1 in a separate variable 
                % otherwise loading the data for 2nd electrode will
                % overwrite analogData
                
                clear analogData analogDataNotched
                load([folderLFP analogChannelString2]);
                
                if notchData
                    analogData2 = analogDataNotched;
                else
                    analogData2 = analogData; % Vinay - these are the 
                    % freshly loaded values for electrode 2
                end
                analogData = analogData1 - analogData2;
            end

            if analysisType == 1        % compute ERP
                clear erp
                erp = mean(analogData(goodPos,:),1);
                
                plot(plotHandles(k,j),timeVals,erp,'color',plotColor,'Linewidth',plotLineWidth);
                
                % Vinay - plotting SEM
                if plotSEM
                    thisPlotColor = get(plot(plotHandles(k,j),timeVals,erp,'color',plotColor),'Color');
                    thisPlotColor(thisPlotColor==0) = 0.75;
                    erpSEM = std(analogData(goodPos,:),1)./sqrt(length(goodPos)); % Vinay - SEM = std/sqrt(n)
                    
                    set(plotHandles(k,j),'Nextplot','add');
                    plot(plotHandles(k,j),timeVals,(erp+erpSEM),'color',thisPlotColor);
                    plot(plotHandles(k,j),timeVals,(erp-erpSEM),'color',thisPlotColor);
                    plot(plotHandles(k,j),timeVals,erp,'color',plotColor,'Linewidth',plotLineWidth);
                    
                    if ~holdOnState
                        set(plotHandles(k,j),'Nextplot','replace');
                    end
                    
                end

%                 if isempty(o) % Orientation tuning
%                     computationVals(j) = abs(min(erp(xsComputation)));
%                 end

                if paramNum1==6 % Orientation tuning
                    computationVals(j) = abs(min(erp(xsComputation)));
                end
            
            elseif analysisType == 2 || analysisType == 3   % compute Firing rates
                disp('Use plotSpikeData instead of plotLFPData...');
            else

                fftBL = abs(fft(analogData(goodPos,BLPos),[],2));
                fftST = abs(fft(analogData(goodPos,STPos),[],2));

                if analysisType == 4
                    plot(plotHandles(k,j),xsBL,log10(mean(fftBL)),'g','Linewidth',plotLineWidth);
                    set(plotHandles(k,j),'Nextplot','add');
                    plot(plotHandles(k,j),xsST,log10(mean(fftST)),'k','Linewidth',plotLineWidth);
                    
                    if ~holdOnState
                        set(plotHandles(k,j),'Nextplot','replace');
                    end
                end

                if analysisType == 5
                    if xsBL == xsST %#ok<BDSCI>
                        plot(plotHandles(k,j),xsBL,log10(mean(fftST))-log10(mean(fftBL)),'color',plotColor,'Linewidth',plotLineWidth);
                    else
                        disp('Choose same baseline and stimulus periods..');
                    end
                end

%                 if isempty(o) % Orientation tuning
%                     computationVals(j) = max(mean(fftST(:,freqComputation),1));
%                 end
                if paramNum1==6 % Orientation tuning
                    computationVals(j) = max(mean(fftST(:,freqComputation),1));
                end
            end

            % Display title
            if (j==1 && k==1)
                title(plotHandles(k,j),[titleParam1 'vs' titleParam2 ':' num2str(titleList1(k)) ',' num2str(titleList2(j))],'FontSize',titleFontSize);
            else
                title(plotHandles(k,j),[num2str(titleList1(k)) ',' num2str(titleList2(j))],'FontSize',titleFontSize);
            end
        end
    end
end

% Orientation tuning
if paramNum1==6
    disp(['o: ' num2str(computationVals)]);
    [prefOrientation,orientationSelectivity] = getOrientationTuning(computationVals,oValsUnique);
    disp(['prefOri: ' num2str(round(prefOrientation)) ', sel: ' num2str(orientationSelectivity)]);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSpikeData1Channel(plotHandles,channelNumber,s,f,o,c,t,folderSpikes,...
analysisType,timeVals,plotColor,unitID,folderName)
titleFontSize = 12;

if ispc
    folderExtract = [folderName 'extractedData\'];
    folderSegment = [folderName 'segmentedData\'];
else
    folderExtract = [folderName 'extractedData/'];
    folderSegment = [folderName 'segmentedData/'];
end

[parameterCombinations,aValsUnique,eValsUnique] = loadParameterCombinations(folderExtract);
[numRows,numCols] = size(plotHandles);

% Get the data
clear signal spikeData
load([folderSpikes 'elec' num2str(channelNumber) '_SID' num2str(unitID)]);

% Get bad trials
badTrialFile = [folderSegment 'badTrials.mat'];
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end

for i=1:numRows
    e = numRows-i+1;
    for j=1:numCols
        a = j;
        clear goodPos
        goodPos = parameterCombinations{a,e,s,f,o,c,t};
        goodPos = setdiff(goodPos,badTrials);

        if isempty(goodPos)
            disp('No entries for this combination..')
        else
            disp(['pos=(' num2str(i) ',' num2str(j) ') ,n=' num2str(length(goodPos))]);
            
            if analysisType == 2
                [psthVals,xs] = getPSTH(spikeData(goodPos),10,timeVals(1),timeVals(end));
                plot(plotHandles(i,j),xs,psthVals,'color',plotColor);
            else
                X = spikeData(goodPos);
                axes(plotHandles(i,j)); %#ok<LAXES>
                rasterplot(X,1:length(X),plotColor);
            end
        end
        
        % Display title
        if (i==1)
            if (j==1)
                title(plotHandles(i,j),['Azi: ' num2str(aValsUnique(a))],'FontSize',titleFontSize);
            else
                title(plotHandles(i,j),num2str(aValsUnique(a)),'FontSize',titleFontSize);
            end
        end

        if (j==numCols)
            if (i==1)
                title(plotHandles(i,j),[{'Ele'} {num2str(eValsUnique(e))}],'FontSize',titleFontSize,...
                    'Units','Normalized','Position',[1.25 0.5]);
            else
                title(plotHandles(i,j),num2str(eValsUnique(e)),'FontSize',titleFontSize,...
                    'Units','Normalized','Position',[1.25 0.5]);
            end
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSpikeData1Parameter1Channel(plotHandles,channelNumber,a,e,s,f,o,c,t,folderSpikes,...
analysisType,timeVals,plotColor,unitID,folderName)
titleFontSize = 12;

if ispc
    folderExtract = [folderName 'extractedData\'];
    folderSegment = [folderName 'segmentedData\'];
else
    folderExtract = [folderName 'extractedData/'];
    folderSegment = [folderName 'segmentedData/'];
end

[parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract);
[~,numCols] = size(plotHandles);

% Get the data
clear signal spikeData
load([folderSpikes 'elec' num2str(channelNumber) '_SID' num2str(unitID)]);

% Get bad trials
badTrialFile = [folderSegment 'badTrials.mat'];
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end

% Out of a,e,s,f,o,c ant t, only one parameter is empty
if isempty(a)
    aList = 1:length(aValsUnique); 
    eList = e+zeros(1,numCols); 
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    cList = c+zeros(1,numCols); 
    tList = t+zeros(1,numCols);
    titleParam = 'Azi: ';
    titleList = aValsUnique;
end

if isempty(e) 
    aList = a+zeros(1,numCols);
    eList = 1:length(eValsUnique);
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    cList = c+zeros(1,numCols); 
    tList = t+zeros(1,numCols);
    titleParam = 'Ele: ';
    titleList = eValsUnique;
end

if isempty(s) 
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = 1:length(sValsUnique);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    cList = c+zeros(1,numCols); 
    tList = t+zeros(1,numCols);
    titleParam = 'Sigma: ';
    titleList = sValsUnique;
end

if isempty(f)
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = s+zeros(1,numCols);
    fList = 1:length(fValsUnique); 
    oList = o+zeros(1,numCols);
    cList = c+zeros(1,numCols); 
    tList = t+zeros(1,numCols);
    titleParam = 'SF: ';
    titleList = fValsUnique;
end

if isempty(o) 
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = 1:length(oValsUnique);
    cList = c+zeros(1,numCols); 
    tList = t+zeros(1,numCols);
    titleParam = 'Ori: ';
    titleList = oValsUnique;
end

if isempty(c) 
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    cList = 1:length(cValsUnique);
    tList = t+zeros(1,numCols);
    titleParam = 'Con: ';
    titleList = cValsUnique;
end

if isempty(t) 
    aList = a+zeros(1,numCols); 
    eList = e+zeros(1,numCols);
    sList = s+zeros(1,numCols);
    fList = f+zeros(1,numCols); 
    oList = o+zeros(1,numCols);
    cList = c+zeros(1,numCols);
    tList = 1:length(tValsUnique);
    titleParam = 'TF: ';
    titleList = tValsUnique;
end

% Plot

for j=1:numCols
    %a = j;
    clear goodPos
    goodPos = parameterCombinations{aList(j),eList(j),sList(j),fList(j),oList(j),cList(j),tList(j)};
    goodPos = setdiff(goodPos,badTrials);

    if isempty(goodPos)
        disp('No entries for this combination..')
    else
        disp(['pos=' num2str(j) ',n=' num2str(length(goodPos))]);
        if analysisType == 2
            [psthVals,xs] = getPSTH(spikeData(goodPos),10,timeVals(1),timeVals(end));
            plot(plotHandles(j),xs,psthVals,'color',plotColor);
        else
            X = spikeData(goodPos);
            axes(plotHandles(j)); %#ok<LAXES>
            rasterplot(X,1:length(X),plotColor);
        end
    end

    % Display title
    if (j==1)
        title(plotHandles(j),[titleParam num2str(titleList(j))],'FontSize',titleFontSize);
    else
        title(plotHandles(j),num2str(titleList(j)),'FontSize',titleFontSize);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSTA1Channel(plotHandles,analogChannelString,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
                s,f,o,c,t,timeVals,plotColors,BLMin,BLMax,STMin,STMax,folderName,staLen,removeMeanSTA)

titleFontSize = 12;

if ispc
    folderExtract = [folderName 'extractedData\'];
    folderSegment = [folderName 'segmentedData\'];
else
    folderExtract = [folderName 'extractedData/'];
    folderSegment = [folderName 'segmentedData/'];
end

[parameterCombinations,aValsUnique,eValsUnique] = loadParameterCombinations(folderExtract);
[numRows,numCols] = size(plotHandles);

% Get the analog data
clear signal analogData
load([folderLFP analogChannelString]);

% Get the spike data
clear signal spikeData
load([folderSpikes 'elec' num2str(spikeChannelNumber) '_SID' num2str(unitID)]);

% Get bad trials
badTrialFile = [folderSegment 'badTrials.mat'];
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end

staTimeLims{1} = [BLMin BLMax];
staTimeLims{2} = [STMin STMax];

for i=1:numRows
    e = numRows-i+1;
    for j=1:numCols
        a = j;
        clear goodPos
        goodPos = parameterCombinations{a,e,s,f,o,c,t};
        goodPos = setdiff(goodPos,badTrials);

        if isempty(goodPos)
            disp('No entries for this combination..')
        else
            goodSpikeData = spikeData(goodPos);
            goodAnalogSignal = analogData(goodPos,:);
            [staVals,numberOfSpikes,xsSTA] = getSTA(goodSpikeData,goodAnalogSignal,staTimeLims,timeVals,staLen,removeMeanSTA);
            
            disp([i j ', numStim: ' length(goodPos) ', numSpikes: ' num2str(numberOfSpikes)]);
            plot(plotHandles(i,j),xsSTA,staVals{1},'color',plotColors{1});
            set(plotHandles(i,j),'Nextplot','add');
            plot(plotHandles(i,j),xsSTA,staVals{2},'color',plotColors{2});
            set(plotHandles(i,j),'Nextplot','replace');
        end
        
        % Display title
        if (i==1)
            if (j==1)
                title(plotHandles(i,j),['Azi: ' num2str(aValsUnique(a))],'FontSize',titleFontSize);
            else
                title(plotHandles(i,j),num2str(aValsUnique(a)),'FontSize',titleFontSize);
            end
        end

        if (j==numCols)
            if (i==1)
                title(plotHandles(i,j),[{'Ele'} {num2str(eValsUnique(e))}],'FontSize',titleFontSize,...
                    'Units','Normalized','Position',[1.25 0.5]);
            else
                title(plotHandles(i,j),num2str(eValsUnique(e)),'FontSize',titleFontSize,...
                    'Units','Normalized','Position',[1.25 0.5]);
            end
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function plotSTA1Parameter1Channel(hOrientationPlot,analogChannelNumber,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
%                 a,e,s,f,o,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotRFMaps(hRFMapPlot,hcenterRFMapPlot,rfMapVals,aValsUnique,eValsUnique,plotColor,holdOnState)

plotSimple = 0;% Just compute the mean and Variance

if plotSimple
    [aziCenter,eleCenter] = getRFcenterSimple(aValsUnique,eValsUnique,rfMapVals); %#ok<*UNRCH>
else
    outParams = getRFcenter(aValsUnique,eValsUnique,rfMapVals);
    aziCenter = outParams(1); eleCenter = outParams(2);
    RFSize = sqrt((outParams(3)^2+outParams(4)^2)/2);
end

if plotSimple
    set(hRFMapPlot,'visible','off')
else
    % Plot the gaussian
    dX = (aValsUnique(end)-aValsUnique(1))/100;
    dY = (eValsUnique(end)-eValsUnique(1))/100;
    [~,outVals,boundaryX,boundaryY] = gauss2D(outParams,aValsUnique(1):dX:aValsUnique(end),eValsUnique(1):dY:eValsUnique(end));
    pcolor(hRFMapPlot,aValsUnique(1):dX:aValsUnique(end),eValsUnique(1):dY:eValsUnique(end),outVals);
    shading(hRFMapPlot,'interp'); 
    set(hRFMapPlot,'Nextplot','add');
    plot(hRFMapPlot,boundaryX,boundaryY,'k');
    set(hRFMapPlot,'Nextplot','replace');
end

% Plot the center only
if holdOnState
    set(hcenterRFMapPlot,'Nextplot','add');
else
    set(hcenterRFMapPlot,'Nextplot','replace');
end
plot(hcenterRFMapPlot,aziCenter,eleCenter,[plotColor '+']);
title(hcenterRFMapPlot,['Loc: (' num2str(aziCenter) ',' num2str(eleCenter) '), Size: ' num2str(RFSize)]);
axis(hcenterRFMapPlot,[aValsUnique(1) aValsUnique(end) eValsUnique(1) eValsUnique(end)]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yLims = getYLims(plotHandles)

[numRows,numCols] = size(plotHandles);
% Initialize
yMin = inf;
yMax = -inf;

for row=1:numRows
    for column=1:numCols
        % get positions
        axis(plotHandles(row,column),'tight');
        tmpAxisVals = axis(plotHandles(row,column));
        if tmpAxisVals(3) < yMin
            yMin = tmpAxisVals(3);
        end
        if tmpAxisVals(4) > yMax
            yMax = tmpAxisVals(4);
        end
    end
end

yLims=[yMin yMax];
end
function rescaleData(plotHandles,xMin,xMax,yLims)

[numRows,numCols] = size(plotHandles);
labelSize=12;
for i=1:numRows
    for j=1:numCols
        axis(plotHandles(i,j),[xMin xMax yLims]);
        if (i==numRows && rem(j,2)==1)
            if j~=1
                set(plotHandles(i,j),'YTickLabel',[],'fontSize',labelSize);
            end
        elseif (rem(i,2)==0 && j==1)
            set(plotHandles(i,j),'XTickLabel',[],'fontSize',labelSize);
        else
            set(plotHandles(i,j),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
        end
    end
end

% Remove Labels on the four corners
%set(plotHandles(1,1),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(1,numCols),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(numRows,1),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(numRows,numCols),'XTickLabel',[],'YTickLabel',[]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outString = getStringFromValues(valsUnique,decimationFactor)

if length(valsUnique)==1
    outString = convertNumToStr(valsUnique(1),decimationFactor);
else
    outString='';
    for i=1:length(valsUnique)
        outString = cat(2,outString,[convertNumToStr(valsUnique(i),decimationFactor) '|']);
    end
    outString = [outString 'all'];
end

    function str = convertNumToStr(num,f)
        if num > 16384
            num=num-32768;
        end
        str = num2str(num/f);
    end
end
function [outString,outArray] = getAnalogStringFromValues(analogChannelsStored,analogInputNums)
outString='';
count=1;
for i=1:length(analogChannelsStored)
    outArray{count} = ['elec' num2str(analogChannelsStored(i))]; %#ok<AGROW>
    outString = cat(2,outString,[outArray{count} '|']);
    count=count+1;
end
if ~isempty(analogInputNums)
    for i=1:length(analogInputNums)
        outArray{count} = ['ainp' num2str(analogInputNums(i))]; %#ok<AGROW>
        outString = cat(2,outString,[outArray{count} '|']);
        count=count+1;
    end
end
end
function outString = getNeuralStringFromValues(neuralChannelsStored,SourceUnitIDs)
outString='';
for i=1:length(neuralChannelsStored)
    outString = cat(2,outString,[num2str(neuralChannelsStored(i)) ', SID ' num2str(SourceUnitIDs(i)) '|']);
end 
end
function [colorString, colorNames] = getColorString

colorNames = 'brkgcmy';
colorString = 'blue|red|black|green|cyan|magenta|yellow';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%c%%%%%%%%%
%%%%%%%%%%%%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Data
function [analogChannelsStored,timeVals,goodStimPos,analogInputNums] = loadlfpInfo(folderLFP) %#ok<*STOUT>
load([folderLFP 'lfpInfo']);
if ~exist('analogInputNums','var')
    analogInputNums=[];
end
end
function [neuralChannelsStored,SourceUnitID] = loadspikeInfo(folderSpikes)
fileName = [folderSpikes 'spikeInfo.mat'];
if exist(fileName,'file')
    load(fileName);
else
    neuralChannelsStored=[];
    SourceUnitID=[];
end
end
function [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract)

load([folderExtract 'parameterCombinations.mat']);

if ~exist('sValsUnique','var')
    sValsUnique=rValsUnique;
end

if ~exist('cValsUnique','var')
    cValsUnique=[];
end

if ~exist('tValsUnique','var')
    tValsUnique=[];
end
end
% function stimResults = loadStimResults(folderExtract)
% load ([folderExtract 'stimResults']);
% end
function badTrials = loadBadTrials(badTrialFile)
load(badTrialFile);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [aziCenter,eleCenter] = getRFcenterSimple(aValsUnique,eValsUnique,rfMapValues)

rfMapValues = rfMapValues-mean(mean(rfMapValues));
% Compute the mean and variance
[maxEnonUnique,maxAnonUnique] = (find(rfMapValues==max(max(rfMapValues))));
maxE = maxEnonUnique(ceil(length(maxEnonUnique)/2));
maxA = maxAnonUnique(ceil(length(maxAnonUnique)/2));

%[maxE maxA]
aziCenter = sum(rfMapValues(maxE,:).*aValsUnique)/sum(rfMapValues(maxE,:));
eleCenter = sum(rfMapValues(:,maxA)'.*eValsUnique)/sum(rfMapValues(:,maxA));
end
function [prefOrientation,orientationSelectivity] = getOrientationTuning(computationVals,oValsUnique)
num=0;
den=0;

for j=1:length(oValsUnique)
    num = num+computationVals(j)*sind(2*oValsUnique(j));
    den = den+computationVals(j)*cosd(2*oValsUnique(j));
end

prefOrientation = 90*atan2(num,den)/pi;
orientationSelectivity = abs(den+1i*num)/sum(computationVals);

if prefOrientation<0
    prefOrientation = prefOrientation+180;
end
end


%--------- TF plots functions

function tfplotLFPDataPivotParameters1Channel(tfplotHandles,channelString,analogChannelString2,a,e,s,f,o,c,t,paramNum1,paramNum2,folderLFP,...
            timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName,notchData,useBipolar,...
            tfMethod,mtmParams,movingWin,numAtomsMP,plotStyle,spectrumType,cmin,cmax,holdOnState,saveMPFlag,plotLineWidth)
        
        if ispc
            folderExtract = [folderName 'extractedData\'];
            folderSegment = [folderName 'segmentedData\'];
            mpFolder = [folderName 'mpAnalysis\']; % for using stored MP data
            mpFolderTemp = 'C:\Documents\MP\data\'; % for online calculations
        else
            folderExtract = [folderName 'extractedData/'];
            folderSegment = [folderName 'segmentedData/'];
            mpFolder = [folderName 'mpAnalysis/']; % for using stored MP data
            mpFolderTemp = '/home/vinay/Documents/MP/data/'; % for online calculations
        end


        titleFontSize = 9;

%         timeForComputation = [40 100]/1000; % ms
%         freqForComputation = [40 60]; % Hz

        [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
            fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract);
        
        numRows = size(tfplotHandles,1);
        numCols = size(tfplotHandles,2);

        % Get the data
        clear signal analogData
        load([folderLFP channelString]);

        % Get bad trials
        badTrialFile = [folderSegment 'badTrials.mat'];
        if ~exist(badTrialFile,'file')
            disp('Bad trial file does not exist...');
            badTrials=[];
        else
            badTrials = loadBadTrials(badTrialFile);
            disp([num2str(length(badTrials)) ' bad trials']);
        end

        % [Vinay] - repeat the set of parameters for each gabor depending on the
        % number of columns to be drawn
        aList = repmat(a,numRows,numCols);
        eList = repmat(e,numRows,numCols);
        sList = repmat(s,numRows,numCols);
        fList = repmat(f,numRows,numCols);
        oList = repmat(o,numRows,numCols);
        cList = repmat(c,numRows,numCols);
        tList = repmat(t,numRows,numCols);


        % [Vinay] decide the row parameter based on paramNum1
        for row = 1:numRows
            switch (paramNum1-1)
                case 1
                    aList(row,:) = row;                    
                    % Every new row takes a new value of param1.
                    % So assign the 'row' value in the list above
                    titleParam1 = 'Azi: ';
                    titleList1 = aValsUnique;
                case 2
                    eList(row,:) = row; 
                    titleParam1 = 'Ele: ';
                    titleList1 = eValsUnique;
                case 3
                    sList(row,:) = row; 
                    titleParam1 = 'Sigma: ';
                    titleList1 = sValsUnique;
                case 4
                    fList(row,:) = row; 
                    titleParam1 = 'SF: ';
                    titleList1 = fValsUnique;
                case 5
                    oList(row,:) = row; 
                    titleParam1 = 'Ori: ';
                    titleList1 = oValsUnique;
                case 6
                    cList(row,:) = row; 
                    titleParam1 = 'Contr: ';
                    titleList1 = cValsUnique;
                case 7
                    tList(row,:) = row; 
                    titleParam1 = 'TF: ';
                    titleList1 = tValsUnique;
                otherwise
                    disp('No particular gabor or parameter selected');
            end
        end

        % [Vinay] decide the column parameter based on paramNum2
        for row = 1:numRows
            switch (paramNum2-1)
                case 1
                    aList(row,:) = 1:length(aValsUnique);
                    % For every row we now have to assign indices
                    % incrementing along the columns.
                    % Basically the column entries go from 1 to 
                    % length{nValsUnique}
                    titleParam2 = 'Azi: ';
                    titleList2 = aValsUnique;
                case 2
                    eList(row,:) = 1:length(eValsUnique); 
                    titleParam2 = 'Ele: ';
                    titleList2 = eValsUnique;
                case 3
                    sList(row,:) = 1:length(sValsUnique);
                    titleParam2 = 'Sigma: ';
                    titleList2 = sValsUnique;
                case 4
                    fList(row,:) = 1:length(fValsUnique);
                    titleParam2 = 'SF: ';
                    titleList2 = fValsUnique;
                case 5
                    oList(row,:) = 1:length(oValsUnique);
                    titleParam2 = 'Ori: ';
                    titleList2 = oValsUnique;
                case 6
                    cList(row,:) = 1:length(cValsUnique);
                    titleParam2 = 'Contr: ';
                    titleList2 = cValsUnique;
                case 7
                    tList(row,:) = 1:length(tValsUnique);
                    titleParam2 = 'TF: ';
                    titleList2 = tValsUnique;
                otherwise
                    disp('No particular gabor or parameter selected');
            end
        end


    % [Vinay] - Get the lengths of indices in parameterCombinations

        aLen = length(aValsUnique);
        eLen = length(eValsUnique);
        sLen = length(sValsUnique);
        fLen = length(fValsUnique);
        oLen = length(oValsUnique);
        cLen = length(cValsUnique);
        tLen = length(tValsUnique);

        % If more than one value, then length is one greater for all the values
        % together
        if (aLen> 1)           aLen=aLen+1;                    end
        if (eLen> 1)           eLen=eLen+1;                    end
        if (sLen> 1)           sLen=sLen+1;                    end
        if (fLen> 1)           fLen=fLen+1;                    end
        if (oLen> 1)           oLen=oLen+1;                    end
        if (cLen> 1)           cLen=cLen+1;                    end
        if (tLen> 1)           tLen=tLen+1;                    end


    % Main loop
    for k = 1:numRows
        for j=1:numCols
            clear goodPos
            goodPos = parameterCombinations{aList(k,j),eList(k,j),sList(k,j),fList(k,j),oList(k,j),cList(k,j),tList(k,j)};
            
            tagPos = [num2str(aList(k,j)) num2str(eList(k,j)) num2str(sList(k,j))...
                               num2str(fList(k,j)) num2str(oList(k,j)) num2str(cList(k,j))...
                               num2str(tList(k,j))];
            
            goodPos = setdiff(goodPos,badTrials);

            if isempty(goodPos)
                disp('No entries for this combination..')
            else
                disp(['pos=' num2str(k) ',' num2str(j) ',n=' num2str(length(goodPos))]);

                Fs = round(1/(timeVals(2)-timeVals(1)));
                BLRange = uint16((BLMax-BLMin)*Fs);
                STRange = uint16((STMax-STMin)*Fs);
                BLPos = find(timeVals>=BLMin,1)+ (1:BLRange);
                STPos = find(timeVals>=STMin,1)+ (1:STRange);

                xsBL = 0:1/(BLMax-BLMin):Fs-1/(BLMax-BLMin);
                xsST = 0:1/(STMax-STMin):Fs-1/(STMax-STMin);

        %         % Vinay - added next two 'if' conditions to take care of the case
        %         % when the lengths of BLPos(STPos) and xsBL(xsST) differ by 1
        %         if ((abs(length(BLRange)-length(xsBL)))==1)
        %             BLRange = ceil((BLMax-BLMin)*Fs); % Vinay - added ceil so that the Range isn't a non-integer (else there's a mismatch in vector lengths ahead)
        %         end
        %         
        %         if ((abs(length(STRange)-length(xsST)))==1)
        %             STRange = ceil((STMax-STMin)*Fs); % Vinay - added ceil so that the Range isn't a non-integer (else there's a mismatch in vector lengths ahead)
        %         end

%                 xsComputation = intersect(find(timeVals>=timeForComputation(1)),find(timeVals<timeForComputation(2)));
%                 freqComputation = intersect(find(xsST>=freqForComputation(1)),find(xsST<=freqForComputation(2)));

                % Vinay - added this notch data check
                if notchData
                    analogData = analogDataNotched;
                end

                if useBipolar
                    analogChannelString1 = channelString; % first electrode selected
                    clear analogData analogDataNotched
                    load([folderLFP analogChannelString1]);
                
                    if notchData
                        analogData1 = analogDataNotched;
                    else
                        analogData1 = analogData; % Vinay - these are the 
                        % freshly loaded values for electrode 1
                    end
                    % Vinay - store the analogData 
                    % for electrode 1 in a separate variable 
                    % otherwise loading the data for 2nd electrode will
                    % overwrite analogData
                    clear analogData analogDataNotched
                    load([folderLFP analogChannelString2]);
                    if notchData
                        analogData2 = analogDataNotched;
                    else
                        analogData2 = analogData; % Vinay - these are the 
                        % freshly loaded values for electrode 2
                    end
                    analogData = analogData1 - analogData2;
                end
                
                
                %--------------------------------------------------------------
                % Implement the method set by tfMethod varibale

                if (tfMethod == 1) % MTM
                    %----------------------------------------------------------
                    if (plotStyle == 3) % line plot (power vs f)
                        if (spectrumType == 1) % raw
                            [SBL,f]=mtspectrumc((analogData(goodPos,BLPos))',mtmParams); % baseline period
                            [SST,f]=mtspectrumc((analogData(goodPos,STPos))',mtmParams); % stimulus period
                            plot(tfplotHandles(k,j),f,conv2Log(SST),'color',plotColor,'Linewidth',plotLineWidth);
                            set(tfplotHandles(k,j),'Nextplot','add');
                            plot(tfplotHandles(k,j),f,conv2Log(SBL),'color','g','Linewidth',plotLineWidth);
                            if holdOnState
                                set(tfplotHandles(k,j),'Nextplot','add');
                            else
                                set(tfplotHandles(k,j),'Nextplot','replace');
                            end

                        elseif (spectrumType == 2) % difference from baseline
                            if xsST == xsBL 
                                [SBL,f]=mtspectrumc((analogData(goodPos,BLPos))',mtmParams);
                                [SST,f]=mtspectrumc((analogData(goodPos,STPos))',mtmParams);
                                plot(tfplotHandles(k,j),f,10*(conv2Log(SST)-conv2Log(SBL)),'color',plotColor,'Linewidth',plotLineWidth);
                            else
                                disp('Choose same baseline and stimulus periods..');
                            end
                            if holdOnState
                                set(tfplotHandles(k,j),'Nextplot','add');
                            else
                                set(tfplotHandles(k,j),'Nextplot','replace');
                            end
                        end

                    %----------------------------------------------------------    
                    elseif (plotStyle == 1) % pcolor plot

                        if (spectrumType == 1) % raw
                            [S1,t,f]=mtspecgramc((analogData(goodPos,:))',movingWin,mtmParams);
                            t = t + timeVals(1); % shift the t values to the actual time
                            pcolor(tfplotHandles(k,j),t,f,conv2Log(S1')); shading interp;
                            caxis([cmin cmax]);

                        elseif (spectrumType == 2) % difference from baseline
                            [S1,t,f]=mtspecgramc((analogData(goodPos,:))',movingWin,mtmParams);
                            t = t + timeVals(1); % shift the t values to the actual time

                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            S1BL = S1(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogS1BL = mean(conv2Log(S1BL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dS1 = 10*(conv2Log(S1) - repmat(mlogS1BL,size(S1,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            pcolor(tfplotHandles(k,j),t,f,dS1'); shading interp;
                            caxis([cmin cmax]);
                        end

                    %----------------------------------------------------------    
                    elseif (plotStyle == 2) % imagesc plot

                        if (spectrumType == 1) % raw
                            [S1,t,f]=mtspecgramc((analogData(goodPos,:))',movingWin,mtmParams);
                            t = t + timeVals(1); % shift the t values to the actual time
                            imagesc(t,f,conv2Log(S1'),'Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);

                        elseif (spectrumType == 2) % difference from baseline
                            [S1,t,f]=mtspecgramc((analogData(goodPos,:))',movingWin,mtmParams);
                            t = t + timeVals(1); % shift the t values to the actual time

                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            S1BL = S1(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogS1BL = mean(conv2Log(S1BL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dS1 = 10*(conv2Log(S1) - repmat(mlogS1BL,size(S1,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum                       
                            imagesc(t,f,dS1','Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                        end
                    end
                    
                    %======================================================

                elseif (tfMethod == 2) % MP Algorithm
                    
                    disp(['number of atoms for reconstruction:' num2str(numAtomsMP)]);
                    
                    % Vinay - create a directory (if not already present)
                    % to store the saved mp energy data
                    if isunix
                        folderSaveMP = [folderName 'mpSavedSpectrum/'];
                        makeDirectory(folderSaveMP);                        
                    else
                        folderSaveMP = [folderName 'mpSavedSpectrum\'];
                        makeDirectory(folderSaveMP);
                    end
                    
%                     tagPos = (['Ca' num2str(a(3)) 'e' num2str(e(3)) 's' num2str(s(3)) 'f' num2str(f(3)) 'o' num2str(o(3)) 'c' num2str(c(3)) 't' num2str(t(3)) 'p' num2str(p(3)) 'r' num2str(r(3)) ...
%                                    '_Ra' num2str(a(2)) 'e' num2str(e(2)) 's' num2str(s(2)) 'f' num2str(f(2)) 'o' num2str(o(2)) 'c' num2str(c(2)) 't' num2str(t(2)) 'p' num2str(p(2)) 'r' num2str(r(2)) ...
%                                    '_Sa' num2str(a(1)) 'e' num2str(e(1)) 's' num2str(s(1)) 'f' num2str(f(1)) 'o' num2str(o(1)) 'c' num2str(c(1)) 't' num2str(t(1)) 'p' num2str(p(1)) 'r' num2str(r(1))]);
%                     tagPos = ['C' tagPos3 'R' tagPos2 'S' tagPos1];

                                                                          
                    if useBipolar
                        
                        if isunix
                            folderSaveMP = [folderSaveMP 'bipolar/'];
                        else
                            folderSaveMP = [folderSaveMP 'bipolar\'];
                        end
                        makeDirectory(folderSaveMP);
                        
                        mpSavedSpectrumName = (['mpSpectrum' num2str(numAtomsMP) 'atoms' '_notch' num2str(notchData)...
                            'pos' tagPos...
                            '_bipolar' num2str(channelString) num2str(analogChannelString2)...
                            '_pivots' num2str(paramNum1) num2str(paramNum2)]);
                                                    
                        % Vinay - save or load saved MP data only if the
                        % saveMPFlag is set, else just do all the
                        % computations
                        if (exist([folderSaveMP mpSavedSpectrumName '.mat'],'file') && saveMPFlag)
                            disp('loading saved spectrum data...');
                            load([folderSaveMP mpSavedSpectrumName '.mat']);
                            L = size(analogData,2);
                        else     
                            clear tag gaborInfo X
                            tag = 'temp/';

                            disp('Running MP for the bipolar case...');

                            % Import the data
                            X = analogData';
                            L = size(analogData,2);
                            signalRange = [1 L]; % full range
                            importData(X,mpFolderTemp,tag,signalRange,Fs);

                            % perform Gabor decomposition
                            Numb_points = L; % length of the signal
                            Max_iterations = 100; % number of iterations

                            disp(['MP decomposition with max iterations =' num2str(Max_iterations)]);

                            runGabor(mpFolderTemp,tag,Numb_points, Max_iterations);
                            disp('MP decomposition done');

                            gaborInfo = getGaborData(mpFolderTemp,tag,1);

                            gaborInfoGoodPos = gaborInfo(goodPos);

                            wrap = [];
                            atomList = (1:numAtomsMP);
                            mpSpectrum = [];
                            disp(['Reconstructing Energy from:' num2str(numAtomsMP) 'atoms, and'  num2str(length(goodPos)) 'trials']);
                            for m=1:length(goodPos)
                                disp(['trial number:' num2str(m) '(actual trial -)' num2str(goodPos(m))]);
                                rEnergy = reconstructEnergyFromAtomsMPP(gaborInfoGoodPos{m}.gaborData,L,wrap,atomList);
                                if m == 1
                                    mpSpectrum = rEnergy;
                                else
                                    mpSpectrum = mpSpectrum + rEnergy;
                                end
                            end
                            mpSpectrum = mpSpectrum/length(goodPos);
                        
                            if saveMPFlag
                                disp('saving spectrum data...');
                                save([folderSaveMP mpSavedSpectrumName], 'mpSpectrum');
                            end
                        end
                        
                    else
                        
                        if isunix
                            folderSaveMP = [folderSaveMP num2str(channelString) '/'];
                        else
                            folderSaveMP = [folderSaveMP num2str(channelString) '\'];
                        end
                        makeDirectory(folderSaveMP);
                        
                        
                        mpSavedSpectrumName = (['mpSpectrum' num2str(numAtomsMP) 'atoms' '_notch' num2str(notchData)...
                            'pos' tagPos...
                            '_pivots' num2str(paramNum1) num2str(paramNum2)]);
                        
                        if (exist([folderSaveMP mpSavedSpectrumName '.mat'],'file') && saveMPFlag)
                            disp('loading saved spectrum data...');
                            load([folderSaveMP mpSavedSpectrumName '.mat']);
                            L = size(analogData,2);
                        else
                            clear tag gaborInfo
                            tag = channelString;
                            if ispc
                                load([mpFolder tag '\gaborInfo.mat']);
                            else
                                load([mpFolder tag '/gaborInfo.mat']);
                            end
                            gaborInfoGoodPos = gaborInfo(goodPos);

                            L = size(analogData,2);
                            wrap = [];
                            atomList = (1:numAtomsMP);
                            mpSpectrum = [];
                            disp(['Reconstructing Energy from:' num2str(numAtomsMP) 'atoms, and'  num2str(length(goodPos)) 'trials']);
                            for m=1:length(goodPos)
                                disp(['trial number:' num2str(m) '(actual trial -)' num2str(goodPos(m))]);
                                rEnergy = reconstructEnergyFromAtomsMPP(gaborInfoGoodPos{m}.gaborData,L,wrap,atomList);
                                if m == 1
                                    mpSpectrum = rEnergy;
                                else
                                    mpSpectrum = mpSpectrum + rEnergy;
                                end
                            end
                            mpSpectrum = mpSpectrum/length(goodPos);
                            
                            if saveMPFlag
                                disp('saving spectrum data...');
                                save([folderSaveMP mpSavedSpectrumName], 'mpSpectrum');
                            end
                        
                        end
                    end
                    
                    mpSpectrum = mpSpectrum'; % transpose so that first index corresponds to 'time' and second to 'freq'
                    
                    t = timeVals;
                    f = 0:Fs/L:Fs/2; 
                    
                    % plot the MP Spectrum
                    if plotStyle == 3 % line plot
                        disp('No line plot for MP method');
                        
                    elseif plotStyle == 1 % pcolor plot
                        
                        if (spectrumType == 1) % raw
                            pcolor(tfplotHandles(k,j),t,f,conv2Log(mpSpectrum')); 
                            shading interp;
                            caxis([cmin cmax]);
                            
                        elseif(spectrumType == 2) % difference from baseline                            
                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            mpSpectrumBL = mpSpectrum(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogmpSpectrumBL = mean(conv2Log(mpSpectrumBL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dmpSpectrum = 10*(conv2Log(mpSpectrum) - repmat(mlogmpSpectrumBL,size(mpSpectrum,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            pcolor(tfplotHandles(k,j),t,f,dmpSpectrum'); 
                            shading interp;
                            caxis([cmin cmax]);
                        end
                        
                    elseif plotStyle == 2 % imagesc plot
                        
                        if (spectrumType == 1) % raw
                            imagesc(t,f,conv2Log(mpSpectrum'),'Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                            
                        elseif(spectrumType == 2) % difference from baseline                            
                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            mpSpectrumBL = mpSpectrum(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogmpSpectrumBL = mean(conv2Log(mpSpectrumBL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dmpSpectrum = 10*(conv2Log(mpSpectrum) - repmat(mlogmpSpectrumBL,size(mpSpectrum,1),1)); 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            imagesc(t,f,dmpSpectrum','Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                        end
                        
                    end
                    
                end

                % Display title
                if (j==1 && k==1)
                    title(tfplotHandles(k,j),[titleParam1 'vs' titleParam2 ':' num2str(titleList1(k)) ',' num2str(titleList2(j))],'FontSize',titleFontSize);
                else
                    title(tfplotHandles(k,j),[num2str(titleList1(k)) ',' num2str(titleList2(j))],'FontSize',titleFontSize);
                end
                
            end  
        end
   end
        
end


