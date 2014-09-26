% Show EEG activity across all electrodes
%
% Vinay Shirhatti, 16 Sep 2014
%==========================================================================

function showEEGActivity(monkeyName,expDate,protocolName,folderSourceString,gridType,loadProtocolNumber)

if ~exist('folderSourceString','var')  folderSourceString='/media/store/';        end
if ~exist('gridType','var')            gridType='EEG';                  end
if ~exist('loadProtocolNumber','var')  loadProtocolNumber = 11;         end % Vinay - anything greater than 10 will read the protocolNumber from the saved data

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
elseif isunix
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

% Get Combinations
[~,aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,cValsUnique,tValsUnique,rValsUnique,pValsUnique] = loadParameterCombinations(folderExtract);

% Vinay - read the variable loadProtocolNumber as protocolNumber if it is
% less than 11. Otherwise it will be loaded from stimResults. This is so
% that one can directly pass the protocolNumber in case it is not recorded
% in stimResults
if (loadProtocolNumber < 11)
    protocolNumber = loadProtocolNumber;
else
    protocolNumber = getProtocolNumber(folderExtract);
end

% Vinay - get the particular gabors that were displayed in the protocol.
gaborsDisplayed = getGaborsDisplayed(folderExtract);
disp(['Gabors displayed in this protocol:' num2str(gaborsDisplayed)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display main options
% fonts
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16; fontSizeTiny = 8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Panels
panelHeight = 0.3; panelStartHeight = 0.68;
staticPanelWidth = 0.22; staticStartPos = 0.025; % [Vinay] - changed width 0.25 to 0.22
dynamicPanelWidth = 0.23; dynamicStartPos = 0.245; % [Vinay] - changed width 0.25 to 0.23
timingPanelWidth = 0.20; timingStartPos = 0.475; % [Vinay] - changed width 0.25 to 0.20
plotOptionsPanelWidth = 0.15; plotOptionsStartPos = 0.675; % [Vinay] - changed width 0.2 to 0.15
selectParamPanelWidth = 0.15; selectParamStartPos = 0.825; % [Vinay] - adding a panel to select parameters
backgroundColor = 'w';

% Vinay - define a flag for notching the line noise
notchData = 0;
% [Vinay] - define a flag to take bipolar (i.e. diff) signals: V1-V2
useBipolar = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dynamic panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynamicHeight = 0.06; dynamicGap=0.015; dynamicTextWidth = 0.6;
hDynamicPanel = uipanel('Title','Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[dynamicStartPos panelStartHeight dynamicPanelWidth panelHeight]);

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


%[Vinay] - create 3 columns corresponding to the parameters of the 3 gabors
%-Centre, Ring, Surround. Textwidth = 0.4, C,R,S popout - each 0.2 width

labelWidth = 0.4; crsWidth = 0.2; % [Vinay] - labelWidth + 3(crsWidth) = 1

uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[labelWidth 1-3*(dynamicHeight+dynamicGap) (3*crsWidth) dynamicHeight],...
    'Style','text','String','Surround        Ring        Centre','FontSize',fontSizeTiny);

%-----------------------
% numGabors = length(gaborsDisplayed); % Vinay - was using this to read 
% only the parameters for the non-null gabors. Not used anymore

for gaborNum = 1:3
    
    % Azimuth
    azimuthString = getStringFromValues(aValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-4*(dynamicHeight+dynamicGap) labelWidth dynamicHeight],...
    'Style','text','String','Azimuth (Deg)','FontSize',fontSizeTiny);
    hAzimuth(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-4*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',azimuthString,'FontSize',fontSizeTiny);

    % Elevation
    elevationString = getStringFromValues(eValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-5*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
    'Style','text','String','Elevation (Deg)','FontSize',fontSizeTiny);
    hElevation(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position',...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-5*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',elevationString,'FontSize',fontSizeTiny);
    
    % Sigma
    sigmaString = getStringFromValues(sValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-6*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
    'Style','text','String','Sigma (Deg)','FontSize',fontSizeTiny);
    hSigma(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-6*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',sigmaString,'FontSize',fontSizeTiny);

    % Spatial Frequency
    spatialFreqString = getStringFromValues(fValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-7*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
    'Style','text','String','Spatial Freq (CPD)','FontSize',fontSizeTiny);
    hSpatialFreq(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-7*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',spatialFreqString,'FontSize',fontSizeTiny);

    % Orientation
    orientationString = getStringFromValues(oValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-8*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
    'Style','text','String','Orientation (Deg)','FontSize',fontSizeTiny);
    hOrientation(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-8*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',orientationString,'FontSize',fontSizeTiny);

    % Contrast
    contrastString = getStringFromValues(cValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-9*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
    'Style','text','String','Contrast (%)','FontSize',fontSizeTiny);
    hContrast(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-9*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',contrastString,'FontSize',fontSizeTiny);

    % Temporal Frequency
    temporalFreqString = getStringFromValues(tValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-10*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
    'Style','text','String','Temporal Freq (Hz)','FontSize',fontSizeTiny);
    hTemporalFreq(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-10*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',temporalFreqString,'FontSize',fontSizeTiny);
    
    % [Vinay] - adding radius and spatial phase cells
    % radius
    radiusString = getStringFromValues(rValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-11*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
    'Style','text','String','Radius (deg)','FontSize',fontSizeTiny);
    hRadius(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-11*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',radiusString,'FontSize',fontSizeTiny);

    % spatial phase
    spatialPhaseString = getStringFromValues(pValsUnique,1, gaborNum);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-12*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
    'Style','text','String','Spatial Phase (Deg)','FontSize',fontSizeTiny);
    hSpatialPhase(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(labelWidth+(gaborNum-1)*crsWidth) 1-12*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
    'Style','popup','String',spatialPhaseString,'FontSize',fontSizeTiny);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Timing panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timingHeight = 0.1; timingTextWidth = 0.5; timingBoxWidth = 0.20; % [Vinay] - changed width from 0.25 to 0.20 
hTimingPanel = uipanel('Title','Timing','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[timingStartPos panelStartHeight timingPanelWidth panelHeight]);

% Vinay - changed these values for EEG
signalRange = [-0.8 1.1];
fftRange = [0 250];
baseline = [-0.4 0];
stimPeriod = [0.4 0.8];
cmin = -1; cmax = 3;

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

% color axis settings
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-8*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','caxis','FontSize',fontSizeSmall);
hTFcmin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [timingTextWidth 1-8*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',cmin,'FontSize',fontSizeTiny); % initialize cmin = -1
hTFcmax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [timingTextWidth+timingBoxWidth 1-8*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',cmax,'FontSize',fontSizeTiny); % initialize cmax = 3


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
    'Style','text','String','Color','FontSize',fontSizeSmall);

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
analysisTypeString = 'VoltageAmplitude|SpectralPower|InstantaneousPhase|InstantanousFreq|PhaseCoherence|AmplitudeCoherence';
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
    'Style','pushbutton','String','plot/play','FontSize',fontSizeMedium, ...
    'Callback',{@plotData_Callback});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Vinay] - adding a panel to select the parameters for the plots

selectParamHeight = 0.1; selectParamWidth = 0.5; selectParamBoxWidth = 0.20; % [Vinay] - changed width from 0.25 to 0.20
textWidth = 0.4; selGWidth = 0.2; selPWidth = 0.4;
hSelectParamPanel = uipanel('Title','Select Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[selectParamStartPos panelStartHeight selectParamPanelWidth panelHeight]);


parametersString = 'NA|azimuth|elevation|sigma|spatialFreq|orientation|contrast|temporalFreq|radius|spatialPhase';
gaborString = 'NA|C|R|S';

selPosString = 'singleTrial|Avg';

%Trial positions - i.e. good trials for the selected combination of
%parameters
uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'Position',[0 1-selectParamHeight textWidth selectParamHeight], ...
    'Style','pushbutton','String','Get trials','FontSize',fontSizeSmall,...
    'Callback',{@getTrialPos_Callback});

uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'Position',[textWidth 1-selectParamHeight 2*selGWidth selectParamHeight], ...
    'Style','text','String','Sel trial(s)','FontSize',fontSizeTiny);

goodPosString = '-';

hGoodPos = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(textWidth+2*selGWidth) 1-selectParamHeight 0.5*selPWidth selectParamHeight], ...
    'Style','popup','String',goodPosString,'FontSize',fontSizeTiny);


% Play timings
uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'Position',[0 1-3*selectParamHeight 1.5*textWidth selectParamHeight], ...
    'Style','text','String','Play Period','FontSize',fontSizeSmall);

playPeriodString = 'full|base|stim';
hPlayPeriod = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(textWidth+selGWidth) 1-3*selectParamHeight selPWidth selectParamHeight], ...
    'Style','popup','String',playPeriodString,'FontSize',fontSizeTiny);


%Select band
uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'Position',[0 1-4*selectParamHeight 1.5*textWidth selectParamHeight], ...
    'Style','text','String','Select Band','FontSize',fontSizeSmall);

bandString = 'full|custom';
hBand = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(textWidth+selGWidth) 1-4*selectParamHeight selPWidth selectParamHeight], ...
    'Style','popup','String',bandString,'FontSize',fontSizeTiny);

%Band limits
uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'Position',[0 1-5*selectParamHeight 1.5*textWidth selectParamHeight], ...
    'Style','text','String','Band Limits','FontSize',fontSizeSmall);

bandLow = 30; bandHigh = 80; % default limits
hBandLow = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [1.5*textWidth 1-5*selectParamHeight 0.5*selPWidth selectParamHeight], ...
    'Style','edit','String',bandLow,'FontSize',fontSizeTiny);
hBandHigh = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [1.5*textWidth+0.5*selPWidth 1-5*selectParamHeight 0.5*selPWidth selectParamHeight], ...
    'Style','edit','String',bandHigh,'FontSize',fontSizeTiny);


% display the specific protocol being used
protocolNameString = {'NoneProtocol';'RingProtocol';'ContrastRingProtocol';'DualContrastProtocol'; ...
    'DualOrientationProtocol';'DualPhaseProtocol';'OrientationRingProtocol';'PhaseRingProtocol'; ...
    'Drifting Ring Protocol';'CrossOrientationProtocol';'AnnulusRingProtocol'};
% protocolNumber = getProtocolNumber(folderExtract);

uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'Position',[0.1 1-10*selectParamHeight 0.8 selectParamHeight], ...
    'Style','text','String',protocolNameString(protocolNumber+1),'FontSize',fontSizeTiny,'FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get plots and message handles

% Get electrode array information
electrodeGridPos = [staticStartPos panelStartHeight staticPanelWidth panelHeight];
hElectrodes = showElectrodeLocations(electrodeGridPos,analogChannelsStored(get(hAnalogChannel,'val')), ...
    colorNames(get(hChooseColor,'val')),[],1,0,gridType);


% Set up the plot box and its dimensions
startXPos = staticStartPos; endXPos = 0.98; startYPos = 0.05; endYPos = 0.60;
centerGap = 0.04;
plotsWidthX = (endXPos-startXPos-centerGap/2);
plotsWidthY = (endYPos-startYPos-centerGap/4);
gap = 0.02; gapSmall = 0.002;

% play grid - to play and show the signal at all the electrodes together
playGridPos=[(startXPos+0.3*plotsWidthX) startYPos 0.4*plotsWidthX endYPos];

hPlayGrid = getPlotHandles(1,1,playGridPos);

% analog signal grid 1 and 2

analogSignalGridPos1 = [startXPos (startYPos+0.5*plotsWidthY+2*gap) 0.3*plotsWidthX-gap 0.5*plotsWidthY];

hAnalogGrid1 = getPlotHandles(1,1,analogSignalGridPos1);

analogSignalGridPos2 = [startXPos+0.7*plotsWidthX+gap (startYPos+0.5*plotsWidthY+2*gap) 0.3*plotsWidthX-gap 0.5*plotsWidthY];

hAnalogGrid2 = getPlotHandles(1,1,analogSignalGridPos2);


% spectral grid 1 and 2

spectralGridPos1 = [startXPos (startYPos) 0.3*plotsWidthX-gap 0.5*plotsWidthY];

hSpectralGrid1 = getPlotHandles(1,1,spectralGridPos1);

spectralGridPos2 = [startXPos+0.7*plotsWidthX+gap (startYPos) 0.3*plotsWidthX-gap 0.5*plotsWidthY];

hSpectralGrid2 = getPlotHandles(1,1,spectralGridPos2);

% Title
uicontrol('Unit','Normalized','Position',[0 0.975 1 0.025],...
    'Style','text','String',[monkeyName expDate protocolName],'FontSize',fontSizeLarge);



%__________________________________________________________________________
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%----------------------CALLBACK FUNCTIONS----------------------------------
%==========================================================================

function plotData_Callback(~,~)
        
        a=cell2mat(get(hAzimuth,'val'));
        e=cell2mat(get(hElevation,'val'));
        s=cell2mat(get(hSigma,'val'));
        f=cell2mat(get(hSpatialFreq,'val'));
        o=cell2mat(get(hOrientation,'val'));
        c=cell2mat(get(hContrast,'val'));
        t=cell2mat(get(hTemporalFreq,'val'));
        r=cell2mat(get(hRadius,'val'));
        p=cell2mat(get(hSpatialPhase,'val'));
        
        analysisType = get(hAnalysisType,'val');
        plotColor = colorNames(get(hChooseColor,'val'));
        BLMin = str2double(get(hBaselineMin,'String'));
        BLMax = str2double(get(hBaselineMax,'String'));
        STMin = str2double(get(hStimPeriodMin,'String'));
        STMax = str2double(get(hStimPeriodMax,'String'));
        
        holdOnState = get(hHoldOn,'val');
        
        notchData = get(hNotchData, 'val');
        
        plotLineWidth = linewidths(get(hChooseWidth,'val')); % Vinay - for plot line width
        
        
        % [Vinay] - get the plot handles here again
        hPlayGrid = getPlotHandles(1,1,playGridPos);

        % analog signal grid 1 and 2

        hAnalogGrid1 = getPlotHandles(1,1,analogSignalGridPos1);

        hAnalogGrid2 = getPlotHandles(1,1,analogSignalGridPos2);

        % spectral grid 1 and 2

        hSpectralGrid1 = getPlotHandles(1,1,spectralGridPos1);

        hSpectralGrid2 = getPlotHandles(1,1,spectralGridPos2);
        
        
        labelMatrix = gridMatrix(gridType);
        
        playPos = get(hGoodPos,'val');
        
        cmin = str2double(get(hTFcmin,'String'));
        cmax = str2double(get(hTFcmax,'String'));
        
        playPeriod = get(hPlayPeriod,'val');


        if analysisType==1 % plot analog signal recorded
            
            analogChannelPos = get(hAnalogChannel,'val');
            analogChannelString = analogChannelStringArray{analogChannelPos};
            analogChannelsStored1 = analogChannelsStored(analogChannelPos);
            
            % [Vinay] - read the analog channel 2 string for bipolar case
            analogChannelPos2 = get(hAnalogChannel2,'val');
            analogChannelString2 = ('none');
            analogChannelsStored2 = [];
            if analogChannelPos2 ~= 1
                analogChannelString2 = analogChannelStringArray{analogChannelPos2-1};
                analogChannelsStored2 = analogChannelsStored(analogChannelPos2-1);
            end
            
            playLFPData(hPlayGrid,hAnalogGrid1,hAnalogGrid2,hSpectralGrid1,hSpectralGrid2,analogChannelStringArray,analogChannelsStored1,analogChannelsStored2,labelMatrix,playPos,a,e,s,f,o,c,t,r,p,folderLFP,...
                analysisType,timeVals,plotColor,plotLineWidth,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar,cmin,cmax,playPeriod);
            
            

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

        rescaleData(plotHandles,xMin,xMax,getYLims(plotHandles));
        rescaleData(hParam1Plot,xMin,xMax,getYLims(hParam1Plot));
        rescaleData(hParam2Plot,xMin,xMax,getYLims(hParam2Plot));
        rescaleData(hParam3Plot,xMin,xMax,getYLims(hParam3Plot));
        rescaleData(hParam4Plot,xMin,xMax,getYLims(hParam4Plot));
        rescaleData(hParam5Plot,xMin,xMax,getYLims(hParam5Plot));
        rescaleData(hParam6Plot,xMin,xMax,getYLims(hParam6Plot));
        showElectrodeLocations(electrodeGridPos,channelNumber,plotColor,hElectrodes,holdOnState,0,gridType);
end

%--------------------------------------------------------------------------

function getTrialPos_Callback(~,~)
    
if ispc
    folderExtract = [folderName 'extractedData\'];
    folderSegment = [folderName 'segmentedData\'];
else
    folderExtract = [folderName 'extractedData/'];
    folderSegment = [folderName 'segmentedData/'];
end

% Get parameters' combinations
parameterCombinations = loadParameterCombinations(folderExtract);

% Get bad trials
badTrialFile = [folderSegment 'badTrials.mat'];
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end

a=cell2mat(get(hAzimuth,'val'));
e=cell2mat(get(hElevation,'val'));
s=cell2mat(get(hSigma,'val'));
f=cell2mat(get(hSpatialFreq,'val'));
o=cell2mat(get(hOrientation,'val'));
c=cell2mat(get(hContrast,'val'));
t=cell2mat(get(hTemporalFreq,'val'));
r=cell2mat(get(hRadius,'val'));
p=cell2mat(get(hSpatialPhase,'val'));
    
clear goodPos
% [Vinay] - goodPos will be obtained by taking the goodPos for all
% the three gabors and finding the common trials in them
clear pos1 pos2 pos3
pos1 = parameterCombinations{a(1),e(1),s(1),f(1),o(1),c(1),t(1),p(1),r(1),1}; % good positions for gabor 1 i.e. S
pos2 = parameterCombinations{a(2),e(2),s(2),f(2),o(2),c(2),t(2),p(2),r(2),2}; % for R
pos3 = parameterCombinations{a(3),e(3),s(3),f(3),o(3),c(3),t(3),p(3),r(3),3}; % for C
goodPos = intersect(pos1,pos2);
goodPos = intersect(goodPos,pos3);
goodPos = setdiff(goodPos,badTrials);

goodPosString = [];
for i = 1:length(goodPos)
    if i == 1
        goodPosString = num2str(goodPos(i));
    else
        goodPosString = [goodPosString '|' num2str(goodPos(i))];
    end
end
goodPosString = ['all|' goodPosString];


hGoodPos = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(textWidth+2*selGWidth) 1-selectParamHeight 0.5*selPWidth selectParamHeight], ...
    'Style','popup','String',goodPosString,'FontSize',fontSizeTiny);

writeGoodPos = ['No. of good trials:' num2str(length(goodPos))];

uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
    'Position',[0 1-2*selectParamHeight 2*textWidth selectParamHeight], ...
    'Style','text','String',writeGoodPos,'FontSize',fontSizeSmall);

end

%--------------------------------------------------------------------------

function bipolar_Callback(hObject,~,~)
   if (get(hObject,'val') ~= 1)
       useBipolar = 1;
   else
       useBipolar = 0;
   end
end

%--------------------------------------------------------------------------

function notch_Callback(hObject,~,~)
       notchData = get(hObject,'val');
end

%--------------------------------------------------------------------------




end


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Supporting functions
%--------------------------------------------------------------------------
function [analogChannelsStored,timeVals,goodStimPos,analogInputNums] = loadlfpInfo(folderLFP) %#ok<*STOUT>
load([folderLFP 'lfpInfo']);
if ~exist('analogInputNums','var')
    analogInputNums=[];
end
end
%--------------------------------------------------------------------------

function [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique, rValsUnique, pValsUnique] = loadParameterCombinations(folderExtract)
% [Vinay] - added rValsUnique and pValsUnique for radius and spatial phase
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

%--------------------------------------------------------------------------

function protocolNumber = getProtocolNumber(folderExtract)
load ([folderExtract 'stimResults']);
protocolNumber = stimResults.protocolNumber;
end

%--------------------------------------------------------------------------

function gaborsDisplayed = getGaborsDisplayed(folderExtract)
load ([folderExtract 'stimResults']);
gaborsDisplayed = stimResults.side;
end

%--------------------------------------------------------------------------

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

%--------------------------------------------------------------------------

function outString = getStringFromValues(valsUnique,decimationFactor, gaborNum)

if length(valsUnique{gaborNum})==1
    outString = convertNumToStr(valsUnique{gaborNum}(1),decimationFactor);
else
    outString='';
    for i=1:length(valsUnique{gaborNum})
        outString = cat(2,outString,[convertNumToStr(valsUnique{gaborNum}(i),decimationFactor) '|']);
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

%--------------------------------------------------------------------------

function [colorString, colorNames] = getColorString

colorNames = 'brkgcmy';
colorString = 'blue|red|black|green|cyan|magenta|yellow';

end

%--------------------------------------------------------------------------

function goodPos = getTrialPositions(folderName,a,e,s,f,o,c,t,r,p)
    
if ispc
    folderExtract = [folderName 'extractedData\'];
    folderSegment = [folderName 'segmentedData\'];
else
    folderExtract = [folderName 'extractedData/'];
    folderSegment = [folderName 'segmentedData/'];
end

% Get parameters' combinations
parameterCombinations = loadParameterCombinations(folderExtract);

% Get bad trials
badTrialFile = [folderSegment 'badTrials.mat'];
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end
    
clear goodPos
% [Vinay] - goodPos will be obtained by taking the goodPos for all
% the three gabors and finding the common trials in them
clear pos1 pos2 pos3
pos1 = parameterCombinations{a(1),e(1),s(1),f(1),o(1),c(1),t(1),p(1),r(1),1}; % good positions for gabor 1 i.e. S
pos2 = parameterCombinations{a(2),e(2),s(2),f(2),o(2),c(2),t(2),p(2),r(2),2}; % for R
pos3 = parameterCombinations{a(3),e(3),s(3),f(3),o(3),c(3),t(3),p(3),r(3),3}; % for C
goodPos = intersect(pos1,pos2);
goodPos = intersect(goodPos,pos3);
goodPos = setdiff(goodPos,badTrials);

end

%--------------------------------------------------------------------------

function badTrials = loadBadTrials(badTrialFile)
load(badTrialFile);
end

%--------------------------------------------------------------------------

function labelMatrix = gridMatrix(gridType)

if strcmp(gridType,'Microelectrode')

    labelMatrix = ...
       [00 02 01 03 04 06 08 10 14 00;
        65 66 33 34 07 09 11 12 16 18;
        67 68 35 36 05 17 13 23 20 22;
        69 70 37 38 48 15 19 25 27 24;
        71 72 39 40 42 50 54 21 29 26;
        73 74 41 43 44 46 52 62 31 28;
        75 76 45 47 51 56 58 60 64 30;
        77 78 82 49 53 55 57 59 61 32;
        79 80 84 86 87 89 91 94 63 95;
        00 81 83 85 88 90 92 93 96 00];
      
elseif strcmp(gridType,'ECoG')

    labelMatrix = ...
       [00 89 81 73 65 25 17 09 01 00;
        00 90 82 74 66 26 18 10 02 00;
        00 91 83 75 67 27 19 11 03 00;
        34 92 84 76 68 28 20 12 04 33;
        00 93 85 77 69 29 21 13 05 00;
        00 94 86 78 70 30 22 14 06 00;
        00 95 87 79 71 31 23 15 07 00;
        00 96 88 80 72 32 24 16 08 00];

elseif strcmp(gridType,'EEG') % [Vinay] - adding the EEG 10-20 system grid
    
    labelMatrix = ...
       [00 01 00 02 00;
        03 15 17 16 04;
        05 13 18 14 06;
        07 11 19 12 08;
        00 09 00 10 00];
end    


end

%--------------------------------------------------------------------------

function goodPosAnalogData = loadAnalogData(folderLFP,channelString,notchData,goodPos,time)

clear analogData analogDataNotched
load([folderLFP 'elec' channelString '.mat']);

if notchData
    goodPosAnalogData = mean(analogDataNotched(goodPos,:),1);
else
    goodPosAnalogData = mean(analogData(goodPos,:),1);
end

if exist('time','var')
    goodPosAnalogData = goodPosAnalogData(:,time);
end

end

%--------------------------------------------------------------------------

function greyColor = getGrey(plotHandles,saturation)

        greyColor = get(plotHandles,'Color');
        greyColor(greyColor==0) = saturation;
end

%--------------------------------------------------------------------------

function playLFPData(hPlayGrid,hAnalogGrid1,hAnalogGrid2,hSpectralGrid1,hSpectralGrid2,analogChannelStringArray,analogChannelsStored1,analogChannelsStored2,labelMatrix,playPos,a,e,s,f,o,c,t,r,p,folderLFP,...
    analysisType,timeVals,plotColor,plotLineWidth,BLMin,BLMax,STMin,STMax,folderName,protocolNumber, notchData,useBipolar,cmin,cmax,playPeriod)

if ispc
    folderExtract = [folderName 'extractedData\'];
    folderSegment = [folderName 'segmentedData\'];
else
    folderExtract = [folderName 'extractedData/'];
    folderSegment = [folderName 'segmentedData/'];
end

timePoints = length(timeVals);

% Get the play duration etc
Fs = round(1/(timeVals(2)-timeVals(1)));
BLRange = uint16((BLMax-BLMin)*Fs);
STRange = uint16((STMax-STMin)*Fs);
BLPos = find(timeVals>=BLMin,1)+ (1:BLRange);
STPos = find(timeVals>=STMin,1)+ (1:STRange);

xsBL = 0:1/(BLMax-BLMin):Fs-1/(BLMax-BLMin);
xsST = 0:1/(STMax-STMin):Fs-1/(STMax-STMin);

switch playPeriod
    case 1
        startTime = 1; endTime = timePoints;
    case 2
        startTime = BLPos(1); endTime = BLPos(end);
    case 3
        startTime = STPos(1); endTime = STPos(end);
end

% Plot the main play screen
numRows = size(labelMatrix,1);
numCols = size(labelMatrix,2);
timePoints = endTime - startTime + 1;

% initialize the matrix that will hold the signal values across electrodes
playMatrix = zeros(numRows,numCols,timePoints);


% Get the good trials
goodPos = getTrialPositions(folderName,a,e,s,f,o,c,t,r,p);

% Get the particular selected good trial if it is selected
if playPos ~= 1
    goodPos = goodPos(playPos-1);
end

% construct the signal matrix
    for row = 1:numRows
        for col = 1:numCols
            channelString = num2str(labelMatrix(row,col));
            if str2num(channelString) == 0
                playData = zeros(1,timePoints);
            else
                goodPosAnalogData = loadAnalogData(folderLFP,channelString,notchData,goodPos,startTime:endTime);
                playData = mean(goodPosAnalogData,1); % if it is a single trial then the mean doesn't change the signal, otherwise it averages across trials
            end
            playMatrix(row,col,:) =  playData;
        end
    end
%     set(hPlayGrid,'Nextplot','replace');

% Plot the signal
if useBipolar
    playSignal1 = loadAnalogData(folderLFP,num2str(analogChannelsStored1),notchData,goodPos,startTime:endTime);
    
    playSignal2 = loadAnalogData(folderLFP,num2str(analogChannelsStored2),notchData,goodPos,startTime:endTime);
    
    playSignalTrial = playSignal1 - playSignal2;
    
else
    playSignalTrial = loadAnalogData(folderLFP,num2str(analogChannelsStored1),notchData,goodPos,startTime:endTime);
end


% Plot the spectrum for the chosen trial(s)




% Get the good trials again to calculate the average plots
goodPos = getTrialPositions(folderName,a,e,s,f,o,c,t,r,p);
    
% Plot the signal
if useBipolar
    playSignal1 = loadAnalogData(folderLFP,num2str(analogChannelsStored1),notchData,goodPos,startTime:endTime);
    
    playSignal2 = loadAnalogData(folderLFP,num2str(analogChannelsStored2),notchData,goodPos,startTime:endTime);
    
    playSignal = playSignal1 - playSignal2;
    
else
    playSignal = loadAnalogData(folderLFP,num2str(analogChannelsStored1),notchData,goodPos,startTime:endTime);
end


% Play all the plots now
for t = 1:endTime-startTime+1
    imagesc(playMatrix(:,:,t),'Parent',hPlayGrid); colormap('copper');
    axis(hPlayGrid,'xy'); % for imagesc every 
    % parameter has to be set using the plot handles 
    % explicitly 
    set(hPlayGrid,'CLim',[cmin cmax]);
%     set(hPlayGrid,'Nextplot','add');
    text(0.2, 0.9, ['timeVal:' num2str(timeVals(t+startTime-1))],'FontSize',14,'Parent',hPlayGrid);
%     pause(0.01);

    set(hAnalogGrid1,'Nextplot','replace');
    plot(hAnalogGrid1,timeVals(startTime:endTime),playSignalTrial,'color',plotColor,'Linewidth',plotLineWidth);
    % Generate a time marker line to be plotted on the analog signal plots
    markerLine1 = get(hAnalogGrid1,'ylim');
    line([timeVals(t+startTime-1) timeVals(t+startTime-1)],markerLine1,'Color','k','Parent',hAnalogGrid1);
    axis('tight');
    
    set(hAnalogGrid1,'Nextplot','add');
    plot(hAnalogGrid1,timeVals(startTime:endTime),playSignal,'color',plotColor,'Linewidth',plotLineWidth);
    markerLine2 = get(hAnalogGrid2,'ylim');
    line([timeVals(t+startTime-1) timeVals(t+startTime-1)],markerLine2,'Color','k','Parent',hAnalogGrid2);
    axis('tight');
    
%     hAllFigs = findall(0,'Type','figure');
    
    pause(0.001);

%     mGrid(t) = getframe(hAllFigs);
    disp(['frame:' num2str(t)]);
end

% disp('saving movie...');
% movie2avi(mGrid,'eegplay.avi','compression','none');
% disp('movie saved...');

end

