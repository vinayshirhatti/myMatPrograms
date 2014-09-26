% List of all protocols

function [expDates,protocolNames,stimTypes] = allProtocolsTestMicroelectrode

expDates{1} = '220513'; protocolNames{1} = 'GRF_001'; stimTypes{1} = 2;

expDates{2} = '270813'; protocolNames{2} = 'GRF_001'; stimTypes{2} = 2;

expDates{3} = '170214'; protocolNames{3} = 'GRF_001'; stimTypes{3} = 2;

expDates{4} = '290514'; protocolNames{4} = 'GRF_001'; stimTypes{4} = 2;

expDates{5} = '300514'; protocolNames{5} = 'GRF_001'; stimTypes{5} = 2; % Recorded from Rig1Display using Blackrock
expDates{6} = '300514'; protocolNames{6} = 'GRF_002'; stimTypes{6} = 2; % Recorded from Rig2Display using Labjack @ 5kHz
expDates{7} = '300514'; protocolNames{7} = 'GRF_003'; stimTypes{7} = 2; % Recorded from Rig2Display using Blackrock
expDates{8} = '300514'; protocolNames{8} = 'GRF_004'; stimTypes{8} = 2; % Recorded from Rig2Display using Labjack @ 20 kHz

expDates{9} = '090614'; protocolNames{9} = 'GRF_002'; stimTypes{9} = 2; % Recorded from Rig1Display using Blackrock, Green ITC18
expDates{10} = '090614'; protocolNames{10} = 'GRF_003'; stimTypes{10} = 2; % Recorded from Rig1Display using Blackrock, Green ITC18 w/o copper gnd
expDates{11} = '090614'; protocolNames{11} = 'GRF_004'; stimTypes{11} = 2; % Recorded from Rig1Display using Blackrock, Red ITC18

expDates{12} = '090614'; protocolNames{12} = 'GRF_005'; stimTypes{12} = 2; % Recorded from Rig1Display using Blackrock, Green ITC18, after fixing the wires. Also recording from photodiode 

expDates{13} = '100614'; protocolNames{13} = 'GRF_001'; stimTypes{13} = 2;
expDates{14} = '100614'; protocolNames{14} = 'GRF_002'; stimTypes{14} = 2; 
expDates{15} = '100614'; protocolNames{15} = 'GRF_003'; stimTypes{15} = 2; % Channel 1-64: simulated data, 65-66: EMG, Ainp3 - photodiode

expDates{16} = '110614'; protocolNames{16} = 'GRF_001'; stimTypes{16} = 2; % Bank C,Channel 1:Rt.occipital,approx O2,Ref:scalp centre,Ground:headpost clamp;Rig1Display,Blackrock+Photodiode;stim:{Flashing stimuli:full screen,varying Contrast & Tf[0.25:40}}
expDates{17} = '110614'; protocolNames{17} = 'GRF_002'; stimTypes{17} = 2; % As above; Stimulus{Rt side gabor}
expDates{18} = '110614'; protocolNames{18} = 'GRF_003'; stimTypes{18} = 2; % As above; Stimulus{lt side gabor}

expDates{19} = '130614'; protocolNames{19} = 'GRF_001'; stimTypes{19} = 2; %
expDates{20} = '130614'; protocolNames{20} = 'GRF_002'; stimTypes{20} = 2;
expDates{21} = '130614'; protocolNames{21} = 'GRF_003'; stimTypes{21} = 2; % gabors shown at 100% contrast; 4 locations (4 quadrants) ch 65 - right, ch 69 - left Occipital

expDates{22} = '170614'; protocolNames{22} = 'GRF_001'; stimTypes{22} = 2; % gabors shown at 100% contrast; 25 locations (5 azi x 5 ele) ch 65 - right, ch 66 - left Occipital

expDates{23} = '180614'; protocolNames{23} = 'GRF_001'; stimTypes{23} = 2; % monitor testing: BENQ; 6 contrasts and 6 temporal freq
expDates{24} = '180614'; protocolNames{24} = 'GRF_002'; stimTypes{24} = 2; % monitor testing: DELL; 6 contrasts and 6 temporal freq

expDates{25} = '190614'; protocolNames{25} = 'GRF_006'; stimTypes{25} = 2; % monitor testing: DELL - 6 contrasts, 6 temporal frequencies - 1.25 to 30 Hz - problem of all contrasts being 100%
expDates{26} = '190614'; protocolNames{26} = 'GRF_007'; stimTypes{26} = 2; % monitor testing: BENQ - 6 contrasts, 7 temporal frequencies - 1.25 to 50 Hz - problem of all contrasts being 100%
expDates{27} = '200614'; protocolNames{27} = 'GRF_001'; stimTypes{27} = 2; % Ashutosh EEG - 6 contrasts, 6 temporal frequencies - 1.25 to 30 Hz
expDates{28} = '200614'; protocolNames{28} = 'GRF_003'; stimTypes{28} = 2; % Ashutosh EEG - RF map; 5 azi 5 ele
expDates{29} = '200614'; protocolNames{29} = 'GRF_005'; stimTypes{29} = 2; % Ashutosh EEG - Eye Open/Close task: Stim on right - close, left - fixate

expDates{30} = '270614'; protocolNames{30} = 'GRF_001'; stimTypes{30} = 2; % monitor testing: DELL - 6 contrasts, 6 temporal frequencies - 1.25 to 30 Hz
expDates{31} = '270614'; protocolNames{31} = 'GRF_002'; stimTypes{31} = 2; % monitor testing: BENQ - 6 contrasts, 7 temporal frequencies - 1.25 to 50 Hz


expDates{32} = '070714'; protocolNames{32} = 'CRS_001'; stimTypes{32} = 2;

expDates{33} = '140714'; protocolNames{33} = 'CRS_001'; stimTypes{33} = 2; % Ring Protocol

expDates{34} = '230714'; protocolNames{34} = 'GRF_001'; stimTypes{34} = 2; % monitor testing: DELL - 6 contrasts, 6 temporal frequencies - 1.25 to 30 Hz, inside Faraday cage
expDates{35} = '230714'; protocolNames{35} = 'GRF_002'; stimTypes{35} = 2; % monitor testing: BENQ - 6 contrasts, 7 temporal frequencies - 1.25 to 50 Hz, inside Faraday cage

expDates{36} = '040814'; protocolNames{36} = 'GRF_001'; stimTypes{36} = 2; % Siddhesh; O1, O2, PO3, PO4; 4 ORI, 2 CONT, 2 TF; STIM TIME: 1000ms
expDates{37} = '040814'; protocolNames{37} = 'GRF_002'; stimTypes{37} = 2; % Siddhesh; O1, O2, PO3, PO4; 4 ORI, 2 CONT, 2 TF; STIM TIME: 1000ms
expDates{38} = '040814'; protocolNames{38} = 'GRF_003'; stimTypes{38} = 3; % Siddhesh; O1, O2, PO3, PO4; 4 ORI, 2 CONT, 2 TF; STIM TIME: 3000ms
expDates{39} = '040814'; protocolNames{39} = 'GRF_004'; stimTypes{39} = 3; % Siddhesh; O1, O2, PO3, PO4; 4 ORI, 2 CONT, 2 TF; STIM TIME: 3000ms; 
% stimTypes = 3 because the stim time is more here i.e. 3 secs, so need to
% extract a bigger data chunk, related adjustments made in
% runExtractAllData.m

expDates{40} = '060814'; protocolNames{40} = 'GRF_001'; stimTypes{40} = 2; % Alpa; protocol same as 36-37 (Siddhesh's)
expDates{41} = '060814'; protocolNames{41} = 'CRS_001'; stimTypes{41} = 2; % Alpa; EEG Electrodes same as above. Ring Protocol - 1 ori, 2 cont, 2 tf, 4 ring radii

expDates{42} = '070814'; protocolNames{42} = 'GRF_001'; stimTypes{42} = 2; % Alpa, Gabors flashed on either left or right half of the screen (3 cont, 2 ori, 4 tf)

expDates{43} = '130814'; protocolNames{43} = 'CRS_001'; stimTypes{43} = 4; % Alpa; EEG Electrodes. Annulus Fixed Protocol (pn = 10)
expDates{44} = '130814'; protocolNames{44} = 'CRS_002'; stimTypes{44} = 4; % Alpa; EEG Electrodes. Annulus Fixed Protocol (10)
expDates{45} = '130814'; protocolNames{45} = 'CRS_003'; stimTypes{45} = 4; % Alpa; EEG Electrodes. Annulus Fixed Protocol (10)
expDates{46} = '130814'; protocolNames{46} = 'CRS_004'; stimTypes{46} = 4; % Alpa; EEG Electrodes. Contrast Ring Protocol (2)

expDates{47} = '080914'; protocolNames{47} = 'CRS_001'; stimTypes{47} = 4; % Divija; EEG Electrodes: 10-20 system 19 electrodes. 
% Dual contrast protocol: 3x3 (CR) contrasts, 3 centre radii, 2 TFs. 2
% blocks (pn = 3)
expDates{48} = '080914'; protocolNames{48} = 'CRS_002'; stimTypes{48} = 4; % Divija; EEG Electrodes: 10-20 system 19 electrodes. 
% Dual contrast protocol: 3x3 (CR) contrasts, 3 centre radii, 2 TFs. 9
% blocks (pn = 3)
expDates{49} = '080914'; protocolNames{49} = 'CRS_003'; stimTypes{49} = 4; % Divija; EEG Electrodes: 10-20 system 19 electrodes. 
% Dual orientation protocol: 1x3 (CR) ori, 1 centre radius, 2 TFs, 3
% contrasts (pn = 4)


expDates{50} = '090914'; protocolNames{50} = 'CRS_001'; stimTypes{50} = 4; % Siddhesh; EEG Electrodes: 10-20 system 19 electrodes. 
% Dual orientation protocol: 1x3 (CR) ori, 1 centre radius, 2 TFs, 3
% contrasts. Same as 49 (pn = 4)
expDates{51} = '090914'; protocolNames{51} = 'CRS_002'; stimTypes{51} = 4; % Siddhesh; EEG Electrodes: 10-20 system 19 electrodes. 
% Dual contrast protocol: 3x3 (CR) contrasts, 3 centre radii, 2 TFs. Same
% as 48 (pn = 3)
expDates{52} = '090914'; protocolNames{52} = 'GRFAUD_001'; stimTypes{52} = 4; % Siddhesh; EEG Electrodes: 10-20 system 19 electrodes. 
% GaborRFAudio protocol, 8 sounds


expDates{53} = '100914'; protocolNames{53} = 'GRFAUD_001'; stimTypes{53} = 4; % Subhash; EEG Electrodes: 10-20 system 19 electrodes. 
% GaborRFAudio protocol, 8 sounds, similar to 52
expDates{54} = '100914'; protocolNames{54} = 'GRF_002'; stimTypes{54} = 4; % Subhash; EEG Electrodes: 10-20 system 19 electrodes
% GRF protocol with variations in contrast, sf, tf


expDates{55} = '180914'; protocolNames{55} = 'GRF_001'; stimTypes{55} = 4; % Murty; EEG Electrodes: 10-20 system 19 electrodes
% GRF protocol with variations in contrast, tf
expDates{56} = '180914'; protocolNames{56} = 'GRF_002'; stimTypes{56} = 4; % Murty; EEG Electrodes: 10-20 system 19 electrodes
% GRF protocol with variations in contrast, ori
expDates{57} = '180914'; protocolNames{57} = 'GRFAUD_003'; stimTypes{57} = 4; % Murty; EEG Electrodes: 10-20 system 19 electrodes
% GaborRFAudio protocol, 6 (3x2repeats) noise sounds
expDates{58} = '180914'; protocolNames{58} = 'CRS_004'; stimTypes{58} = 4; % Murty; EEG Electrodes: 10-20 system 19 electrodes
% Dual Orientation protocol, 3 contrasts, 4 orientations(ring)


expDates{59} = '230914'; protocolNames{59} = 'GRF_001'; stimTypes{59} = 4; % Sidrat; EEG Electrodes: 10-20 system 21 electrodes
% GRF protocol with variations in size
expDates{60} = '230914'; protocolNames{60} = 'CRS_002'; stimTypes{60} = 4; % Sidrat; EEG Electrodes: 10-20 system 21 electrodes
% CRS protocol with variations in size


