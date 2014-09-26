% This program is derived from runGabor. This program, however, does not
% run the gabord code. 

% foldername:  The folder where MP directories are created. 
% tag: An additional string to distinguish the data folder.  
% Numb_points is the length of the signal in points (default = 1024);
% Max_iterations = maximum number of Gabor atoms for each signal. (default 500)

% Note:
% Please change the sourcefolder variable on line 21 to the folder where
% the executable gabord is located.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray, 2008 
% Distributed under the General Public License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function prepareMPForOrchestra(foldername,tag,Numb_points, Max_iterations,folderNameOrchestra)


if ~exist('Numb_points','var')            Numb_points=1024;       end
if ~exist('Max_iterations','var')         Max_iterations=500;     end
if ~exist('folderNameOrchestra','var')
    folderNameOrchestra= foldername;       
end

% foldername=appendIfNotPresent(foldername,'\');
foldername=appendIfNotPresent(foldername,'/'); % Vinay - changed for linux
tagLinux=appendIfNotPresent(removeIfPresent(tag,'\'),'/');
% tag=appendIfNotPresent(tag,'\');
tag=appendIfNotPresent(tag,'/'); % Vinay - for linux
folderNameOrchestra = appendIfNotPresent(folderNameOrchestra,'/');

% fnin = [foldername tag];
fnin = [foldername tagLinux]; % Vinay - for linux
fninOrchestra = [folderNameOrchestra tagLinux];

% fn = [fnin 'ImportData_SIG\GaborMP\'];
fn = [fnin 'ImportData_SIG/GaborMP/']; % Vinay - for linux
makeDirectory(fn);
% fnout = [fnin 'GaborMP\'];
fnout = [fnin 'GaborMP/']; % Vinay - for linux
makeDirectory(fnout);
fnoutOrchestra = [fninOrchestra 'GaborMP/'];


Shift_points = Numb_points;

% Get information about the number of channels from the header file
% filename = [fnin 'ImportData_SIG\sig.hdr'];
filename = [fnin 'ImportData_SIG/sig.hdr'];
Numb_chans = getFromFile(filename,'Numb_chans');
All_chans = Numb_chans;


% Create the control file (called local.ctl)
filename = [fn 'local.ctl'];
fp = fopen(filename,'w');

fprintf(fp,'%s\n','# Template for running Gabor MP Analysis');
fprintf(fp,'%s\n','%INPUT_OUTPUT');
fprintf(fp,'%s\n','Numb_inputs=1');
fprintf(fp,'%s\n','Numb_outputs=1');
fprintf(fp,'%s\n','Mode=parallel');
fprintf(fp,'%s\n','%INPUT');
fprintf(fp,'%s%s\n','Path=',[fninOrchestra 'ImportData_SIG/']);
fprintf(fp,'%s\n','Header_file=sig.hdr');
fprintf(fp,'%s\n','Calibrate=1');
fprintf(fp,'%s%d\n','Numb_points=',Numb_points);
fprintf(fp,'%s%d\n','Shift_points=',Shift_points);
fprintf(fp,'%s\n','%OUTPUT');
fprintf(fp,'%s%s\n','Path=',fnoutOrchestra);
fprintf(fp,'%s%d\n','All_chans=',All_chans);
fprintf(fp,'%s%d\n','Numb_chans=',Numb_chans);
fprintf(fp,'%s\n','Start_chan=1');
fprintf(fp,'%s\n','Start_chan_no=0');
fprintf(fp,'%s\n','Header_file=book.hdr');
fprintf(fp,'%s\n','Type=book');
fprintf(fp,'%s\n','File_format=double');
fprintf(fp,'%s\n','Name_template=mp#.bok');
fprintf(fp,'%s\n','Max_len=600');
fprintf(fp,'%s\n','Chans_per_file=-1');
fprintf(fp,'%s\n','%GABOR_DECOMPOSITION');
%fprintf(fp,'%s\n','#Energy=50& int & &&Table=Gabor_MP&');
%fprintf(fp,'%s\n','#Energy=50& int & &&Table=Gabor_MP&');
fprintf(fp,'%s%d\n','Max_Iterations=',Max_iterations);
%fprintf(fp,'%s\n','#Coherence=1& int & &&Table=Gabor_MP&');
fclose(fp);

end