% runEnergy(foldername,tag,DecimationFactor,Exclude_Freq)
% This program writes the control file for running enerygy and then runs it.

function runEnergy(foldername,tag,DecimationFactor,Exclude_Freq)

fnin = [foldername tag];
fn = [fnin 'GaborMP/Energy1a/'];
makeDirectory(fn);
fnout = [fnin 'Energy1a/'];
makeDirectory(fnout);

if nargin == 2
	DecimationFactor = 1;
	L_Exclude_Freq = 0;
elseif nargin == 3
	L_Exclude_Freq = 0;
else
	[L_Exclude_Freq, junk] = size(Exclude_Freq);
end

% Get information about the number of channels from the header file (book.hdr)
filename = [fnin 'GaborMP/book.hdr'];
Numb_chans = getFromFile(filename,'Numb_chans');
All_chans = Numb_chans;
Max_Iterations = getFromFile(filename,'Max_Iterations');

% Create the control file (called local.ctl)
filename = [fn 'local.ctl'];
fp = fopen(filename,'w');

fprintf(fp,'%s\n','# Template for doing EMP energy analysis');
fprintf(fp,'%s\n','%INPUT_OUTPUT');
fprintf(fp,'%s\n','Numb_inputs=1');
fprintf(fp,'%s\n','Numb_outputs=1');
fprintf(fp,'%s\n','Mode=parallel');

fprintf(fp,'%s\n','%INPUT');
fprintf(fp,'%s%s\n','Path=',[fnin 'GaborMP/']);
fprintf(fp,'%s\n','Header_file=book.hdr');
fprintf(fp,'%s\n','Type=Book');

fprintf(fp,'%s\n','%OUTPUT');
fprintf(fp,'%s%s\n','Path=',fnout);
fprintf(fp,'%s%d\n','All_chans=',All_chans);
fprintf(fp,'%s%d\n','Numb_chans=',Numb_chans);
fprintf(fp,'%s\n','Start_chan=1');
fprintf(fp,'%s\n','Start_chan_no=0');
fprintf(fp,'%s\n','Type=Energy');
fprintf(fp,'%s\n','Header_file=emp.hdr');
fprintf(fp,'%s\n','Name_template=emp#.en');
fprintf(fp,'%s\n','Chans_per_file=-1');

fprintf(fp,'%s\n','%EMP_DATA');
fprintf(fp,'%s\n','Precision=1.e-10');
fprintf(fp,'%s%d\n','Max_Iterations=',Max_Iterations);
fprintf(fp,'%s\n','Wrap=1');
fprintf(fp,'%s\n','Average=0');
fprintf(fp,'%s\n','Square=1');
fprintf(fp,'%s%d\n','Decimate=',DecimationFactor);
fprintf(fp,'%s\n','Energy_frac=1.5');

if L_Exclude_Freq > 0
	fprintf(fp,'%s\n','%EMP_EXCLUDE');
	fprintf(fp,'%s%d\n','Frequency=',L_Exclude_Freq);
	for i=1:L_Exclude_Freq
		fprintf(fp,'%s%d%s%d%s\n','(',Exclude_Freq(i,1),',',Exclude_Freq(i,2),')');
	end
end

fclose(fp);

% Run the Energy program (Call energy)
sourcefolder = '/Users/sray/Research/timeFrequency/mpp/programs/source/';
commandline = [sourcefolder 'empdat ' filename];
unix(commandline);