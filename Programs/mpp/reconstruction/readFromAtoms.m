% This program reads the gabor atom file and returns the desired values

function [Book, Natoms, Nwin] = readFromAtoms(foldername,tag,Tno,Cno)

[X,Y] = compareinout(foldername,tag);
[Lsig,Ntrials,Nchans] = size(X);
%disp(['Number of trials: ' num2str(Ntrials) ', Number of channels: ' num2str(Nchans)]);


fn = [foldername tag 'Atoms.txt'];
Data = textread(fn);

[M,N] = size(Data);
if M==6
    Data = Data';
end


ActualRawSignal = X(:,Tno,Cno);
ActualFilteredSignal = Y(:,Tno,Cno);
        
% Find the correct location
Startpoints = find(Data(:,4) == Lsig);
DesiredPoint = Nchans*(Tno-1)+Cno;
Stp = Startpoints(DesiredPoint);
A = Data(Stp,:);
Etp = Stp+A(3);

Book = Data(Stp+1:Etp,:)';
Natoms = A(3);
Nwin = A(4);