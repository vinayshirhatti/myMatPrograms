function atoms = reconstructFromAtoms(GaborData,N,Fs,wrap,Natoms)

% GaborData(1,:) - atom octave
% GaborData(2,:) - atom frequency in Hz
% GaborData(3,:) - atom time in sec
% GaborData(4,:) - atom modulus
% GaborData(5,:) - atom phase


for i=1:Natoms
    
    oct = GaborData(1,i);
    ksi = GaborData(2,i);
    u   = GaborData(3,i);
    mod = GaborData(4,i);
    fi  = GaborData(5,i);
    
    [signal,xs] = reconstructG(oct,ksi,u,fi,N,Fs,wrap);
    
    atoms(i,:) = signal*mod;
end