function [Energy,xs,ys,sumEnergy] = reconstructEnergyFromAtoms(GaborData,N,Fs,wrap,Natoms,DecimationFactor)

% GaborData(1,:) - atom octave
% GaborData(2,:) - atom frequency in Hz
% GaborData(3,:) - atom time in sec
% GaborData(4,:) - atom modulus
% GaborData(5,:) - atom phase

if DecimationFactor==1
    sumEnergy = zeros(N,N);
else
    sumEnergy = zeros(2*N/DecimationFactor,N/DecimationFactor);
end

for i=1:Natoms
    
    oct = GaborData(1,i);
    ksi = GaborData(2,i);
    u   = GaborData(3,i);
    mod = GaborData(4,i);
    fi  = GaborData(5,i);
    
    if DecimationFactor==1
        [E,xs,ys] = reconstructWg(oct,ksi,u,N,Fs,wrap);
    else
        [E,xs,ys] = reconstructWg(oct,ksi,u,N,Fs,wrap,DecimationFactor,DecimationFactor/2);
    end
    Energy{i} = E*mod^2;
    sumEnergy = sumEnergy+Energy{i};
end