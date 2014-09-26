% This program reconstructs the Energy (Gabor/Fourier/Dirac depending
% on the value of oct. It returns the Energy and normalization value.

% This program is similar to ReconstructG, which reconstructs the waveform
% instead of the energy.

% Inputs:
% N - size of window
% oct - octave (log2(scale)      0,..., log2(N)
% u - position                   %seconds    is converted to 0,..N-1
% ksi - frequency                %Hz         is converted to 0,..N-1
% fi - angle                    -pi,... pi
% wrap - whether the functions should be wrapped. Default = 1;


function [Energy,xs,ys] = reconstructWg(oct,ksi,u,N,Fs,wrap,DecimateTime,DecimateFreq)

if ~exist('Fs')             Fs=1000;            end
if ~exist('wrap')           wrap=1;             end
if ~exist('DecimateTime')   DecimateTime=1;     end
if ~exist('DecimateFreq')   DecimateFreq=1;     end

u   = u*Fs; 
ksi = round(ksi*N/Fs);          % Takes values between 0 and N-1
    
if oct==0           % Dirac
    Energy = zeros(N,N);
    Energy(:,u+1) = 1/N;

elseif oct==nextpow2(N)       % Fourier
    Energy = zeros(N,N);
    
    if ksi==0 | ksi==N/2
        Energy(ksi+1,:) = 1/N;
    else
        Energy(ksi+1,:) = 1/(2*N);
        Energy(N-ksi+1,:) = 1/(2*N);
    end

else                            % Gabor
    
    if ~wrap
        Energy = WignerVille(oct,u,ksi,N,Fs);
    else % Wrap in time
        
        Energy0 = WignerVille(oct,u,ksi,N,Fs);
        Energy1 = WignerVille(oct,u-N,ksi,N,Fs);
        Energy_1 = WignerVille(oct,u+N,ksi,N,Fs);
        
        Energy  = Energy0+Energy1+Energy_1;
    end
end

% Decimate
xs = ((DecimateTime-1)/2:DecimateTime:N)/Fs; xs=xs(1:N/DecimateTime);
Energy_T = squeeze(sum(reshape(Energy,[N DecimateTime N/DecimateTime]),2))';
ys = ((DecimateFreq-1)/2: DecimateFreq: N)*Fs/N; ys=ys(1:N/DecimateFreq);
clear Energy
Energy = squeeze(sum(reshape(Energy_T,[N/DecimateTime DecimateFreq N/DecimateFreq]),2))';

% Normalize
Energy = Energy/sum(sum(Energy));

end


function [X,xs,ys] = WignerVille(oct,u,ksi,N,Fs)

if ~exist('Fs')     Fs = 1000;      end

s = 2^oct;

if ksi > N/2
    disp('Signal is undersampled. ksi should not exceed N/2');
end


% Main window
t = 0:N-1;            
w = 0:N-1;  

X1 = WG((t-u)/s, (2*pi/N)*(w-ksi)*s);   % Main 
X2 = WG((t-u)/s, (2*pi/N)*(w-N-ksi)*s); % Aliasing
Xp = X1+X2;                     

% Computation of -ksi
X3 = WG((t-u)/s, (2*pi/N)*(w+ksi)*s);   % Aliasing
X4 = WG((t-u)/s, (2*pi/N)*(w-N+ksi)*s); % Main
Xn = X3+X4;                     

X = 0.5*(Xp+Xn);
xs = t/Fs;
ys = w*Fs/N;

end

function S = WG(t,w)

    TimeAxis = exp(-2*pi*t.*t);
    FreqAxis = exp(-w.*w/(2*pi))'; % equation 60 in mallat and Zhang
    S = FreqAxis*TimeAxis;
end