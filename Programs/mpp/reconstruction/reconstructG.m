% This program reconstructs the basis signal (Gabor/Fourier/Dirac depending
% on the value of oct. It returns the signal after normalization.

% Inputs:
% N - size of window
% oct - octave (log2(scale)      0,..., log2(N)
% u - position                   %seconds    is converted to 0,..N-1
% ksi - frequency                %Hz         is converted to 0,..N-1
% fi - angle                    -pi,... pi
% wrap - whether the functions should be wrapped. Default = 0;


function [signal,xs] = reconstructG(oct,ksi,u,fi,N,Fs,wrap,decimateTime)

if ~exist('Fs','var')             Fs=1000;            end
if ~exist('wrap','var')           wrap=1;             end
if ~exist('decimateTime','var')   decimateTime=1;     end

u   = u*Fs; 
ksi = ksi*N/Fs;
t = 0:N-1;

if oct == 0                     % Dirac
    signal = zeros(1,N);
    signal(u+1) = sign(cos(fi));

elseif oct == nextpow2(N)       % Fourier
    
    signal = cos(2*pi*ksi*t/N+fi);
    normSignal = sqrt(sum(signal.*signal));
    signal = signal/normSignal;

else                            % Gabor
    
    s = 2^oct;
    cosine_part = cos(2*pi*ksi*t/N+fi);
    
    if ~wrap
        exp_part = G((t-u)/s);
        signal = exp_part.*cosine_part;
        normSignal = sqrt(sum(signal.*signal));
        signal = signal/normSignal;
    else

        exp_part0 = G((t-u)/s);
        exp_part1 = G((t-N-u)/s);
        exp_part_1 = G((t+N-u)/s);

        exp_part = exp_part0 + exp_part1 + exp_part_1;
        signal = exp_part.*cosine_part;
        
        normSignal = sqrt(sum(signal.*signal));
        signal = signal/normSignal;
    end
end

% Decimate - note that here instead of summing the values (as is done for
% energy, we're just taking the mean. As a result, the energy of the
% decimated signal is less than that of the original signal
xs = t/Fs;
if decimateTime>1
    xs = ((decimateTime-1)/2:decimateTime:N)/Fs; xs=xs(1:N/decimateTime);
    signal = squeeze(mean(reshape(signal,[decimateTime N/decimateTime]),1));
end

end

function S = G(t)
    S = exp(-pi*t.*t);  % The constant 2^(1/4) is not needed because normalization takes care of it
end