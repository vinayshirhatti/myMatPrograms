%==========================================================================
% Vinay Shirhatti
% Time Frequency Analysis
% 
% Short-time Fourier transform
%==========================================================================

function [WindowStft,SPhalf,xTime] = stftTFA(data,lenWindow,trans_step)
%% Load the data and define the parameters
% load('mSTAData.mat'); % The signal data
% 
% lenWindow=16; % Window length
% 
% data = meanSTAVals(1,:);
% 
% trans_step = 2;

if ~exist('trans_step','var')
    trans_step = lenWindow/2;
else if (trans_step>lenWindow ||trans_step<1)
        trans_step = lenWindow/2;
    end
end

Fs = 2000; % Sampling frequency
N = size(data,2); % Length of the signal
F = Fs/N*(0:(N/2)-1); % Span of frequencies to be considered
xTime = (0:1/Fs:(N-1)/Fs);

window = hamming(lenWindow); % Define the window
% window = dpss(lenWindow,1,1);

figure; subplot(211); plot(data); axis('tight');
subplot(212); plot(window,'k');

%% STFT - Windowing and fft. And then sliding the window to next segment
lenW = length(window);
% trans_step=N/2;
for step=1:((N-(lenW))/trans_step)
    SlideWindow = zeros(N,1); % Initialize window as zeros of signal length
    SlideWindow((step-1)*trans_step+1:(step-1)*trans_step+lenW) = window; % Slide the hamming window to the appropriate section of the zero window
    WindowedSignal = SlideWindow'.*data; % Multiply signal with the window
    WindowStft(step,:) = fft(WindowedSignal); % Take fft
end

%% Combining the segment stft to compute the spectrogram and plotting it
SP = (abs(WindowStft)).^2;
SPhalf = SP(:,(1:N/2));
xlims=-200:((step+1)*(trans_step)/N)*400:200; ylims=F;
figure;
% imagesc(xlims,ylims',SPhalf');figure(gcf); % For plotting STFT
imagesc(xlims,ylims',log10(SPhalf'));figure(gcf); % For plotting log(STFT)
% pcolor(xlims,ylims,conv2Log(SPhalf')); shading interp;
set(gca,'YDir','normal');
colorbar;
title('STFT Spectrogram');
xlabel('Time (ms)');ylabel('Frequency (Hz)');
