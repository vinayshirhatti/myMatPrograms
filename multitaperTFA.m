% MULTITAPER CODE
% 
%==========================================================================

function [S S1] = multitaperTFA(data)

% load('mSTAData.mat');
% 
% data = meanSTAVals(4,:);

params.Fs=2000;
params.tapers=[5 9];
params.trialave=1;
params.err=0;
params.pad=0;

[S,f]=mtspectrumc(data,params);

figure;
% subplot(121);
plot_vector(S,f);


movingwin=[0.015 0.003];

[S1,t,f]=mtspecgramc(data,movingwin,params);
% subplot(122)
figure;
% t = -200:400/(length(data)-1):200;
plot_matrix(S1,t,f);
% xlim([-200:400/(length(data)-1):200]);
% set(gca,'XTickLabel',[-200;100;200])
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('MT Spectrogram');
% caxis([8 28]); 
colorbar;
% plot_matrix(S1,t,f);
end

