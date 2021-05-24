function [t,ibi,t_int,ibi_int, F, PSD] = IBI_PSD(hb_sample,fs,fs_heart,window,overlap,nfft)

% INPUTS
% hb_samples: samples where are located the R peaks
% fs: sample rate of the array with R peaks
% fs_heart: sample rate to interpolate IBI - USE MINIMUM 1 HZ
% window: size in time to compute PSD
% overlap: window overlap in %
% nfft: spectral resolution 

% OUTPUTS
% t: timestamps of IBI (x-axis)
% ibi: interbeat intervals (y axis)
% t_int: interpolated timestamps using fs_heart
% ibi_int: interpolated ibi
% F: frequency stamps spectrum
% PSD: frequency spectrm

% EXAMPLE
% [t,ibi,t_int,ibi_int, F, PSD] = IBI_PSD(hb_sample,fs,1,100,50,1024)
%% compute IBI
% hb_sample = Rpeak_sample

hb_time = hb_sample/fs;
ibi = [];
for i = 1:length(hb_time)-1 % to compute interbeat intervals is substracted the timestamps of consecutive heartbeats
    hr = hb_time(i+1) - hb_time(i);
    ibi = [ibi hr];
end

%% compute timestamps of IBI
% starting at previous heartbeat
for m=1:length(ibi)
    t(m) = sum(ibi(1:m));
end

% starting at 0
% t(1) = 0; 
% for m=2:length(ibi)
%     t(m) = sum(ibi(1:m));
% end

%% linspace timestamps
t_int = t(1):1/fs_heart:t(end);

%% interpolated ibi
ibi_int = interp1(t,ibi,t_int,'spline')'; 

%% PSD
samplewindow = round(window*fs_heart);
noverlap = round(samplewindow*overlap/100);
ibi_int=detrend(ibi_int,'linear'); %optional
% ibi_int=ibi_int-mean(ini_int); %optional
[PSD,F] = pwelch(t_int,window,noverlap,(nfft*2)-1,fs_heart,'onesided');
end

