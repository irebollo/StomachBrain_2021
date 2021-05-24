function [ibi, indx, v, M] = detectpeaks(ECG, ECG_corr, time, freq, ratemean)
% This function does the R-peaks detection using a sliding time window. 
% The time window length has to be entered in the command window when is required
%
% inputs
% ECGnew: struct with the EKG
% time: time array of the whole dataset
% freq: sample frequency
%
% outputs
% HR: heart rate (interbeat intervals)
% indx: samples where are located the heartbeats
% ECGnew: struct with the EKG
% peak_indx: logical array with 1 where is located the heartbeat
% M: max amplitude uV
%
% Author: Diego Candia Rivera
% Email: diego.candia.r@ing.uchile.cl

M = 2*max(ECG_corr);
ECGm = ECG_corr;
midwind = round(freq*(ratemean+0.2)/2);
wind = 2*midwind;


delta = 5;
% find local maxima
for i = 1:length(ECGm)
    if i - midwind <= 1 % this takes care of the start of the window
        o1 = i;
        o2 = o1 + wind;
        ECG_window = ECGm(o1:o2);
        [~, Mi] = max(ECG_window);
        j = Mi(1) + o1 - 1;
        if j - delta <= 0 % this takes care of the start of the window
            [~,mx] = max(ECG(1:j+delta));  
        else
            [~,mx] = max(ECG(j-delta:j+delta));
            mx = mx(1) + j - delta - 1;
        end        
        [~,mx] = max(ECG(1:j+delta));
    elseif i + midwind >= length(ECG) % this takes care of the end of the window
        o1 = i - midwind ;
        o2 = length(ECG);
        ECG_window = ECGm(o1:o2);
        [~, Mi] = max(ECG_window);
        j = Mi(1) + o1 - 1;
        if j + delta > length(ECG) % this takes care of the end of the window
            [~,mx] = max(ECG(j-delta:end));  
        else
            [~,mx] = max(ECG(j-delta:j+delta));
            mx = mx(1) + j - delta - 1;
        end        
    else
        o1 = i-midwind;
        o2 = i+midwind;
        ECG_window = ECGm(o1:o2);
        [~, Mi] = max(ECG_window);
        j = Mi(1) + o1 -1;
        [~,mx] = max(ECG(j-delta:j+delta));
        mx = mx(1) + j - delta - 1;
    end
    ECGm(mx) = M;
end

% store indexes of possible R peaks
indx = [];
for i = 1:length(ECGm)-1
    if ECGm(i) == M
        indx = [indx i];
    end
end

% discard consecutive peaks
m = ECGm(indx);
d = diff(indx);
jump = (d> (ratemean-0.2)*freq);
p = [];

sect=1;
while sect<=length(d)
    if jump(sect)
        p = [p indx(sect)];
    else
        s = [];
        while ~jump(sect) & sect<length(d)
            s = [s sect];
            sect = sect + 1;
        end
        [lm, li] = max(m(s));
        p = [p indx(s(li))];
    end
    sect = sect+1;
end
indx = p;

v = ECG_corr(indx);
ibi = calc_heartrate(indx,time);


figure
hist(ibi,20)
title('Heart-Rate Histogram')
xlabel('Time [s]')
ylabel('Frequency')
figure
ECGplot1 = ECG_corr*15;
ECGplot1(ECGplot1>5) = 5;
plot(ECGplot1)
hold on
plot(indx,ECGplot1(indx),'ok')
hold on
ECGplot2 = zscore(ECG)+10;
ECGplot2(ECGplot2>15) = 15;
plot(ECGplot2)
hold on
plot(indx,ECGplot2(indx),'ok')
title('EEG signal with ECG artifacts')
xlabel('[samples]')
ylabel('[uV]')
ylim([0 15])
set(gcf,'units','points','position',[10,10,1200,300])
