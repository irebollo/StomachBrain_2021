function [ECG_corr, temp, ratemean] = calc_template(ECG, freq, time)
%CALC_TEMPLATE Calcuates a template of an r-peak based on an ecg signal.
%   Requests from the user to select a clean time window within which
%   r-peaks can be easily selected based on an amplitude-threshold and
%   computes a template, which is an average of the r-peaks detected based
%   on the threshold. Returns this template, which can be used to find
%   r-peaks in the complete ecg.

% INPUT:
% ecg - the complete ecg data

% OUTPUT:
% temp - the r-peak template to use in the r-peak detection precedure

% References: MODIFIED FROM N. Wolpert and C. Richter CODES.

% standardize ecg
ecg_z = zscore(ECG.^3);

% instruct the user to select the on- and offset of a clean time window
% within which r-peaks are clearly detectable

figure
plot(ecg_z);
set(gcf,'units','points','position',[10,10,1200,300])
title('z-transformed ecg data')

fprintf('Please zoom in to find a suitable time window containing clean r-peaks\n')
s1 = input('Press <y> if you found a suitable time window: ', 's');

if strcmp(s1,'y')
    title('Select time window ONSET with mouse')
    [onset,~] = ginput(1);
    title('Select time window OFFSET with mouse')
    [offset,~] = ginput(1);
end

% extract the ecg data that correspond to the time window selected by the
% user:
% look for the nearest timestamp in the ecg data structure that matches the
% on- and offset (converted into seconds)
ecg_tw = ECG(onset:offset);

% adapt an amplitude threshold to pick the r-peaks in the ecg time window
thresh = mean(ecg_tw) + 2.5*std(ecg_tw); % default z-threshold
s1 = 'y';
s0 = 'y';
while strcmp(s1,'y')
    % detect peaks above threshold separated by at least 0.35 seconds
    [p, v] = peakdetect2(ecg_tw, thresh, freq.*0.35);
    
    close
    figure
    plot(ecg_tw);
    hold on
    plot(p,v,'ok')
    axis tight
    line(get(gca,'xlim'),[thresh thresh],'linestyle','--','color','r')
    xlabel('samples')
    ylabel('z-score')
    zoom
       
    s1 = input('Would you like to change the threshold (y/n)? ', 's');
    if strcmp(s1,'y')
        thresh = input('New threshold (z-score): ');
    end
    close
end

rate = [];
for i = 1:length(p)-1
    rate(i) = (p(i+1)-p(i))/freq;
end
ratemean = mean(rate);

% build template based on identified r-peaks
temp=[];
% if only one peak has been detected, the template is just the peak
if length(p)==1
    temp=ecg_tw;
% in alll other cases, it is the mean of the peaks detected
else
    for ii=1:length(p)
        if p(ii)-200 > 1 && p(ii)+200 < length(ecg_tw)
            temp(ii,:) = ecg_tw(p(ii)-200:p(ii)+200);
        end
    end

    temp = mean(temp);
end

% Modify polarity if needed
if sign(skewness(temp)) == -1
    ecg_tw = -ecg_tw;
    temp = -temp;
end

%% Correlation template and ECG
[ECG, ECG_corr] = calc_correlation(ECG, temp, time, freq);

s1 = input('Would you like to modify template (y/n)? ', 's');
    if strcmp(s1,'y')
        [ECG_corr, temp, rate] = calc_template(ECG, freq, time);
    end

end
