% clc
% clear variables
% close all

%INPUTS: ECG, TIME, FREQ RATE
ECG = dataECG_cutted;
time = timeResampled ;
freq = 5000;
%% STEP 10: Calculate heart rate
% Calculate template and its correlation with ECG
[ECG_corr, temp, ratemean] = calc_template(ECG, freq, time);
% Detect R-peaks on ECG correlation
[ibi, indx, val, M] = detectpeaks(ECG, ECG_corr, time, freq, ratemean);

%% STEP 11: Check heartbeats
close all
figure
hist(ibi,20,'Color',[0.75 0 0.75])
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

% save peaks added and removed
ADD = []; % samples where heartbeats were added
REM = []; % samples where heartbeats were removed
peak_indx = zeros(1,length(ECG));
peak_indx(indx) = 1;

[ibi,indx,ECG,ECG_corr,peak_indx,ADD,REM] = editpeaks(ibi,indx, ECG, ECG_corr,peak_indx,time, freq,ADD, REM,M);
