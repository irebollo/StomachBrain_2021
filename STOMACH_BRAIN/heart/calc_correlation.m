function [ECG, ECG_corr] = calc_correlation(ECG, temp, time, freq)

ECG_pad = [zeros(1,1000) ECG zeros(1,1000)];
cr = zeros(size(ECG_pad));
for i=1:length(ECG_pad)-length(temp)
    cr(i+round(length(temp)/2)-1) = sum(ECG_pad(i:i+length(temp)-1).*temp);
end
ECG_corr = cr(1001:end-1000)/max(cr); % normalize correlation to 1

% optional normalization
ECG_corr = zscore(ECG_corr.^3);
ECG_corr(ECG_corr < 0) = 0;
ECG_corr(ECG_corr > 10) = 10;
ECG_corr = (ECG_corr - min(ECG_corr) )/(max(ECG_corr) - min(ECG_corr));

% plot to evaluate
figure
ECGplot1 = ECG_corr*15;
ECGplot1(ECGplot1>5) = 5;
plot(ECGplot1)
hold on
ECGplot2 = zscore(ECG)+10;
ECGplot2(ECGplot2>15) = 15;
plot(ECGplot2)
title('ICA and Correlation')
xlabel('[samples]')
ylabel('[uV]')
xlim([1 50000])
ylim([0 15])
set(gcf,'units','points','position',[10,10,1200,300])

end

