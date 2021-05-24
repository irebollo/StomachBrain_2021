function sk = skewness(data)
x0 = data - nanmean(data);
s2 = nanmean(data.^2);
m3 = nanmean(data.^3);
sk = m3 ./ s2 .^(1.5);
