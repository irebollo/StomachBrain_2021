function ibi = calc_heartrate(indx,time)

% This function generates the interbeat intervals matrix (HR) and the sample location of each R-peak
%
% inputs
% Peakn: logical array with 1 where is located the heartbeat
% time: time array of the dataset
%
% outputs
% HR: heart rate (interbeat intervals)
% indx: samples where are located the heartbeats
%
% Author: Diego Candia Rivera
% Email: diego.candia.r@ing.uchile.cl


ibi = [];
for i = 1:length(indx)-1 % to compute interbeat intervals is substracted the timestamps of consecutive heartbeats
    hr = time(indx(i+1)) - time(indx(i));
    ibi = [ibi hr];
end
ibi = ibi';

end