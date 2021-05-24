function [HeartBeats, R_time] = heart_peak_detect(cfg,data)

% First input method. Simplified, Not recommended. See second input method below.
%
%[R_sample, R_time] = heart_peak_detection(ECG,fs)
%
% Detect R peaks in ECG data sampled at fs Hz.
% Inputs:
%       ECG         vector of ECG data
%       fs          sampling frequency
% Outputs:
%       R_sample    samples where beats have been detected.
%       R_time      time points where beats have been detected (assuming
%                   time starts at 0)
%
%
% Second input method:
%
% [HeartBeats] = heart_peak_detect(cfg)
% [HeartBeats] = heart_peak_detect(cfg,data)
%
% If you are calling heart_peak_detect with only the configuration as first
% input argument and the data still has to be read from file, you should
% specify
%   cfg.dataset      = string with the filename
%   cfg.channel      = the one channel to be read and taken as ECG channel
%
% If you are calling heart_peak_detect with also the second input argument
% "data", then that should either
%   - be a vector of ECG data points. first input cfg.fsample should
%   contain sampling frequency in this case.
%   - contain data that was read from file in a previous call to
%   FT_PREPROCESSING.
%
% In either case, the configuration options below apply.
%
%
% The options are specified with
%   - Preprocessing options:
%     cfg.hpfilter        = 'no' or 'yes'  highpass filter (default = 'yes')
%     cfg.hpfilttpye      = filter type (see ft_preprocessing, default = 'firws')
%     cfg.hpfreq          = highpass frequency in Hz (default = 1);
%     cfg.lpfilter        = 'no' or 'yes'  lowpass filter (default = 'yes')
%     cfg.lpfilttpye      = filter type (see ft_preprocessing, default = 'firws')
%     cfg.lpfreq          = lowpass frequency in Hz (default = 100);
%
%   - Algorithm options (see description below):
%     cfg.thresh          = z-threshold for 1st step detection of R-peaks (default = 10)
%     cfg.mindist         = minimum inter beat interval in seconds (IBI) (default = 0.35)
%     cfg.corthresh       = proportion of maximum correlation (default = 0.6)
%     cfg.PRmax           = maximum duration between P and R peaks in seconds (default = 0.25)
%     cfg.QRmax           = maximum duration between Q and R peaks in seconds (default = 0.05)
%     cfg.RSmax           = maxumum duration between R and S peaks in seconds (default = 0.1)
%     cfg.QTmax           = maximum QT interval in seconds (default = 0.42)
%
%   - Plotting options:
%     cfg.plotthresh      = open a figure to modify the threshold (default = 'yes')
%     cfg.plotbeat        = open a figure to show the average ECG around R peak (default = 'yes')
%     cfg.plotcorr        = open a figure to show the correlation level along the recording and enable editing threshold correlation (default = 'yes')
%     cfg.plotfinal       = open a figure to show final results with all peaks found (default = 'yes')
%     cfg.plotall         = whether to do all of the above (default = 'yes')
%
%   Algorithm:
%
%   The signal is first high and low pass filtered (default 1-100 Hz). The
%   square of the z-scored ECG is computed and a first detection is
%   performed as the peaks passing the cfg.thresh. Not all R peaks need to
%   be selected at this step. Just enough to create a template heart beat
%   ECG (HB) is necessary.
%   If cfg.plotthresh is true, then a figure is shown, allowing the user to
%   edit the threshold.
%   Then a template HB is computed (shown in a figure if cfg.plotbeat) and
%   convolved with the whole ECG time series. The resulting convolution is
%   normalized to have a maximum of 1 and beats are taken as peaks above
%   cfg.corthresh.
%   In both steps, a minimum distance between beats of cfg.mindist is
%   enforced.
%   Other peaks are found based on each R peak. Q is the minimum within 50
%   ms before R, S is the minimum within 100 ms after R, and T is the
%   maximum between the S peak and a maximum QT interval of 420 ms (a rough
%   standard...).
%

% v0 Maximilien Chaumon November 2016
% based on previously undocumented anonymous version
% v- Anne Buot Mars 2017
% modif visu outliers, ploting default option is yes
% debug addition/removal of new R-peaks: both the R-sample and R-value are updated
% debug update the IBI histogram after addition/removal of R-peaks

narginchk(1,2)
if isnumeric(cfg) && ~isempty(cfg) % first input method
    if not(isvector(cfg))
        error('ECG should be a one channel vector of data');
    end
    if not(exist('data','var'))
        error('Please provide sampling rate');
    end
    fs = data;
    data = [];
    data.fsample = fs;
    data.label = {'ECG'};
    data.trial = {cfg(:)'};
    data.time = {linspace(0,size(data.trial{1},2)/fs,size(data.trial{1},2))};
    cfg = [];
    cfg.structouput = 0;
elseif isstruct(cfg) || isempty(cfg) % second input method
    if nargin == 1
        data = ft_preprocessing(cfg);
    elseif nargin == 2
        if isnumeric(data)% data input as vector
            if ~isvector(data)
                error('ECG should be a one channel vector of data');
            end
            if ~isstruct(ft_checkopt(cfg,'fsample','double'))
                error('Please provide sampling frequency in cfg.fsample')
            end
            tmp = data;
            data = [];
            data.fsample = cfg.fsample;
            data.label = {'ECG'};
            data.trial = {tmp(:)'};
            data.time = {linspace(0,size(data.trial{1},2)/data.fsample,size(data.trial{1},2))};
            cfg = rmfield(cfg,'fsample');
        elseif isstruct(ft_checkdata(data,'datatype','raw'))
            if isfield(cfg,'channel')
                data = ft_preprocessing(cfg,data);
            end
            if not(isvector(data.trial{1}))
                first 
                error('Channel should be specified')
            end
        end
    end
else
    error('Bad input');
end

data_in = data;

def = [];
def.hpfilter        = 'yes';
def.hpfreq          = 1;    % low bound of high pass filter of ECG
def.hpfilttype      = 'firws';
def.lpfilter        = 'yes';
def.lpfilttype      = 'firws';
def.lpfreq          = 100;
def.thresh          = 10;    % z-threshold for 1st step detection of R-peaks
def.mindist         = 0.35; % minimum IBI
def.corthresh       = 0.5;  % proportion of maximum correlation
def.PRmax           = 0.25;
def.QRmax           = 0.05;
def.RSmax           = 0.1;
def.QTmax           = 0.42;
def.structoutput    = 1;
def.plotall         = 0;
def.plotthresh      = 1;
def.plotbeat        = 1;
def.plotcorr        = 1;
def.plotfinal       = 0;

cfg = setdef(cfg,def);

if cfg.lpfreq > data.fsample/2
    error(['Lowpass filter frequency too high. Set cfg.lpfreq below ' num2str(data.fsample/2)])
end
tmp             = [];
tmp.hpfilter    = cfg.hpfilter;
tmp.hpfreq      = cfg.hpfreq;
tmp.hpfilttype  = cfg.hpfilttype;
tmp.lpfilter    = cfg.lpfilter;
tmp.lpfreq      = cfg.lpfreq;
tmp.lpfilttype  = cfg.lpfilttype;
data            = ft_preprocessing(tmp,data);

ECG = data.trial{1};
time = data.time{1};

%% find R peaks

% standardize ecg
ECG2z = nanzscore(ECG).^2;

thresh = cfg.thresh; % default z-threshold
[R_sample, R_value] = peakseek(ECG2z, thresh, data.fsample .* cfg.mindist);
while istrue(cfg.plotthresh) || istrue(cfg.plotall)
    figure(47894);clf
    plot(ECG2z);
    hold on;
    scatter(R_sample,R_value);
    hline(thresh,':r');
    xlabel('samples');
    ylabel('zscore');
    zoom;
    
    title(sprintf(['Creating a template HB\nWe should have enough HB (don''t need all)\nclick to change threshold\nright-click to confirm'], thresh));
    
    [x,y,but] = ginput(1);
    if ~isempty(but) && but > 1
        break
    else
        thresh = y;
        [R_sample, R_value] = peakseek(ECG2z, thresh, data.fsample .* cfg.mindist);
    end
end


%% We now build the template and compute the correlation with the ecg channel

% build template based on identified r-peaks
HB_bound = round(.5 * data.fsample);
HB = NaN(numel(R_sample),2*HB_bound+1);
for ii=1:length(R_sample)
    if R_sample(ii)-HB_bound > 1 && R_sample(ii)+HB_bound < length(ECG)
        HB(ii,:) = ECG(R_sample(ii)-HB_bound:R_sample(ii)+HB_bound);
    end
end

mHB = nanmean(HB,1);

% we'll assume that signal at dt seconds before detected R-peaks should be
% lower than R-peaks. dt = 50 ms seems reasonable.
dt = round(.05 * data.fsample);
if sign(skewness(ECG)) == -1
    % flip
    ECG = -ECG;
    mHB = -mHB;
    HB = -HB;
end
while istrue(cfg.plotbeat) || istrue(cfg.plotall)
    figure(47894);clf
    t = linspace(-HB_bound,HB_bound,numel(mHB));
    hall = plot(t,HB','color',[.8 .8 .8]);
    hold on
    hm = plot(t,mHB);
    vline(0,'r');
    legend([hall(1),hm],{'individual ECG','mean ECG'},'location','best')
    axis tight
    title('r peak points up ? ''left mouse'' = no ''right mouse'' = yes');
    [x,y,but] = ginput(1);
    if ~isempty(but) && but > 1
        break
    else
        ECG = -ECG;
        mHB = -mHB;
        HB = -HB;
    end
end

ecg_n = ECG./max(ECG); % normalized ecg

% pad ecg
ecg_pad = [zeros(1,1000) ECG zeros(1,1000)];
cr = zeros(size(ecg_pad));

% compute correlation
for ii=1:length(ecg_pad)-length(mHB)
    cr(ii+round(length(mHB)/2)-1) = sum(ecg_pad(ii:ii+length(mHB)-1).*mHB);
end
cr = cr./max(cr); % normalize correlation to 1

% find peaks in correlation
thresh = cfg.corthresh;
[R_sample, R_value] = peakseek(cr(1001:end-1000), thresh, data.fsample .* cfg.mindist);

while istrue(cfg.plotcorr) || istrue(cfg.plotall)
    IBI_s = diff(R_sample)/data.fsample;
    plotIBI(ecg_n,cr,R_sample,R_value,thresh,IBI_s);
    s1 = questdlg('Would you like to change the threshold?');
    switch s1
        case 'Yes'
            [dum,thresh] = ginput(1);
            [R_sample, R_value] = peakseek(cr(1001:end-1000), thresh, data.fsample .* cfg.mindist);
        case 'No'
            break
        case 'Cancel'
            error('Don''t click cancel unless you want to cancel...')
    end
end
while istrue(cfg.plotcorr) || istrue(cfg.plotall)
    IBI_s = diff(R_sample)/data.fsample;
    plotIBI(ecg_n,cr,R_sample,R_value,thresh,IBI_s);
    
    s2 = questdlg('Are there still outliers present? ');
    switch s2
        case 'Yes'
            axes(findobj(gcf,'tag','hhist'))
            title('right click to define threshold for low IBI')
            [lw,~] = ginput(1);
            title('right click to define threshold for high IBI')
            [hg,~] = ginput(1);
            % find outlier peaks
            IBI_out = IBI_s < lw | IBI_s > hg;
            for ii = R_sample(IBI_out)
                while istrue(cfg.plotcorr) || istrue(cfg.plotall)
                    h = figure(478490);clf;
                    plot(ecg_n);
                    hold on
                    %                     plot(cr(1001:end-1000),'r')
                    scatter(R_sample,ecg_n(R_sample))
                    xlim([ii-10000 ii+10000]);
                    xlabel('samples')
                    ylabel('a.u.')
                    title('Add with left mouse. Remove with right mouse. Move on with Escape.')
                    [x,y,but] = ginput(1);
                    x = round(x);
                    if but == 27
                        break
                    elseif but == 1
                        % find closest maximum
                        X = [1:numel(ecg_n);ecg_n*numel(ecg_n)]';
                        [closest] = dsearchn(X,[x y*numel(ecg_n)]);
                        scatter(closest,ecg_n(closest))
                        R_sample(end+1) = closest; [R_sample,indR] = sort(R_sample);
                        R_value(end+1) = ecg_n(closest); R_value = R_value(indR);
                        if ~isequal(length(indR),length(R_value))
                            error('Not the same number of elements in each vectors')
                        end
                    else
                        % find closest detected peak
                        [dum,closest] = min(abs(R_sample-x));
                        R_sample(closest) = [];
                        R_value(closest) = [];
                    end
                end
            end
            try
                close(h)
            end
        otherwise
            break
    end
end

if cfg.structoutput
    %% find Q S T
    Q_sample = NaN(size(R_sample));S_sample = NaN(size(R_sample));T_sample = NaN(size(R_sample));
    for i_R = 1:numel(R_sample)
        % Q
        idx = max(1,round(R_sample(i_R) - cfg.QRmax * data.fsample)):R_sample(i_R);
        ECGtmp = ECG(idx);
        [v,p] = min(ECGtmp);
        Q_sample(i_R) = p + idx(1) - 1;
        % P
        idx = max(1,round(R_sample(i_R) - cfg.PRmax * data.fsample)):Q_sample(i_R);
        ECGtmp = ECG(idx);
        [v,p] = max(ECGtmp);
        P_sample(i_R) = p + idx(1) - 1;
        % S
        idx = R_sample(i_R):min(numel(ECG),round(R_sample(i_R) +  cfg.RSmax * data.fsample));
        ECGtmp = ECG(idx);
        [v,p] = min(ECGtmp);
        S_sample(i_R) = p + idx(1) - 1;
        % T
        idx = S_sample(i_R):min(numel(ECG),round(Q_sample(i_R)+cfg.QTmax*data.fsample));
        ECGtmp = ECG(idx);
        [v,p] = max(ECGtmp);
        T_sample(i_R) = p + idx(1) - 1;
    end
    P_time = skipnans(data_in.time{1},P_sample);
    Q_time = skipnans(data_in.time{1},Q_sample);
    R_time = skipnans(data_in.time{1},R_sample);
    S_time = skipnans(data_in.time{1},S_sample);
    T_time = skipnans(data_in.time{1},T_sample);
    
    for i = 1:numel(R_sample)
        HeartBeats(i).P_sample  = P_sample(i);
        HeartBeats(i).P_time    = P_time(i);
        HeartBeats(i).Q_sample  = Q_sample(i);
        HeartBeats(i).Q_time    = Q_time(i);
        HeartBeats(i).R_sample  = R_sample(i);
        HeartBeats(i).R_time    = R_time(i);
        HeartBeats(i).S_sample  = S_sample(i);
        HeartBeats(i).S_time    = S_time(i);
        HeartBeats(i).T_sample  = T_sample(i);
        HeartBeats(i).T_time    = T_time(i);
    end
    if istrue(cfg.plotfinal) || istrue(cfg.plotall)
        figure;
        hold on
        set(gca,'colororder',[0         0    1.0000
            0    0.5000         0
            0.8000         0         0
            0    0.7500    0.7500
            0.7500         0    0.7500]);
        
        plot(data_in.time{1},data_in.trial{1},'k');
        
        todel = isnan(P_sample);
        toplot = [P_time(~todel); data_in.trial{1}(P_sample(~todel))];
        plot(toplot(1,:),toplot(2,:),'.','markersize',10);
        todel = isnan(Q_sample);
        toplot = [Q_time(~todel); data_in.trial{1}(Q_sample(~todel))];
        plot(toplot(1,:),toplot(2,:),'.','markersize',10);
        todel = isnan(R_sample);
        toplot = [R_time(~todel); data_in.trial{1}(R_sample(~todel))];
        plot(toplot(1,:),toplot(2,:),'.','markersize',10);
        todel = isnan(S_sample);
        toplot = [S_time(~todel); data_in.trial{1}(S_sample(~todel))];
        plot(toplot(1,:),toplot(2,:),'.','markersize',10);
        todel = isnan(T_sample);
        toplot = [T_time(~todel); data_in.trial{1}(T_sample(~todel))];
        plot(toplot(1,:),toplot(2,:),'.','markersize',10);
        legend({'ECG','P','Q','R','S','T'});
    end
    clear beats_time
else
    HeartBeats = R_sample;
end
drawnow


















