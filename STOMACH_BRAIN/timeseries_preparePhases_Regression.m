function timeseries_preparePhases_Regression(subj_idx,cfgMain)

%{

Load EGG and FMRI timeseries, filter fMRI timeseries at gastric peak, 
extract phases of fmri timeseries(CSFregressed) and save them into
disk in timeseries folder

inputs:
subj_idx = s number
cfgMain must contain fields
    kernelWidth,Timeseries2Regress,frequencySpread ,fOrder,cfgMain.beginCut,cfgMain.endCut

kernelWidth: with of the smoothing kernel from preprocessing, paper  = 3mm

cfgMain.Timeseries2Regress should be 'csf' to load residuals of csf regression
fOrder : mutiplicative factor for the order of the filter
frequencySpread: spead of the time domain filter in hz * 1000, paper = 0.015 hz = 15,

begin and end cut are the voulmes that are discarded to avoid the filter
ringing artifact

cfgMain.transitionWidth is the transition width of the filter, paper is 15
offset is with respect to EGG peaking filter, only for control analysis.
offset is in hz x 1000 e.g. and offset of 0.006 hz is a value of 6

File input: Fullband bold timeseries
Y:\Subjects\Subject13\Timeseries\MRItimeseries\csfRegressionResiduals_FB_S_13_kw3

Output: saves Phase and Amplitude data separatly in subject timeseries folder
Y:\Subjects\Subject13\Timeseries\MRItimeseries\csfResiduals_FB_phases_s13_kw3_fir2_fspread_015

IR COMMENTED 28/06/2017

%}
%% Pass parameters of cfgMain to function

cfgMain.endCut = cfgMain.endCut;
% randomized  = cfgMain.randomized;
saveFiltered = cfgMain.saveFiltered;
offset = cfgMain.offset;

%% output filename

disp('+++++++++++++++++++++++++++++++ loading data')

% inputs name
% BOLDTimeseriesFilename = global_filename(subj_idx,cfgMain,strcat('filename_',cfgMain.Timeseries2Regress,'_Residuals_FB'));
BOLDTimeseriesFilename = global_filename(subj_idx,cfgMain,strcat('filename_',cfgMain.Timeseries2Regress,'_Residuals_FB'));
% BOLDTimeseriesFilename = global_filename(subj_idx,cfgMain,'BOLD_filtered_fullbandFilename');

timeseries = load(BOLDTimeseriesFilename);

% Load the information about the peaks of the EGG


if strcmp(cfgMain.task,'REST')

EGGPhaseXVolumeFilename = global_filename(subj_idx,cfgMain,'EGGPhaseXVolumeFilename');
load(EGGPhaseXVolumeFilename)
else
    EGGPhaseXVolumeFilename = strcat(global_path2subject(subj_idx,cfgMain.sample),'timeseries',filesep,'EGGTimeseries',filesep,cfgMain.task,filesep,...
    'EGG_data_S_',sprintf('%.2d',subj_idx));
load(EGGPhaseXVolumeFilename)
logEGGpreprocessing.mostPowerfullFrequency = EGG_data.peakFrequency;
end



mostPowerfullFrequency = logEGGpreprocessing.mostPowerfullFrequency;


% nVolumes = size(timeseries.error_csf_z,1)

% Output
BOLDPhasesTimeseriesFilename =  global_filename(subj_idx,cfgMain,strcat('filename_',cfgMain.Timeseries2Regress,'_Residuals_FB_phases'));

plotDir = strcat (global_path2subject(subj_idx,cfgMain.sample),'PreprocessingLog',filesep);
plotFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_',cfgMain.task,'_MRIPhases');

%% Filter-hilbert

disp('+++++++++++++++++++++++++++++++ FILTER HILBERT')

centerFrequency = mostPowerfullFrequency; %
filter_frequency_spread=cfgMain.frequencySpread/1000; % In hz
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units

filteredMRI=tools_bpFilter(eval(strcat('timeseries.',cell2mat(fields(timeseries)))),sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);

% filteredMRI=tools_bpFilter(timeseries.BOLD_filtered_zscored',sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);

phaseMRI = hilbert(filteredMRI);


if strcmp (cfgMain.task,'REST')
    nVolumes = 450
else
    nVolumes = EGG_data.nVolumes;
end
    
phaseMRI = phaseMRI(cfgMain.beginCut:nVolumes-15,:); % Cut data to have the same length as EGG (cut this way to get rid of fmri edge artifact on EGG)
filteredMRI = filteredMRI(cfgMain.beginCut:nVolumes-15,:);



%% Saving

disp('+++++++++++++++++++++++++++++++ SAVING')

save(BOLDPhasesTimeseriesFilename,'phaseMRI','-v7.3')

%% Sanity plot
if cfgMain.savePlots == 1
    insideBrain = tools_getIndexBrain('inside');
    
%     voxelCoordinates = sub2ind([79,95,79],9,45,53); % voxel in somatomotor cortex

    voxelCoordinates = sub2ind([53,63,52],11,30,37); % Right somatosensory cortex
    voxelCoordinates_inside = zeros(53*63*52,1);
    voxelCoordinates_inside(voxelCoordinates)=1;
    voxelCoordinates_inside = voxelCoordinates_inside(insideBrain);
    ind_voxelCoordinates_inside = find(voxelCoordinates_inside);
    
    if cfgMain.plotFigures == 0;
        SanityPlot = figure('visible','off');
    else
        SanityPlot = figure('visible','on');
    end
    
    
    hold on
    subplot(3,1,1)
    plot(timeseries.error_csf_z(:,ind_voxelCoordinates_inside),'r','LineWidth',2)
%     plot(timeseries.BOLD_filtered_zscored(ind_voxelCoordinates_inside,:),'r','LineWidth',2)
    title(['S',sprintf('%.2d',subj_idx),32,'fullband CSF timeseries in rSS'],'fontsize',18)
    grid on
    subplot(3,1,2)
    plot(cfgMain.beginCut:nVolumes-15,filteredMRI(:,ind_voxelCoordinates_inside),'r','LineWidth',4)
    title(['S',sprintf('%.2d',subj_idx),32,'rSS bandpassfiltered at' 32 num2str(mostPowerfullFrequency)],'fontsize',18)
    grid on
    subplot(3,1,3)
    plot(cfgMain.beginCut:nVolumes-15,angle(phaseMRI(:,ind_voxelCoordinates_inside)),'r','LineWidth',4)
    title(['S',sprintf('%.2d',subj_idx),32,' phases rSS bandpassfiltered at' 32 num2str(mostPowerfullFrequency)],'fontsize',18)
    grid on
    
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    
    print ('-dpng', '-painters', eval('plotFilename'))
    print ('-depsc2', '-painters', eval('plotFilename'))
    saveas(SanityPlot,strcat(plotFilename,'.fig'))
end



end