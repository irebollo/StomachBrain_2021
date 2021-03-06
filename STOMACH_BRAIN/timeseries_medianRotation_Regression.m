function timeseries_medianRotation_Regression(subj_idx,cfgMain)
%{

Stores a .nii image of the median PLV obtained with a timeshifted EGG signals (concatenated and rotated).
It does all 360 possible rotations of EGG signal and calculate PLV at each
rotations and takes the median PLV obtained across timeshift across
voxels to be used as surrogate PLV in group level statistics


inputs:
subj_idx = s number
cfgMainMain must contain fields
    kernelWidth,Timeseries2Regress,frequencySpread ,fOrder,cfgMainMain.beginCut,cfgMainMain.endCut

kernelWidth: with of the smoothing kernel from preprocessing, paper  = 3mm

cfgMainMain.Timeseries2Regress should be 'csf' to load residuals of csf regression
fOrder : mutiplicative factor for the order of the filter
frequencySpread: spead of the time domain filter in hz * 1000, paper = 0.015 hz = 15,

begin and end cut are the voulmes that are discarded to avoid the filter
ringing artifact

cfgMainMain.transitionWidth is the transition width of the filter, paper is 15
offset is with respect to EGG peaking filter, only for control analysis.
offset is in hz x 1000 e.g. and offset of 0.006 hz is a value of 6

Input BOLD timeseries
Y:\Subjects\Subject13\Timeseries\MRItimeseries\csfRegressionResiduals_FB_S_13_kw3

% Output: saves data in subject timeseries folder as a 3D .nii image
Y:\Subjects\Subject13\Timeseries\PhasesAnalysis\medianRotation_csfr_S_13_kw3_fir2_fspread_015_fOrder_5_tw_15

IR commented 28/06/2017

%}
%% Import cfgMain parameters
fOrder = cfgMain.fOrder;
frequencySpread = cfgMain.frequencySpread;
kernelWidth= cfgMain.kernelWidth;
offset = cfgMain.offset;

outsideBrain = tools_getIndexBrain('outside');
insideBrain = tools_getIndexBrain('inside');

plotDir = strcat (global_path2subject(subj_idx,cfgMain.sample),'PreprocessingLog',filesep);
plotFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_chancePLVhistogram');
%% Load data and outputfilename

medianRotationFilename = strcat(global_filename(subj_idx,cfgMain,strcat('medianRotationFilename_',cfgMain.Timeseries2Regress))); % output
BOLDPhasesTimeseriesFilename = global_filename(subj_idx,cfgMain,strcat('filename_',cfgMain.Timeseries2Regress,'_Residuals_FB_phases')); % input filename

load (BOLDPhasesTimeseriesFilename)

SubjectDataRoot = strcat(global_path2subject(subj_idx,cfgMain.sample),'Timeseries',filesep);
EGGPhaseXVolumeFilename = global_filename(subj_idx,cfgMain,'EGGPhaseXVolumeFilename');
load(EGGPhaseXVolumeFilename)
phaseXVolume = phaseXVolume(1:420);
phaseMRI = phaseMRI(1:420,:);

% SubjectDataRoot = strcat(global_path2subject(subj_idx),'Timeseries',filesep);
% filenameArtDet =strcat(SubjectDataRoot,'EGGTimeseries',filesep,...
%     'ArtefactDetection_S_',sprintf('%.2d',subj_idx));
% load(filenameArtDet)


% if ~isempty(artdefDownsampled)
% for iEvent=1:size(artdefDownsampled,1)
% phaseMRI(artdefDownsampled(iEvent,3):artdefDownsampled(iEvent,4),:) = [];
% phaseXVolume(artdefDownsampled(iEvent,3):artdefDownsampled(iEvent,4)) = [];
% end
% end


%% Rotate EGG

indexRotations=31:390; % Rotating at least two minute (30 TR = 60s) at the beggining or end
rotatedPhaseEGG = zeros(length(indexRotations),length(phaseXVolume));
for iRotation = 1 : length(indexRotations)
rotatedPhaseEGG(iRotation,:) = circshift(phaseXVolume,[0 indexRotations(iRotation)]);
end


%% Calculate PLV
disp('+++++++++++++++++++++++++++++++ RPLV')

% initialize structure for distribution of rotated PLV values
RPLV = zeros(length(indexRotations),53*63*46); % R from rotated
%iterate through all rotation and calculate PLV

for iRotation = 1 : length(indexRotations)
    

PLV = zeros (53,63,46);
PLV = PLV(:);

currentPhaseEGG = rotatedPhaseEGG (iRotation,:) ;
phaseDifference = bsxfun (@minus , angle(phaseMRI), angle(currentPhaseEGG )');
PLV(insideBrain) =    abs (mean (exp (1i* phaseDifference ) ) ); % 

RPLV(iRotation,:) = PLV;% timeseries_get_PLV(phaseMRI,rotatedPhaseEGG(iRotation,:)'); 

disp('Rotation number for subject:')
disp(iRotation)
disp (subj_idx)
end

%% get and save median rotation

medianPLV= zeros(1,length(RPLV));
medianPLV = median(RPLV,1);

medianPLV(outsideBrain) = 0;
medianPLV = reshape (medianPLV,53,63,46);

tools_writeMri(medianPLV,medianRotationFilename)

%% SanityCheck : % check if value of rotation of 

if cfgMain.savePlots == 1
    
%     voxelCoordinates = sub2ind([79,95,79],9,45,53); % voxel in somatomotor cortex
    voxelCoordinates = sub2ind([53,63,46],11,30,37);

    voxelCoordinates_inside = zeros(53*63*46,1);
voxelCoordinates_inside(voxelCoordinates)=1;
voxelCoordinates_inside = voxelCoordinates_inside(insideBrain);
ind_voxelCoordinates_inside = find(voxelCoordinates_inside);

if cfgMain.plotFigures == 0;
    SanityPlot = figure('visible','off');
else
    SanityPlot = figure('visible','on');
end

% Plot histogram of PLV across the brain

nhist(medianPLV(insideBrain))
xlabel('PLV')
title(['S',sprintf('%.2d',subj_idx),32,'surrogatePLV across bain. Mean:' num2str(mean(medianPLV(insideBrain))) ' rSS voxel:' 32 num2str(medianPLV(insideBrain(ind_voxelCoordinates_inside)))],'fontsize',18)


set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto');

print ('-dpng', '-painters', eval('plotFilename'))
print ('-depsc2', '-painters', eval('plotFilename'))
saveas(SanityPlot,strcat(plotFilename,'.fig'))    
end

end