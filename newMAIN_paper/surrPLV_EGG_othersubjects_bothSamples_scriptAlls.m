% review_control_EGG_othersubjects



cfgMain = global_getcfgmain;
insideBrain = tools_getIndexBrain('inside');
subjectsNavigastric = global_subjectList_navi;
subjectsPhysiens = global_subjectList_phys;

%% First loop navigastric subjects
for iSubj = 1:length(subjectsNavigastric)
subj_idx = subjectsNavigastric(iSubj)
task = 'REST'
cfgMain.sample = 1
TrueSampleOfThisParticipant = cfgMain.sample;
rootDir= strcat(global_path2root(TrueSampleOfThisParticipant));
% subjects which are not this subject
surrogatePLVmatrix = zeros(length(subjectsNavigastric)+length(subjectsPhysiens)-1,length(insideBrain));
%% Load BOLD
% load BOLD of subject unfiltered
cfgMain.task = 'REST'
filename_bold_input = global_filename(subj_idx,cfgMain,'filename_csfr_Residuals_FB');
load(filename_bold_input)
%% Start loop randomization
subjListForAnalysis = subjectsNavigastric~= subj_idx
SurrogateSubjectList = subjectsNavigastric(subjListForAnalysis)
for iSurrogateSubject = 1:(length(subjectsNavigastric)-1)
    tic
 currentSurrogateSubject = SurrogateSubjectList(iSurrogateSubject)
% load EGG first other subject
EGGPhaseXVolumeFilename = global_filename(currentSurrogateSubject,cfgMain,strcat('EGGPhaseXVolumeFilename'));
load(EGGPhaseXVolumeFilename)
mostPowerfullFrequency = logEGGpreprocessing.mostPowerfullFrequency;
%Filter
centerFrequency = mostPowerfullFrequency; %
filter_frequency_spread=cfgMain.frequencySpread/1000; % In hz
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
filteredMRI=tools_bpFilter(error_csf_z,sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);
phaseMRI = hilbert(filteredMRI);
phaseMRI = angle(phaseMRI(cfgMain.beginCut:cfgMain.endCut,:)); % Cut data to have the same length as EGG (cut this way to get rid of fmri edge artifact on EGG)
% error_csf_z
clear filteredMRI
empPLV = abs (mean (exp (1i* (bsxfun (@minus , phaseMRI, angle (phaseXVolume)'))))); % get PLV
% clear other variables
clear  phaseMRI
% store PLV in surrogate PLV matrix
surrogatePLVmatrix(iSurrogateSubject,:) = empPLV;
toc
end
%% now loop through physiens subjects
cfgMain.sample = 2
SurrogateSubjectList = subjectsPhysiens
for iSurrogateSubject = 1:length(subjectsPhysiens)
    tic
  currentSurrogateSubject = SurrogateSubjectList(iSurrogateSubject)
% load EGG first other subject
EGGPhaseXVolumeFilename = global_filename(currentSurrogateSubject,cfgMain,strcat('EGGPhaseXVolumeFilename'));
load(EGGPhaseXVolumeFilename)
mostPowerfullFrequency = logEGGpreprocessing.mostPowerfullFrequency;
%Filter
centerFrequency = mostPowerfullFrequency; %
filter_frequency_spread=cfgMain.frequencySpread/1000; % In hz
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
filteredMRI=tools_bpFilter(error_csf_z,sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);
phaseMRI = hilbert(filteredMRI);
phaseMRI = angle(phaseMRI(cfgMain.beginCut:cfgMain.endCut,:)); % Cut data to have the same length as EGG (cut this way to get rid of fmri edge artifact on EGG)
clear filteredMRI
empPLV = abs (mean (exp (1i* (bsxfun (@minus , phaseMRI, angle (phaseXVolume)'))))); % get PLV
clear  phaseMRI
surrogatePLVmatrix(iSurrogateSubject+length(subjectsNavigastric),:) = empPLV;
toc
end
%% store
filename_allSurrogates = strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'ControlEGG_othersubjectBOTHSAMPLES',filesep,'distribution_sPLV_s',sprintf('%.2d',subj_idx),'.mat')
save(filename_allSurrogates,'surrogatePLVmatrix')
median_sPLV = median(surrogatePLVmatrix);
surrPLV = zeros (53,63,46); % empty 3d Volume for storing empirical PLV
surrPLV = surrPLV(:); % transformed into a vector
surrPLV(insideBrain) = median_sPLV; % get PLV
% bsxfun applies the operation in @minus from the vector of the third input
% to each column of the matrix of the second input
PLV3D = reshape(surrPLV,53,63,46); % reshape it from vector to matrix
filename_medianSurrogate  = strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'ControlEGG_othersubjectBOTHSAMPLES',filesep,'median_sPLV_s',sprintf('%.2d',subj_idx),'.nii')
tools_writeMri(PLV3D,filename_medianSurrogate)

end %loop navigastric

%% now physiens subjects


%% First loop navigastric subjects
for iSubj = 1:length(subjectsPhysiens)
subj_idx = subjectsPhysiens(iSubj)
task = 'REST'
cfgMain.sample = 2
TrueSampleOfThisParticipant = cfgMain.sample;
rootDir= strcat(global_path2root(TrueSampleOfThisParticipant));
% subjects which are not this subject
surrogatePLVmatrix = zeros(length(subjectsNavigastric)+length(subjectsPhysiens)-1,length(insideBrain));
%% Load BOLD
% load BOLD of subject unfiltered
cfgMain.task = 'REST'
filename_bold_input = global_filename(subj_idx,cfgMain,'filename_csfr_Residuals_FB');
load(filename_bold_input)
%% Start loop randomization
cfgMain.sample = 1
SurrogateSubjectList = subjectsNavigastric
for iSurrogateSubject = 1:(length(subjectsNavigastric)-1)
    tic
 currentSurrogateSubject = SurrogateSubjectList(iSurrogateSubject)
% load EGG first other subject
EGGPhaseXVolumeFilename = global_filename(currentSurrogateSubject,cfgMain,strcat('EGGPhaseXVolumeFilename'));
load(EGGPhaseXVolumeFilename)
mostPowerfullFrequency = logEGGpreprocessing.mostPowerfullFrequency;
%Filter
centerFrequency = mostPowerfullFrequency; %
filter_frequency_spread=cfgMain.frequencySpread/1000; % In hz
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
filteredMRI=tools_bpFilter(error_csf_z,sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);
phaseMRI = hilbert(filteredMRI);
phaseMRI = angle(phaseMRI(cfgMain.beginCut:cfgMain.endCut,:)); % Cut data to have the same length as EGG (cut this way to get rid of fmri edge artifact on EGG)
% error_csf_z
clear filteredMRI
empPLV = abs (mean (exp (1i* (bsxfun (@minus , phaseMRI, angle (phaseXVolume)'))))); % get PLV
% clear other variables
clear  phaseMRI
% store PLV in surrogate PLV matrix
surrogatePLVmatrix(iSurrogateSubject,:) = empPLV;
toc
end
%% now loop through physiens subjects
cfgMain.sample = 2
subjListForAnalysis = subjectsPhysiens~= subj_idx
SurrogateSubjectList = subjectsPhysiens(subjListForAnalysis)
for iSurrogateSubject = 1:length(subjectsPhysiens)-1
    tic
  currentSurrogateSubject = SurrogateSubjectList(iSurrogateSubject)
% load EGG first other subject
EGGPhaseXVolumeFilename = global_filename(currentSurrogateSubject,cfgMain,strcat('EGGPhaseXVolumeFilename'));
load(EGGPhaseXVolumeFilename)
mostPowerfullFrequency = logEGGpreprocessing.mostPowerfullFrequency;
%Filter
centerFrequency = mostPowerfullFrequency; %
filter_frequency_spread=cfgMain.frequencySpread/1000; % In hz
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
filteredMRI=tools_bpFilter(error_csf_z,sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);
phaseMRI = hilbert(filteredMRI);
phaseMRI = angle(phaseMRI(cfgMain.beginCut:cfgMain.endCut,:)); % Cut data to have the same length as EGG (cut this way to get rid of fmri edge artifact on EGG)
clear filteredMRI
empPLV = abs (mean (exp (1i* (bsxfun (@minus , phaseMRI, angle (phaseXVolume)'))))); % get PLV
clear  phaseMRI
surrogatePLVmatrix(iSurrogateSubject+length(subjectsNavigastric),:) = empPLV;
toc
end
%% store
filename_allSurrogates = strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'ControlEGG_othersubjectBOTHSAMPLES',filesep,'distribution_sPLV_s',sprintf('%.2d',subj_idx),'.mat')
save(filename_allSurrogates,'surrogatePLVmatrix')
median_sPLV = median(surrogatePLVmatrix);
surrPLV = zeros (53,63,46); % empty 3d Volume for storing empirical PLV
surrPLV = surrPLV(:); % transformed into a vector
surrPLV(insideBrain) = median_sPLV; % get PLV
% bsxfun applies the operation in @minus from the vector of the third input
% to each column of the matrix of the second input
PLV3D = reshape(surrPLV,53,63,46); % reshape it from vector to matrix
filename_medianSurrogate  = strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'ControlEGG_othersubjectBOTHSAMPLES',filesep,'median_sPLV_s',sprintf('%.2d',subj_idx),'.nii');
tools_writeMri(PLV3D,filename_medianSurrogate)

end %loop navigastric

