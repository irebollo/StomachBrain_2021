
%{
paper_2021_main


Scripts that reproduces the analysis on navigastric experiment for
publishing

0 - Set global variables
1 - MRI preprocessing with SPM
2 - preprocessing of EGG data
3 - timeseries analysis
4 - Group level statistics of Gastric-BOLD coupling
5- Figure 1
    a- Project gastric network to cortical surface (Bash & Matlab)
    b- Plot gastric network (Python)
    c- Quantify overlap with Resting state networks
    d- Effect sizes in resting state networks
    e- Effect sizes in Glasser parcellation regions

Ignacio Rebollo 30/04/2021
%}


%% 0 set global variables
cfgMain = global_getcfgmain;
cfgMain.task = 'REST'
subjectsPhysiens= [8 9 10 13 15 16 17 18 19 20 21 22 23 24 26 29 31 32 33 34 35 36 37 38 39 40 41 43 44] % 1st sample
subjectsNavigastric = [1:10 12 13 15:19 21 22 23 25 26 27  29 30 31 32 33 34 35 36 38 91 96] % 2nd sample

%% 1 - MRI preprocessing with SPM


%% 1 MRI preprocessing 


% slice timing correction,  normalization to MNI template, motion correction spatial smoothing

for iSubject=1:length(subjectsNavigastric)
    NAVIGASTRIC_group_preproc(subjectsNavigastric(iSubject),cfgMain.kernelWidth)
end


%% 2 - preprocessing of EGG data

cfgMain.sample=1 % First sample "Aka Navigastric"
subjects = subjectsPhysiens

for iSubj =1:length(subjects)
   suj = subjects(iSubj)
 prepro_egg(suj,cfgMain)
     press2continue = str2double(input('\n PRESS A KEY TO CONTINUE\n' ,'s'));

end

cfgMain.sample=2 % First sample "Aka PHYSIENS"
subjects = subjectsNavigastric
for iSubj =1:length(subjects)
   suj = subjects(iSubj)
 prepro_egg(suj,cfgMain)
     press2continue = str2double(input('\n PRESS A KEY TO CONTINUE\n' ,'s'));
end


%% 3 - Timeseries analysis

cfgMain.sample=2 % First sample "Aka PHYSIENS"
subjects = subjectsPhysiens

for iSubj =1:length(subjects)
subj_idx = subjects(iSubj)
timeseries_prepare_import2matlab(subj_idx,cfgMain) 
timeseries_preprocessBOLD(subj_idx,cfgMain) % preprocess the timeseries (remove polinomial and filter)
timeseries_csfSignal_obtainAndRegress(subj_idx,cfgMain) % obtain the wm and csf signals from the BOLD preprocessed timeseries and stores them
timeseries_preparePhases_Regression(subj_idx,cfgMain) % Filter and hilbert transform residuals of csf regression
timeseries_mapPLV_Regression(subj_idx,cfgMain) % Obtain PLV per voxel
end

cfgMain.sample=1 % First sample "Aka NAVIGASTRIC"
subjects = subjectsNavigastric

for iSubj =1:length(subjects)
subj_idx = subjects(iSubj)
timeseries_prepare_import2matlab(subj_idx,cfgMain) 
timeseries_preprocessBOLD(subj_idx,cfgMain) % preprocess the timeseries (remove polinomial and filter)
timeseries_csfSignal_obtainAndRegress(subj_idx,cfgMain) % obtain the wm and csf signals from the BOLD preprocessed timeseries and stores them
timeseries_preparePhases_Regression(subj_idx,cfgMain) % Filter and hilbert transform residuals of csf regression
timeseries_mapPLV_Regression(subj_idx,cfgMain) % Obtain PLV per voxel
end

% Obtain chain level PLV by computing EGG-BOLD coupling with the BOLD of
% one subjects and the EGG of all the other participants
surrPLV_EGG_othersubjects_bothSamples_scriptAlls

%% 4 - Group level statistics of Gastric-BOLD coupling

cfgMain = global_getcfgmain;
cfgMain.numberofrandomizations = 1000
cfgMain.clusterAlpha = 0.005; % exact one-sided p value that will be used for the cluster forming threshold 
cfgMain.alpha =0.05

timeseries_statsCluster_Regression_surrogateSubject_2samples(subjectsNavigastric,subjectsPhysiens,cfgMain)
 

%% 5a- Figure 1: Make coupling strenght images and compute effectsizes

% Make couping strenght images of all subjects for then computing coupling
% strength
makeCSimages_alls

% Obtain effect sizes (Matlab file)

effecSizes_bootstrap

% 5b- Project gastric network to cortical surface (Bash & Matlab)
% Project gastric network significant voxels
% It might be easier to run directly on the console
! /media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/standalone_scripts_for_MNI_fsaverage_projection/CBIG_RF_projectMNI2fsaverage.sh -s /media/irebollo/storage/projects/Navigastric/ClusterResults/All_63SA0250_final/SurrParBothSample_nR1000_CA0050_kw3_fir2_fspread_015_fOrder_5_tw_15csfr_mask.nii -o /media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/ -m /home/irebollo/matlab/bin/
! gzip -d /media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/lh.SurrParBothSample_nR1000_CA0050_kw3_fir2_fspread_015_fOrder_5_tw_15csfr_mask.FUSPROJ.nii.gz
! gzip -d /media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/rh.SurrParBothSample_nR1000_CA0050_kw3_fir2_fspread_015_fOrder_5_tw_15csfr_mask.FUSPROJ.nii.gz

% Note: I have modified the standalone_scripts_for_MNI_fsaverage_projection
% script to change the output filename from
% _RF_ANTs_MNI152_orig_to_fsaverage.nii to _FUSRPROJ. 

% Project effect sizes
! /media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/standalone_scripts_for_MNI_fsaverage_projection/CBIG_RF_projectMNI2fsaverage.sh -s /media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/cd.nii -o /media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/ -m /home/irebollo/matlab/bin/
! gzip -d /media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/rh.cd.FUSPROJ.nii.gz
! gzip -d /media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/lh.cd.FUSPROJ.nii.gz


% make figure 1a
% run in console/python
% python /media/irebollo/storage/projects/Navigastric/scripts_4_github/pysurf/plotCDGasnetRSN.py

%5c Effect sizes and overlap per RSN
EffectSizes_perRSN
OverlapRSNs

%5d bayes stats on group difference
bayes2groups

%5e Figure 1g 
EffectSizes_AllROIS

%control on piriform
control_piriform

% table
MakeTable3
%% 6- Figure 2: Position of gastric network in cortical gradients
% This barch resamples gastric network to margules voxel size

Figure_2_Gradients

%% 7 - Overlaps and effect size in each region 

% Figure 3
SpiderPlot_LateralCingulateSulci
% run this script in the console/python for the overlay in the brain surface 
    %python /media/irebollo/storage/projects/Navigastric/scripts_4_github/pysurf/plotGasnetPosterior_borders.py
% Figure 4
SpiderPlot_Operculum
% run this script in the console/python 
    %python /media/irebollo/storage/projects/Navigastric/scripts_4_github/pysurf/plotGasnetOperculum_borders
    
% figure 5
SpiderPlot_Posterior
% run this script in the console/python 
    %python /media/irebollo/storage/projects/Navigastric/scripts_4_github/pysurf/plotGasnetPosterior_borders.py  
    
% Figure 6
SpiderPlot_OtherRegions
% run this script in the console/python 
    %python /media/irebollo/storage/projects/Navigastric/scripts_4_github/pysurf/plotGasnetTransmodalRegions4PAPER.py

%% Figure 7 Demographics
%1-  Whole brain level demographics
%run r script
% /media/irebollo/storage/projects/Navigastric/scripts_4_github/Rfun/Script_demos.R

%% ROI level demographics
% First run
EffectSizes_AllROIS_export2R
% then run in R
% /media/irebollo/storage/projects/Navigastric/scripts_4_github/Rfun/Script_demos_Glasser.R

%%  voxel level demographics

% These batchs run the SPM second level models for the different
% explanatory variables

% Model including all the covariates without outliers
model101_BASIC_gender_FWD_BMI_Freq_Day_Group

% Model for stai
model102_BASIC_stai

% Model for EGG power
model103_BASIC_power

%Model for last meal
model104_BASIC_LastMeal

% LFHRV
model110_genderStaiBMIGroup_LFHRV

% HFHRV
model110_genderStaiBMIGroup_HFHRV

% HRRatio
model110_genderStaiBMIGroup_RATIOHRV

%% 

control_piriform
EffectSizes_AllROIS
EffectSizes_AllROIS_export2R
EffectSizes_perRSN
MakeTable3
OverlapRSNs
SpiderPlot_LateralCingulateSulci
SpiderPlot_Operculum
SpiderPlot_OtherRegions
SpiderPlot_Posterior
