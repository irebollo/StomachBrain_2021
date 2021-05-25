# StomachBrain_2021
Scripts for reproducing the results of Rebollo &amp; Tallon-Baudry 2021



The scripts are organized into folder-specific subroutines. 
The main script to run all the analysis is paper_main_2021.m, and is located in the root folder. This scripts contain step by step procedure to reproduce the results. The scripts include the filenames to to python, r or bash scripts in their corresponding order

The rest of the scripts are organized into folders. Matlab scripts for the preprocesing and the phase timeseries analysis 
(From EGG peak channel selection up to group level statistical analysis), are located in STOMACH_BRAIN. These scripts are based on our previous repository from Rebollo et al 2018 Elife paper, and further documentation can be found here https://github.com/irebollo/stomach_brain_Scripts

The folder newMAIN_paper contains all the analysis specific to this article, including the comparison with cortical gradients, the effect sizes in Glasser regions.
 
The folder bashcripts contains a slightly modified version of the registration fusion stand alone scripts (https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/registration/Wu2017_RegistrationFusion). However, the command line to convert files from volume to surface can be found in paper_main_2021.m

The folder pysurf contains the plotting routines to produce the brain surface figures of the paoer

FInally the folder RFun and demographics contains R scripts and SPM batcher to perform correlation analysis at the ROI and voxel level respectively 
_ 
Dependencies
These scripts require the following dependencies:
-Raincloud plots
-Colorbrewer
-Fieldtrip
-SPM 12
-Rstudio
-Python 3 & Pysurfer
-Freesurfer
