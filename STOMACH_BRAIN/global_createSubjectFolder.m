function global_createSubjectFolder(subj_idx,sample)

rootDir= strcat(global_path2root(sample));
% 

mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Brainamp',filesep,'Tasks'))
% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Brainamp',filesep,'Tasks',filesep,'EBA'))
mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Brainamp',filesep,'Tasks',filesep,'MOTOR'))
mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Brainamp',filesep,'Tasks',filesep,'SOMAFLOW'))
% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Brainamp',filesep,'Tasks',filesep,'MT'))


% rootDir= strcat('D:\NAVIGASTRIC\test2pipelines',filesep);
% 
% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx)))
% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Brainamp'))
% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'fMRI'))
% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'taskfMRI'))
% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'taskfMRI',filesep,'EBA'))
% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'taskfMRI',filesep,'SOMAFLOW'))
% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'taskfMRI',filesep,'MOTOR'))
% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'taskfMRI',filesep,'MT'))

mkdir(strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),filesep,'Regress'))
mkdir(strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),filesep,'stats'))

% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'anatomy'))
% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries'))
% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'PreprocessingLog'))
% mkdir(strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'EGGTimeseries'))

mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'EGGTimeseries',filesep,'EBA'))
mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'EGGTimeseries',filesep,'MOTOR'))
mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'EGGTimeseries',filesep,'SOMAFLOW'))
% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'EGGTimeseries',filesep,'MT'))
mkdir(strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'EGGTimeseries',filesep,'REST'))


% mkdir(strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'MRItimeseries'))
% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'MRItimeseries',filesep,'EBA'))
mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'MRItimeseries',filesep,'MOTOR'))
mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'MRItimeseries',filesep,'SOMAFLOW'))
mkdir(strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'MRItimeseries',filesep,'REST'))
% 
% 
mkdir(strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'PhasesAnalysis'))
% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'PhasesAnalysis',filesep,'EBA'))
mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'PhasesAnalysis',filesep,'MOTOR'))
mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'PhasesAnalysis',filesep,'SOMAFLOW'))
mkdir(strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'PhasesAnalysis',filesep,'REST'))

% navigastric
mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'ControlEGG_othersubjectNavigastric'))
% mkdir(strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'Heart'))


%physiens
% mkdir(strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'ControlEGG_othersubjectBOTHSAMPLES'))
% mkdir(strcat(rootDir,'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'GlobalSignal'))





