%{
This script will produice coupling strenght nifti images
by substracting chance from empiricaö PLV volumes

Ignacio Rebollo
03-05-2021

%}

% store CS images all subjects

subjectsPhysiens = [8 9 10 13 15 16 17 18 19 20 21 22 23 24 26 29 31 32 33 34 35 36 37 38 39 40 41 43 44];
subjectsNavigastric = [1:10 12 13 15:19 21 22 23 25 26 27  29 30 31 32 33 34 35 36 38 91 96];
cfgMain.sample = 1
cfgMain = global_getcfgmain


for iS=1:length(subjectsNavigastric)
    subj_idx = subjectsNavigastric(iS);

    % empirical and chance PLV filenames
    
    filenamePLV = strcat(global_filename(subj_idx,cfgMain,strcat('PLVXVoxelFilename_',cfgMain.Timeseries2Regress)),'.nii');
       filenamePLVSurrogate  = strcat(global_path2root(1),'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'ControlEGG_othersubjectBOTHSAMPLES',filesep,'median_sPLV_s',sprintf('%.2d',subj_idx),'.nii');
    
       
      PLVGroupEmpirical= ft_read_mri(filenamePLV); % Put into cell

    
    % Load surrogate PLV and prepare structure for surrogate PLV
    
    PLVGroupSurrogate = ft_read_mri(filenamePLVSurrogate);

    
    empirical = PLVGroupEmpirical.anatomy(:);
    surrogate= PLVGroupSurrogate.anatomy(:);
    
cs = empirical-surrogate;

        couplingStrenghtFilename = global_filename(subj_idx,cfgMain,'couplingStrenghtFilename')

tools_writeMri(cs,couplingStrenghtFilename)

    
end

    cfgMain.sample = 2
for iS=1:length(subjectsPhysiens)
    
    subj_idx = subjectsPhysiens(iS);
    
    % empirical and chance PLV filenames
    
    filenamePLV = strcat(global_filename(subj_idx,cfgMain,strcat('PLVXVoxelFilename_',cfgMain.Timeseries2Regress)),'.nii');
       filenamePLVSurrogate  = strcat(global_path2root(2),'subj',sprintf('%.2d',subj_idx),filesep,'Timeseries',filesep,'ControlEGG_othersubjectBOTHSAMPLES',filesep,'median_sPLV_s',sprintf('%.2d',subj_idx),'.nii');
    
    
           PLVGroupEmpirical= ft_read_mri(filenamePLV); % Put into cell

    
    % Load surrogate PLV and prepare structure for surrogate PLV
    
    PLVGroupSurrogate = ft_read_mri(filenamePLVSurrogate);

    
    empirical = PLVGroupEmpirical.anatomy(:);
    surrogate= PLVGroupSurrogate.anatomy(:);
    
    cs = empirical-surrogate;

        couplingStrenghtFilename = global_filename(subj_idx,cfgMain,'couplingStrenghtFilename')

tools_writeMri(cs,couplingStrenghtFilename)



end


