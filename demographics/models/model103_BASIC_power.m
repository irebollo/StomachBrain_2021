%-----------------------------------------------------------------------
% Job saved on 09-Oct-2018 16:38:17 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6906)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------


modelName = 'model103_BASIC_Power'
models_dir = [global_path2root_folder(),'\ClusterResults\demographics\models\']
if ~exist([models_dir modelName])

   mkdir([models_dir modelName]) 
end
matlabbatch{1}.spm.stats.factorial_design.dir = {[models_dir modelName]};
%% input images
        
subjectsPhysiens = global_subjectList_phys;
subjectsNavigastric = global_subjectList_navi;

subjects_all = [subjectsNavigastric,subjectsPhysiens];


root4csimages=[global_path2root_folder,'ClusterResults',filesep,'demographics',filesep,'cs_4spmGLM'];

powerDataOutlier = [0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	1	0	0	0	0]

subjectImages_temp=[]
subjectCounter = 0
for iS=1:length(powerDataOutlier)
    
    if powerDataOutlier(iS)==0
        subjectCounter = subjectCounter+1
    
    if iS<=length(subjectsNavigastric)
        subj_idx = subjectsNavigastric(iS)
    else
        subj_idx = subjectsPhysiens(iS-length(subjectsNavigastric))
    end
    filename = [root4csimages,filesep,'s6_csi_',sprintf('%.2d',iS),'_subj_idx_',num2str(subj_idx),'.nii']
    subjectImages_temp{subjectCounter,1} = [filename,',1']
    
    end
end

matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans=subjectImages_temp


%%
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov = struct('c', {}, 'cname', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov.files = {[global_path2root_folder(),'\scripts\demographics\covariates\model103_covariates.txt']};
matlabbatch{1}.spm.stats.factorial_design.multi_cov.iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.multi_cov.iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {[global_path2root_folder(),'\scripts\STOMACH_BRAIN\files\SPM_mask_THRESHOLDED_physiens.img,1']};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat = {[models_dir modelName '\SPM.mat']};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{3}.spm.stats.con.spmmat = {[models_dir modelName '\SPM.mat']};
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'intercept';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Female';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'malee+';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 -1];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = ' BMI+';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 0 1];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'BMI-';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 0 -1];
matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'FWD+';
matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 1];
matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'FWD -';
matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 -1];
matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'EGGfreq +';
matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [0 0 0 0 1];
matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'EGGfreq -';
matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [0 0 0 0 -1];
matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = 'TOD+';
matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = [0 0 0 0 0 1];
matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{11}.tcon.name = 'TOD-';
matlabbatch{3}.spm.stats.con.consess{11}.tcon.weights = [0 0 0 0 0 -1];
matlabbatch{3}.spm.stats.con.consess{11}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{12}.tcon.name = 'Group+';
matlabbatch{3}.spm.stats.con.consess{12}.tcon.weights = [0 0 0 0 0 0 1];
matlabbatch{3}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{13}.tcon.name = 'Group-';
matlabbatch{3}.spm.stats.con.consess{13}.tcon.weights = [0 0 0 0 0 0 -1];
matlabbatch{3}.spm.stats.con.consess{13}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{14}.tcon.name = 'Power+';
matlabbatch{3}.spm.stats.con.consess{14}.tcon.weights = [0 0 0 0 0 0 0 1];
matlabbatch{3}.spm.stats.con.consess{14}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{15}.tcon.name = 'Power-';
matlabbatch{3}.spm.stats.con.consess{15}.tcon.weights = [0 0 0 0 0 0 0 -1];
matlabbatch{3}.spm.stats.con.consess{15}.tcon.sessrep = 'none';



matlabbatch{3}.spm.stats.con.delete = 1;
matlabbatch{4}.spm.stats.results.spmmat = {[models_dir modelName '\SPM.mat']};
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.png = true;


%% run job

nrun = 1; % enter the number of runs here
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch, inputs{:});
