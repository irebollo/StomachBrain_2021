%-----------------------------------------------------------------------
% Job saved on 05-May-2021 19:05:57 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.util.imcalc.input = {
                                        '/media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/otherParcellLocalizers/gradientsMargulies/volume.all.1.nii,1'
                                        '/media/irebollo/storage/projects/Navigastric/ClusterResults/All_63SA0250_final/SurrParBothSample_nR1000_CA0050_kw3_fir2_fspread_015_fOrder_5_tw_15csfr_mask.nii,1'
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'gasnet_marguliesSpace.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {'/media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/otherParcellLocalizers/gradientsMargulies'};
matlabbatch{1}.spm.util.imcalc.expression = 'i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
