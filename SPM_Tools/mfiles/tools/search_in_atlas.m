atlas = ft_read_atlas('C:\MatlabToolboxes\fieldtrip-20140901\template\atlas\aal\ROI_MNI_V4.nii')
label = atlas_lookup(atlas, [33 -22 70], 'queryrange',1 , 'inputcoord','mni')