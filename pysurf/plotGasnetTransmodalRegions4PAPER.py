"""
Basic Visualization
===================

Initialize a basic visualization session.

"""

from surfer import Brain
from mayavi import mlab
import numpy as np
from surfer import io
import mne

"""
Define the three important variables.
Note that these are the first three positional arguments
in tksurfer (and pysurfer for that matter).
"""



subject_id = 'fsaverage'
surf = 'inflated'
subjects_dir="/media/irebollo/storage/projects/Navigastric/freesurferData/"

# Define ROIs and colors of interest

LOS_all = ['a9-46v','p10p','9-46d','a10p','10pp','11l','s6-8','8BL',
           '10d','10r','10v',
           '6r','IFJa','IFSp','p9-46v','44',
           'TE1m','TE1p','TE1a','STSvp','STSva']

colors_all = ['green','green','green','green','green','green','green','green',
              'goldenrod','goldenrod','goldenrod',
              'blue','blue','blue','blue','blue',
              'pink','pink','pink','pink','pink']


# %% PLot labels LH

overlay_gasnetDL = "/media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/lh.gasnet.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii"
gasnetES_L = "/media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/lh.cd.FUSPROJ.nii"


# overlay_gasnetDL = "/media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/lh.gasnet.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii"
gasnet_l = io.read_scalar_data(overlay_gasnetDL)
gasnet_l[gasnet_l > 0.001] = 1


# gasnetES_L = "/media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/lh.cd.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii"
sig1_l = io.read_scalar_data(gasnetES_L)
sig1_l[gasnet_l < 1 ] = 0

hemi="lh"
labels = mne.read_labels_from_annot('/media/irebollo/storage/projects/Navigastric/freesurferData/fsaverage', parc='HCPMMP1', hemi=hemi)
labels = {l.name: l for l in labels}


LOI_l = [None] *len(LOS_all)
count=0
for name in LOS_all: 
    LOI_l[count] = 'L_'+LOS_all[count]+'_ROI-lh'
    count= count+1
    


brain = Brain(subject_id,hemi,surf,subjects_dir="/media/irebollo/storage/projects/Navigastric/freesurferData/",background="w")
V = brain.geo[hemi].x.size
data = np.zeros((V))

count=0
for name in LOS_all: 
    LOI_l = 'L_'+LOS_all[count]+'_ROI'
    label_name = LOI_l
    brain.add_label(label_name,borders=True,color=colors_all[count],alpha=0.85)
    count = count +1

'''
count=0
for name in LOS_all: 
    LOI_l = 'L_'+LOS_all[count]+'_ROI'
    label_name = LOI_l
    brain.add_label(label_name,borders=3,color=colors_all[count],alpha=0.85)
    count = count +1
  '''  

# Compute overlap with all parcellations

brain.add_data(sig1_l, min=0.3, max=0.5,thresh=0.001,colormap="autumn", alpha=1, colorbar=True)

brain.toggle_toolbars(True)


# %% Visualization


brain.show_view(dict(azimuth=180, elevation=100, distance=450))
brain.save_image("/media/irebollo/storage/projects/Navigastric/scripts/plots&figures/Other_LL.png")


brain.show_view(dict(azimuth=83, elevation=95, distance=400))
brain.save_image("/media/irebollo/storage/projects/Navigastric/scripts/plots&figures/Other_LF.png")


brain.show_view(dict(azimuth=85, elevation=50, distance=500))
brain.save_image("/media/irebollo/storage/projects/Navigastric/scripts/plots&figures/Other_LS.png")



# %% PLot labels RH

overlay_gasnetDR = "/media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/rh.gasnet.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii"
gasnetES_R = "/media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/rh.cd.FUSPROJ.nii"

# overlay_gasnetDR = "/media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/rh.gasnet.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii"
gasnet_r = io.read_scalar_data(overlay_gasnetDR)
gasnet_r[gasnet_r > 0.001] = 1


# gasnetES_R = "/media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/rh.cd.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii"
sig1_r = io.read_scalar_data(gasnetES_R)
sig1_r[gasnet_r < 1 ] = 0

hemi="rh"
labels = mne.read_labels_from_annot('/media/irebollo/storage/projects/Navigastric/freesurferData/fsaverage', parc='HCPMMP1', hemi=hemi)
labels = {l.name: l for l in labels}


LOI_r = [None] *len(LOS_all)
count=0
for name in LOS_all: 
    LOI_r[count] = 'R_'+LOS_all[count]+'_ROI-rh'
    count= count+1
    


brain = Brain(subject_id,hemi,surf,subjects_dir="/media/irebollo/storage/projects/Navigastric/freesurferData/",background="w")
V = brain.geo[hemi].x.size
data = np.zeros((V))

count=0
for name in LOS_all: 
    LOI_r = 'R_'+LOS_all[count]+'_ROI'
    label_name = LOI_r
    brain.add_label(label_name,borders=True,color=colors_all[count],alpha=0.85)
    count = count +1


# Compute overlap with all parcellations

brain.add_data(sig1_r, min=0.3, max=0.5,thresh=0.001,colormap="autumn", alpha=1, colorbar=True)

brain.toggle_toolbars(True)

# %% Visualize


brain.show_view(dict(azimuth=97, elevation=95, distance=400))
brain.save_image("/media/irebollo/storage/projects/Navigastric/scripts/plots&figures/Other_RF.png")


brain.show_view(dict(azimuth=10, elevation=100, distance=450))
brain.save_image("/media/irebollo/storage/projects/Navigastric/scripts/plots&figures/Other_RL.png")

brain.show_view(dict(azimuth=95, elevation=50, distance=500))
brain.save_image("/media/irebollo/storage/projects/Navigastric/scripts/plots&figures/Other_RS.png")


