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
hemi = 'lh'


subjects_dir="/media/irebollo/storage/projects/Navigastric/freesurferData/"


"""
LEFT
"""

overlay_gasnetDL = "/media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/lh.gasnet.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii"
gasnetES_L = "/media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/lh.cd.FUSPROJ.nii"




labels = mne.read_labels_from_annot('/media/irebollo/storage/projects/Navigastric/freesurferData/fsaverage', parc='Yeo2011_7Networks_N1000', hemi=hemi)
labels = {l.name: l for l in labels}

              
#brain.add_annotation("aparc.a2009s")


# overlay_gasnetDL = "/media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/lh.gasnet.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii"
gasnet_l = io.read_scalar_data(overlay_gasnetDL)
gasnet_l[gasnet_l > 0.001] = 1


# gasnetES_L = "/media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/lh.cd.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii"
sig1_l = io.read_scalar_data(gasnetES_L)
sig1_l[gasnet_l < 1 ] = 0



hemi="lh"



labels = mne.read_labels_from_annot('/media/irebollo/storage/projects/Navigastric/freesurferData/fsaverage', parc='Yeo2011_7Networks_N1000', hemi=hemi)
labels = {l.name: l for l in labels}



brain = Brain(subject_id,hemi,surf,subjects_dir="/media/irebollo/storage/projects/Navigastric/freesurferData/",background="w")
V = brain.geo[hemi].x.size
data = np.zeros((V))


LOS_all = ['7Networks_1','7Networks_2','7Networks_3','7Networks_4','7Networks_5','7Networks_6','7Networks_7']
LOI_r = [None] *len(LOS_all)
LOI_l = [None] *len(LOS_all)
count=0
for name in LOS_all: 
    LOI_r[count] = 'R_'+LOS_all[count]+'_ROI-rh'
    LOI_l[count] = LOS_all[count]+'-lh'
    count= count+1

count=0
for name in LOI_l: 
    count= count+1
    data[labels[name].vertices] = count
#brain.add_data(data, min=1, max=count,thresh=0.1,colormap="RdBu",alpha=1, colorbar=False)

LOS_all_colors = ['fuchsia','cyan','lime','yellow','lightyellow','blue','red']

count=0
for name in LOS_all: 
    #LOI_l = 'L_'+LOS_all[count]+'_ROI'
    LOI_l = LOS_all[count]
    label_name = LOI_l
    brain.add_label(label_name,borders=False,color=LOS_all_colors[count],alpha=0.25)
    count = count +1

count=0
for name in LOS_all: 
    #LOI_l = 'L_'+LOS_all[count]+'_ROI'
    LOI_l = LOS_all[count]
    label_name = LOI_l
    brain.add_label(label_name,borders=1,color='w',alpha=0.95)
    count = count +1
    

brain.add_data(sig1_l, min=0.3, max=0.5,thresh=0.001,colormap="autumn", alpha=1, colorbar=True)
brain.show_view(dict(azimuth=0, elevation=-90, distance=500))
brain.save_image("/media/irebollo/storage/projects/Navigastric/scripts/plots&figures/RSNCD_LL.png")
brain.show_view(dict(azimuth=0, elevation=90, distance=500))
brain.save_image("/media/irebollo/storage/projects/Navigastric/scripts/plots&figures/RSNCD_LM.png")

brain.toggle_toolbars(True)
#Needed to change VTK Data - /7networks_1, Actor, Mapper more options, resolve coincindent topology to default
##### Right hemisphere


overlay_gasnetDR = "/media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/rh.gasnet.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii"
gasnetES_R = "/media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/rh.cd.FUSPROJ.nii"

# overlay_gasnetDR = "/media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/rh.gasnet.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii"
gasnet_r = io.read_scalar_data(overlay_gasnetDR)
gasnet_r[gasnet_r > 0.001] = 1


# gasnetES_R = "/media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/rh.cd.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii"
sig1_r = io.read_scalar_data(gasnetES_R)
sig1_r[gasnet_r < 1 ] = 0


hemi="rh"


brain = Brain(subject_id,hemi,surf,subjects_dir="/media/irebollo/storage/projects/Navigastric/freesurferData/",background="w")
V = brain.geo[hemi].x.size
data = np.zeros((V))


labels = mne.read_labels_from_annot('/media/irebollo/storage/projects/Navigastric/freesurferData/fsaverage', parc='Yeo2011_7Networks_N1000', hemi=hemi)
labels = {l.name: l for l in labels}

LOS_all = ['7Networks_1','7Networks_2','7Networks_3','7Networks_4','7Networks_5','7Networks_6','7Networks_7']
LOI_r = [None] *len(LOS_all)
LOI_l = [None] *len(LOS_all)

count=0
for name in LOS_all: 
    LOI_r[count] = LOS_all[count]+'-rh'
    LOI_l[count] = LOS_all[count]+'-lh'
    count= count+1

count=0
for name in LOI_r: 
    count= count+1
    data[labels[name].vertices] = count
    
    
#brain.add_data(data, min=1, max=count,thresh=0.1,colormap="RdBu",alpha=1, colorbar=False)

LOS_all_colors = ['fuchsia','cyan','lime','yellow','lightyellow','blue','red']

count=0
for name in LOS_all: 
    #LOI_l = 'L_'+LOS_all[count]+'_ROI'
    LOI_l = LOS_all[count]
    label_name = LOI_l
    brain.add_label(label_name,borders=False,color=LOS_all_colors[count],alpha=0.15)
    count = count +1

count=0
for name in LOS_all: 
    #LOI_l = 'L_'+LOS_all[count]+'_ROI'
    LOI_l = LOS_all[count]
    label_name = LOI_l
    brain.add_label(label_name,borders=1,color='w',alpha=0.7)
    count = count +1

    
brain.add_data(sig1_r, min=0.3, max=0.5,thresh=0.001,colormap="autumn", alpha=1, colorbar=True)

brain.show_view(dict(azimuth=0, elevation=90, distance=500))
brain.save_image("/media/irebollo/storage/projects/Navigastric/scripts/plots&figures/RSNCD_RL.png")
brain.show_view(dict(azimuth=0, elevation=-90, distance=500))
brain.save_image("/media/irebollo/storage/projects/Navigastric/scripts/plots&figures/RSNCD_RM.png")





