#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 09:24:20 2019

@author: irebollo
"""

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
import os

"""
Define the three important variables.
Note that these are the first three positional arguments
in tksurfer (and pysurfer for that matter).
"""


subject_id = 'fsaverage'
surf = 'inflated'
subjects_dir="/media/irebollo/storage/projects/Navigastric/freesurferData/"

LOS_all = ['V1','V2','V3','V3A','V3B','VMV1','V4','V6','V6A','V7','VMV2','V8','VVC','VMV3','ProS','DVT','POS1','POS2','VIP','PCV','7PL','7Am','7Pm','7m','31a','31pd','31pv','PIT','LO2','PH','FFC','TF','MT','MST','PGp','TPOJ3','LO3','IP0','IPS1','RSC','v23ab','d23ab']

#

"""
Plot Left
"""



overlay_gasnetDL = "/media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/lh.gasnet.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii"
gasnet_l = io.read_scalar_data(overlay_gasnetDL)
gasnet_l[gasnet_l > 0.001] = 1


gasnetES_L = "/media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/lh.cd.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii"
sig1_l = io.read_scalar_data(gasnetES_L)
sig1_l[gasnet_l < 1 ] = 0

hemi="lh"


    

brain.add_data(sig1_l, min=0.3, max=0.5,thresh=0.001,colormap="autumn", alpha=1, colorbar=True,hemi="lh")




brain = Brain('fsaverage', 'split', 'inflated', views=['lat', 'med'],background='white')

brain.add_data(sig1_l, min=0.3, max=0.5,thresh=0.001,colormap="autumn", alpha=1, colorbar=True,hemi="lh")

brain.add_data(sig1_r, min=0.01, max=0.15,colormap="plasma", alpha=1, colorbar=True,hemi='rh')






##### Right hemisphere

overlay_gasnetDR = "/media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/rh.gasnet.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii"
gasnet_r = io.read_scalar_data(overlay_gasnetDR)
gasnet_r[gasnet_r > 0.001] = 1


gasnetES_R = "/media/irebollo/storage/projects/Navigastric/ClusterResults/effecsizes/rh.cd.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii"
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
for name in LOI_r: 
    data[labels[name].vertices] = count+1
    count= count+1

brain.add_data(data, min=1, max=count+1,thresh=1,colormap="RdBu", alpha=1, colorbar=False)


count=0
for name in LOS_all: 
    
    LOI_r = 'R_'+LOS_all[count]+'_ROI'
    label_name = LOI_r
    brain.add_label(label_name,borders=True,color="k")
    count = count +1
    
# plot
    
brain.add_data(sig1_r, min=0.3, max=0.5,thresh=0.001,colormap="autumn", alpha=1, colorbar=True)


brain.show_view(dict(azimuth=40, elevation=-100, distance=400))
brain.save_image("/media/irebollo/storage/projects/Navigastric/scripts/plots&figures/VIS_RM.png")

brain.show_view(dict(azimuth=115, elevation=-100, distance=400))

brain.save_image("/media/irebollo/storage/projects/Navigastric/scripts/plots&figures/VIS_RL.png")


brain.toggle_toolbars(True)

