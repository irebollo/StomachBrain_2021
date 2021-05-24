% MakeTable3

subjectsPhysiens = [8 9 10 13 15 16 17 18 19 20 21 22 23 24 26 29 31 32 33 34 35 36 37 38 39 40 41 43 44];
subjectsNavigastric = [1:10 12 13 15:19 21 22 23 25 26 27  29 30 31 32 33 34 35 36 38 91 96];


load('/media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/files/ROIs2Networks')


dataAll_phys_R = zeros(length(subjectsPhysiens),180)
dataAll_phys_L = zeros(length(subjectsPhysiens),180)


filename = '//media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/files/lh.HCPMMP1.annot'
[verticesHCP, labelHCP, colortableHC_LH] = read_annotation(filename);
filename = '//media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/files/rh.HCPMMP1.annot'
[verticesHCP, labelHCP, colortableHC_RH] = read_annotation(filename);

[vertex_coordsGasnetL, facesGasnetL] = read_surf('/media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/lh.gasnet.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii')



load('/media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/files/ROIs2Networks7net')

%%


%% Order ROI according to RSN

filename = '/media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/files/rh.Yeo2011_7Networks_N1000.annot'
[verticesYeo, labelYeo, colortableYeo] = read_annotation(filename);

filename = '/media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/files/rh.HCPMMP1.annot'
[verticesHCP, labelHCP, colortableHC_LH] = read_annotation(filename);


% for each HCP parcellation ROI we want to know which RSNs it belongs to
% 1 - make a list of all the vexels belonging to each HCPlabel
% 2 - Iterate through each ROI and obtain for each vertex in each ROI which
% Yeo rsn it belongs too
% 3 - Take for each ROI the RSN with more VXL to it and assign it to
% ROI_rsn_label

nROI = length(unique(labelHCP));
% allROI = unique(labelHCP);

allROI = colortableHC_LH.table(:,5)

ROI_vertex = cell(nROI,1);
for iROI=1:nROI
ROI_vertex{iROI}=find(labelHCP ==allROI(iROI));
end


ROI_rsn = zeros(nROI,1);
for iROI=1:nROI
RSN_vertex_iROI = labelYeo(ROI_vertex{iROI});
ROI_rsn(iROI) = mode(RSN_vertex_iROI);  
end

unique(ROI_rsn)

ROI_rsn_label_lh = cell(1,nROI);


% Find the corresponding name of each StrucIndex
labelsRSN = {'subcortical','Visual','SomatoMotor','DA','Saliency','Limbic','Control','DMN'}
for iROI=2:nROI
index_label = colortableYeo.table(:,5) == ROI_rsn(iROI);
ROI_rsn_label_lh{iROI}=labelsRSN{index_label};
end
ROI_rsn_label_lh{1}='NA'
ROI_rsn_label_lh = ROI_rsn_label_lh'
% Sort each ROI according to which RSN it belongs too

sorted_ROIs_lh = []
index_sorted_ROIs_lh = []
for iNET=1:8
    indexTmp = find(colortableYeo.table(iNET,5) == ROI_rsn)
    sorted_ROIs_lh = [sorted_ROIs_lh;indexTmp]
    index_sorted_ROIs_temp = ones(length(indexTmp),1)*iNET;
    index_sorted_ROIs_lh = [index_sorted_ROIs_lh;index_sorted_ROIs_temp];
end

labels_HCP_RH = colortableHC_RH.struct_names(2:181);

%% Compute Overlaps

[vol, M, mr_parms, volsz]=load_mgh('/media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/lh.gasnetSURF.mgh')
filename = '//media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/files/lh.HCPMMP1.annot'
[verticesHCP, labelHCP, colortableHC_LH] = read_annotation(filename);

vol_t = vol>0.01;

nROI = length(unique(labelHCP));
% allROI = unique(labelHCP);

allROI = colortableHC_LH.table(:,5)

ROI_vertex = cell(nROI,1);
for iROI=1:nROI
ROI_vertex{iROI}=find(labelHCP ==allROI(iROI));
end

OverlapLH = zeros(nROI,1)
for iROI=2:nROI
    % find index of all vertex of that ROI
  OverlapLH(iROI)=(sum(vol_t(ROI_vertex{iROI}))/length(ROI_vertex{iROI}))*100
  
end


[vol, M, mr_parms, volsz]=load_mgh('/media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/rh.gasnetSURF.mgh')
filename = '//media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/files/rh.HCPMMP1.annot'
[verticesHCP, labelHCP, colortableHC_RH] = read_annotation(filename);

vol_t = vol>0.01;

nROI = length(unique(labelHCP));
% allROI = unique(labelHCP);

allROI = colortableHC_RH.table(:,5)

ROI_vertex = cell(nROI,1);
for iROI=1:nROI
ROI_vertex{iROI}=find(labelHCP ==allROI(iROI));
end

OverlapRH = zeros(nROI,1)
for iROI=2:nROI
    % find index of all vertex of that ROI
  OverlapRH(iROI)=(sum(vol_t(ROI_vertex{iROI}))/length(ROI_vertex{iROI}))*100
  
end



