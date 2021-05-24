%{
This script will 

1- Load the regions from glasser et al 2017 parcellation
2- Load individual coupling strenght maps
3- recompute effect sizes
4- Calculare mean effect size per glasser roi
4- Plot


Ignacio Rebollo
04-05-2021

%}

%% Settings

subjectsPhysiens = [8 9 10 13 15 16 17 18 19 20 21 22 23 24 26 29 31 32 33 34 35 36 37 38 39 40 41 43 44];
subjectsNavigastric = [1:10 12 13 15:19 21 22 23 25 26 27  29 30 31 32 33 34 35 36 38 91 96];


% file with labels of HCP-atlas ROI
load('/mnt/RAIDSTORAGE/SADREST/scripts_matlab/files/ROIs2Networks')

% Load HCP atlas
filename = '/mnt/RAIDSTORAGE/SADREST/files/lh.HCPMMP1.annot'
[verticesHCP, labelHCP, colortableHC_LH] = read_annotation(filename);
filename = '/mnt/RAIDSTORAGE/SADREST/files/rh.HCPMMP1.annot'
[verticesHCP, labelHCP, colortableHC_RH] = read_annotation(filename);


% Preallocate ROIs
dataAll_phys_R = zeros(length(subjectsPhysiens),180)
dataAll_phys_L = zeros(length(subjectsPhysiens),180)



% Store coupling strength persubject
% Sample 2
for iSubj = 1:length(subjectsPhysiens)
    subj_idx=subjectsPhysiens(iSubj)
    subj_path= ['/media/irebollo/storage/projects/Navigastric/dataPhysiens/subj',sprintf('%.2d',subj_idx),'/Timeseries/fsparcels/']
%     data_lh = load ([subj_path '/cs_hcp_par_RFANTS.lh.dat']);
%     data_rh = load ([subj_path '/cs_hcp_par_RFANTS.rh.dat']);
    data_lh = load ([subj_path '/cs_hcp_par_RFANTS_nt.lh.dat']);
    data_rh = load ([subj_path '/cs_hcp_par_RFANTS_nt.rh.dat']);
dataAll_phys_R(iSubj,:) = [data_rh];
dataAll_phys_L(iSubj,:) = [data_lh];

end

% sample 1
dataAll_navi_R = zeros(length(subjectsNavigastric),180)
dataAll_navi_L = zeros(length(subjectsNavigastric),180)


for iSubj = 1:length(subjectsNavigastric)
    subj_idx=subjectsNavigastric(iSubj)
    subj_path= ['/media/irebollo/storage/projects/Navigastric/subjects/subj',sprintf('%.2d',subj_idx),'/Timeseries/fsparcels/']
%     data_lh = load ([subj_path '/cs_hcp_par_RFANTS.lh.dat']);
%     data_rh = load ([subj_path '/cs_hcp_par_RFANTS.rh.dat']);
        data_lh = load ([subj_path '/cs_hcp_par_RFANTS_nt.lh.dat']);
    data_rh = load ([subj_path '/cs_hcp_par_RFANTS_nt.rh.dat']);
dataAll_navi_R(iSubj,:) = [data_rh];
dataAll_navi_L(iSubj,:) = [data_lh];

end


% concatenate data from both samples
data_All_RH = [dataAll_phys_R;dataAll_navi_R]
data_All_LH = [dataAll_phys_L;dataAll_navi_L]


% compute effect size from coupling strength
[h p ci statsAll_L] = ttest(data_All_LH)
[h p ci statsAll_R] = ttest(data_All_RH)

cd_allROIS_R = statsAll_R.tstat/sqrt(length(subjectsNavigastric)+length(subjectsPhysiens))
cd_allROIS_L = statsAll_L.tstat/sqrt(length(subjectsNavigastric)+length(subjectsPhysiens))



CDbothBothSides = [cd_allROIS_R;cd_allROIS_L];
nansEitherSide = isnan(CDbothBothSides);
nansEitherSide = sum(nansEitherSide)>0;
CDbothBothSides_nonanEitherside = CDbothBothSides;
CDbothBothSides_nonanEitherside(:,nansEitherSide)=[];


CDbothAVG = mean([CDbothBothSides_nonanEitherside])


[B,Idx] =sort(CDbothAVG,'descend')





namesR = colortableHC_RH.struct_names(2:181);
namesR(nansEitherSide)=[]
namesR(:) = namesR(Idx);

for in=1:length(namesR)
    currentName = namesR{in}
    namesRegions{in}=currentName(3:length(currentName)-4)
end


%% Figure 1G 25 regions with largest effect sizes

N_regions = 25;
LOS_all = namesRegions(1:N_regions)

AllES = figure

ROIS_L = zeros (1, length(LOS_all))
ROIS_R=zeros(1, length(LOS_all))
for i =1:length(LOS_all) 
    ROIS_L(i) = find(strcmp(colortableHC_LH.struct_names,strcat('L_',LOS_all{i},'_ROI')))-1
        ROIS_R(i) = find(strcmp(colortableHC_RH.struct_names,strcat('R_',LOS_all{i},'_ROI')))-1
end

[h p ci statsAll_L] = ttest(data_All_LH(:,ROIS_L))
t_roi_all_L = statsAll_L.tstat
[h p ci statsAll_R] = ttest(data_All_RH(:,ROIS_R))
t_roi_all_R = statsAll_R.tstat


cd_all_R = t_roi_all_R/sqrt(length(subjectsNavigastric)+length(subjectsPhysiens))
cd_all_L = t_roi_all_L/sqrt(length(subjectsNavigastric)+length(subjectsPhysiens))

cd_all_both = nanmean([cd_all_R;cd_all_L])'

bar(cd_all_both(1:N_regions),'FaceColor','b')
hold on;
% if we want each bar to have a different color, loop
% for b = 1:length(dat)
%     bar(b*1.2, bR(b), 'FaceColor',  colors(idxR(b), : ), 'EdgeColor', 'k', 'BarWidth', 1);
% end


labelsPlot =namesRegions(1:N_regions) 
set(gca,'XTick',[1:N_regions])
set(gca,'XTickLabel',LOS_all)
grid on
title('Top 25 regions with largest effect sizes')
set(gca,'FontSize',14)
set(gca,'XTickLabelRotation',45)

hold on


txtBg = 'none';
ax = gca;
 fz = 10; fontweight = 'bold';
 txt='>'
 xpos = [1:1:N_regions]
%  valuesRight = cd_allROIS_R_nn(Idx(1:N_regions))
 color =  [0.5 0.5 0.5]
for iP = 1:N_regions


    if ~isnan(cd_all_R(iP))
    % draw the stars in the bar
    h = text(xpos(iP), cd_all_R(iP), txt, ...
        'horizontalalignment', 'center', 'backgroundcolor', ...
        txtBg, 'margin', 1, 'fontsize', fz, 'fontweight', fontweight, 'color', color, 'Parent', ax);
    else
    end
end


txtBg = 'none';
ax = gca;
 fz = 10; fontweight = 'bold';
 txt='<'
 xpos = [1:1:N_regions]
%  valuesRight = cd_allROIS_R_nn(Idx(1:N_regions))
 color =  [0.5 0.5 0.5]
for iP = 1:N_regions


    if ~isnan(cd_all_L(iP))
    % draw the stars in the bar
    h = text(xpos(iP), cd_all_L(iP), txt, ...
        'horizontalalignment', 'center', 'backgroundcolor', ...
        txtBg, 'margin', 1, 'fontsize', fz, 'fontweight', fontweight, 'color', color, 'Parent', ax);
    else
    end
end

ylim([0.4 0.85])
ylabel('Cohen D')
figure(AllES)

h= gcf
set(h,'PaperOrientation','landscape');
set(h,'PaperPositionMode','auto');  

    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.3])
    set(gca,'FontSize',14)
grid on
    print ('-dpdf', '-painters', '/media/irebollo/storage/projects/Navigastric/scripts/plots&figures/EeffectSizeTOP25')

  
    labelsPlot= labelsPlot'