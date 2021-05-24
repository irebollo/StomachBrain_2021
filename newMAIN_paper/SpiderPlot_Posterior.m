% Effect Sizes in other brain regions

%% First get effect sizes of all glasser regions

subjectsPhysiens = [8 9 10 13 15 16 17 18 19 20 21 22 23 24 26 29 31 32 33 34 35 36 37 38 39 40 41 43 44];
subjectsNavigastric = [1:10 12 13 15:19 21 22 23 25 26 27  29 30 31 32 33 34 35 36 38 91 96];


load('/media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/files/ROIs2Networks')


dataAll_phys_R = zeros(length(subjectsPhysiens),180)
dataAll_phys_L = zeros(length(subjectsPhysiens),180)


filename = '/media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/files/lh.HCPMMP1.annot'
[verticesHCP, labelHCP, colortableHC_LH] = read_annotation(filename);
filename = '/media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/files/rh.HCPMMP1.annot'
[verticesHCP, labelHCP, colortableHC_RH] = read_annotation(filename);


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

data_All_RH = [dataAll_phys_R;dataAll_navi_R]
data_All_LH = [dataAll_phys_L;dataAll_navi_L]

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






%% Then Get Cohen D data in ROI: 

% LOS_all = {'9-46d','a9-46v','a10p','10pp','10r','IFJa','s6-8','TE1m','TE1p','STSvp'}
% LOS_all = {'a9-46v','p10p','9-46d','a10p','10pp','10d','10r','10v','9m','6r','IFJa','IFSp','s6-8','8BL','TE1m','TE1p','STSvp'}
LOS_all = {'V1','V2','V3','V3A','V3B','VMV1','V4','V6','V6A','V7','VMV2','V8','VVC','VMV3','ProS','DVT','POS1','POS2','FFC','PIT','LO2','MT','PGp','TPOJ3','LO3','IPS1','RSC','v23ab','d23ab','PCV','7Am','7Pm','7PL','7Am','7m','31a','31pd','31pv'}

ROIS_L=zeros(1, length(LOS_all))
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

%% Then get overlap data in ROI


[vol, M, mr_parms, volsz]=load_mgh('/media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/lh.gasnetSURF.mgh');
filename = '//media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/files/lh.HCPMMP1.annot';
[verticesHCP, labelHCP, colortableHC_LH] = read_annotation(filename);

vol_t = vol>0.01;

nROI = length(unique(labelHCP));
% allROI = unique(labelHCP);

allROI = colortableHC_LH.table(:,5);

ROI_vertex = cell(nROI,1);
for iROI=1:nROI
    ROI_vertex{iROI}=find(labelHCP ==allROI(iROI));
end

OverlapLH = zeros(nROI,1)
for iROI=2:nROI
    % find index of all vertex of that ROI
    OverlapLH(iROI)=(sum(vol_t(ROI_vertex{iROI}))/length(ROI_vertex{iROI}))*100
    
end


[vol, M, mr_parms, volsz]=load_mgh('/media/irebollo/storage/projects/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/rh.gasnetSURF.mgh');
filename = '//media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/files/rh.HCPMMP1.annot';
[verticesHCP, labelHCP, colortableHC_RH] = read_annotation(filename);

vol_t = vol>0.01;

nROI = length(unique(labelHCP));
% allROI = unique(labelHCP);

allROI = colortableHC_RH.table(:,5);

ROI_vertex = cell(nROI,1);
for iROI=1:nROI
    ROI_vertex{iROI}=find(labelHCP ==allROI(iROI));
end

OverlapRH = zeros(nROI,1)
for iROI=2:nROI
    % find index of all vertex of that ROI
    OverlapRH(iROI)=(sum(vol_t(ROI_vertex{iROI}))/length(ROI_vertex{iROI}))*100
    
end


%% make polar plot witrh regions of interest % of overlap


categ = char(LOS_all);

dataOverlapR=OverlapRH(ROIS_R+1)' % this +1 is because of the first ROI is subcortical/ outside atlas
dataOverlapR(isnan(dataOverlapR))=0
dataOverlapL= OverlapLH(ROIS_L+1)'
dataOverlapL(isnan(dataOverlapL))=0

dim=size(dataOverlapR,2);          %set dimension or 'verteces'
data_num=size(dataOverlapR,1);
dat=zeros(data_num,dim+1);
dat=dataOverlapR;
dat(:,dim+1)=dataOverlapR(:,1);    %connect last point to 1st
data_num=size(dataOverlapR,1);
dat2=zeros(data_num,dim+1);
dat2=dataOverlapL;
dat2(:,dim+1)=dataOverlapL(:,1);    %connect last point to 1st
R=max(max(dataOverlapR));


t = (2*pi/dim);theta = 0:t:2*pi;k=1;
figure
p1=polarplot(theta,dat(:),'b','Linewidth',2);
hold on
p2=polarplot(theta,dat2(:),'r','Linewidth',2);
hold on


ax=gca;
ax.GridAlpha=1;
ax.GridColor='k';
ax.MinorGridColor='k';
ax.ThetaColor='k';
ax.RColor='k';
ax.RAxis.Color='k';
ax.RAxisLocation=90;
ax.RLim=[0 100];
ax.ThetaZeroLocation='top';
ax.ThetaDir='clockwise';
ax.ThetaTick=[0:360/dim:360];
ax.ThetaTickLabel={categ};
ax.FontSize=12;
ax.FontWeight='bold';

legend('Right','Left')


print ('-dpdf', '-painters', '/media/irebollo/storage/projects/Navigastric/scripts/plots&figures/OverlapSpiderPosterior')

%% Effect sizes


categ = char(LOS_all);


dataESR=cd_all_R
dataESR(isnan(dataESR))=0
dataESR(dataOverlapR<3)=0
dataESL= cd_all_L
dataESL(isnan(dataESL))=0
dataESL(dataOverlapL<3)=0

dim=size(dataESR,2);          %set dimension or 'verteces'
data_num=size(dataESR,1);
dat=zeros(data_num,dim+1);
dat=dataESR;
dat(:,dim+1)=dataESR(:,1);    %connect last point to 1st
data_num=size(dataESR,1);
dat2=zeros(data_num,dim+1);
dat2=dataESL;
dat2(:,dim+1)=dataESL(:,1);    %connect last point to 1st
R=max(max(dataESR));


t = (2*pi/dim);theta = 0:t:2*pi;k=1;
figure
p1=polarplot(theta,dat(:),'b','Linewidth',2);
hold on
p2=polarplot(theta,dat2(:),'r','Linewidth',2);
hold on


ax=gca;
ax.GridAlpha=1;
ax.GridColor='k';
ax.MinorGridColor='k';
ax.ThetaColor='k';
ax.RColor='k';
ax.RAxis.Color='k';
ax.RAxisLocation=90;
ax.RLim=[0 0.75];
ax.ThetaZeroLocation='top';
ax.ThetaDir='clockwise';
ax.ThetaTick=[0:360/dim:360];
ax.ThetaTickLabel={categ};
ax.FontSize=12;
ax.FontWeight='bold';

legend('Right','Left')


print ('-dpdf', '-painters', '/media/irebollo/storage/projects/Navigastric/scripts/plots&figures/ESSpiderPosterior')
