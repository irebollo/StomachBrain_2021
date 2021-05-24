subjectsPhysiens = [8 9 10 13 15 16 17 18 19 20 21 22 23 24 26 29 31 32 33 34 35 36 37 38 39 40 41 43 44];
subjectsNavigastric = [1:10 12 13 15:19 21 22 23 25 26 27  29 30 31 32 33 34 35 36 38 91 96];


load('/media/irebollo/DATAPART1/SADREST/scripts_matlab/files/ROIs2Networks')


dataAll_phys_R = zeros(length(subjectsPhysiens),180)
dataAll_phys_L = zeros(length(subjectsPhysiens),180)


filename = '//media/irebollo/DATAPART1/SADREST/files/lh.HCPMMP1.annot'
[verticesHCP, labelHCP, colortableHC_LH] = read_annotation(filename);
filename = '//media/irebollo/DATAPART1/SADREST/files/rh.HCPMMP1.annot'
[verticesHCP, labelHCP, colortableHC_RH] = read_annotation(filename);


for iSubj = 1:length(subjectsPhysiens)
    subj_idx=subjectsPhysiens(iSubj)
    subj_path= ['/media/Navigastric/dataPhysiens/subj',sprintf('%.2d',subj_idx),'/Timeseries/fsparcels/']
%     data_lh = load ([subj_path '/cs_hcp_par_RFANTS.lh.dat']);
%     data_rh = load ([subj_path '/cs_hcp_par_RFANTS.rh.dat']);
    data_lh = load ([subj_path '/csALLBRAIN_hcp_par_RFANTS.lh.dat']);
    data_rh = load ([subj_path '/csALLBRAIN_hcp_par_RFANTS.rh.dat']);
dataAll_phys_R(iSubj,:) = [data_rh];
dataAll_phys_L(iSubj,:) = [data_lh];

end


dataAll_navi_R = zeros(length(subjectsNavigastric),180)
dataAll_navi_L = zeros(length(subjectsNavigastric),180)

for iSubj = 1:length(subjectsNavigastric)
    subj_idx=subjectsNavigastric(iSubj)
    subj_path= ['/media/Navigastric/subjects/subj',sprintf('%.2d',subj_idx),'/Timeseries/fsparcels/']
%     data_lh = load ([subj_path '/cs_hcp_par_RFANTS.lh.dat']);
%     data_rh = load ([subj_path '/cs_hcp_par_RFANTS.rh.dat']);
        data_lh = load ([subj_path '/csALLBRAIN_hcp_par_RFANTS.lh.dat']);
    data_rh = load ([subj_path '/csALLBRAIN_hcp_par_RFANTS.rh.dat']);
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







namesR = colortableHC_RH.struct_names(2:181);
namesR(nansEitherSide)=[]
namesR(:) = namesR(Idx);

for in=1:length(namesR)
    currentName = namesR{in}
    namesRegions{in}=currentName(3:length(currentName)-4)
end


%% Figure top 1/3 effect sizes
LOS_all = {'4','V1','7m','31pd','31pv','LO3','TPOJ3','PGp','MT','Ig','AVI','Pir'}


ROIS_L = zeros (1, length(LOS_all))
ROIS_R=zeros(1, length(LOS_all))
for i =1:length(LOS_all) 
    ROIS_L(i) = find(strcmp(colortableHC_LH.struct_names,strcat('L_',LOS_all{i},'_ROI')))-1
        ROIS_R(i) = find(strcmp(colortableHC_RH.struct_names,strcat('R_',LOS_all{i},'_ROI')))-1
end


N_regions = length(LOS_all)

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

ylim([0.0 0.8])
ylabel('Cohen D')
figure(AllES)


%% get effect sizes 

cs_LH=data_All_LH(:,ROIS_L)
cs_RH=data_All_RH(:,ROIS_R)
ttest()

%% bayes
nobs = 63; % number of observations
xref = 2; % significant t value for one-sided ttest 29 degrees of freedom reference effect size (significant effect)
% source: http://www.ttable.org/uploads/2/1/7/9/21795380/9754276.png?852

for iRoi =1:12
xdat = cs_LH(:,iRoi)'

[bf_log10(iRoi)]= my_ttest_bayes(xdat, xref);
[res_Bayes(iRoi,:), bf(iRoi)] = interpret_Bayes(bf_log10(iRoi))
BF(iRoi)=  1/(10^bf(iRoi))

end

%% Bayes against median gray cs

median_cs_gray = empPLV_medianGray - chancePLV_medianGray;
median_cs_gastric = empPLV_medianGastric - chancePLV_medianGastric;
difference_gray_gastric = (median_cs_gastric - median_cs_gray)/2

coupling_strenght_roi = empPLV_ROI - chancePLV_ROI


for iRoi=1:6
   [h p ci stats]= ttest2(coupling_strenght_roi(iRoi,:),difference_gray_gastric)
   pp(iRoi)= p  
   tt(iRoi) = stats.tstat
end

for iRoi=1:6
t1smpbf(tt(iRoi),30)
end

for iRoi=1:6
   mean(coupling_strenght_roi(iRoi,:)-difference_gray_gastric)
end




nobs = 30; % number of observations
xref = +1.699; % significant t value for one-sided ttest 29 degrees of freedom reference effect size (significant effect)
% xref = +2.042; % significant t value for two-sided ttest 29 degrees of freedom reference effect size (significant effect)

% source: http://www.ttable.org/uploads/2/1/7/9/21795380/9754276.png?852

clear bf_log10 res_Bayes bf bayesFactor
for iRoi =1:6
xdat = coupling_strenght_roi(iRoi,:)-median_cs_gray

xref = mean(xdat)/std(xdat)
xrefIroi(iRoi) = xref
[bf_log10(iRoi)]= my_ttest_bayes(xdat, xref);
[res_Bayes{iRoi}, bf(iRoi)] = interpret_Bayes(bf_log10(iRoi))
bayesFactor(iRoi) = 1/(10^bf(iRoi))

end




nobs = 30; % number of observations
% xref = +1.699; % significant t value for one-sided ttest 29 degrees of freedom reference effect size (significant effect)
xref = +2.042; % significant t value for two-sided ttest 29 degrees of freedom reference effect size (significant effect)


figure
bar(mean(coupling_strenght_roi,2)-mean(median_cs_gray))
legend

xticks([1:6])
xticklabels({'rPI','rDAI','rVAI','lPI','lDAI','lVAI'})
ylabel('Coupling strenght')


clear bf_log10 res_Bayes bf bayesFactor

% xref = +2.042; % significant t value for two-sided ttest 29 degrees of freedom reference effect size (significant effect)
% xref = +1.699; % significant t value for one-sided ttest 29 degrees of freedom reference effect size (significant effect)

for iRoi =1:6
xdat = coupling_strenght_roi(iRoi,:)-median_cs_gray

xref = mean(xdat)/std(xdat)
xrefIroi(iRoi) = xref
[bf_log10(iRoi)]= my_ttest_bayes(xdat, xref);
[res_Bayes{iRoi}, bf(iRoi)] = interpret_Bayes(bf_log10(iRoi))
bayesFactor(iRoi) = 1/(10^bf(iRoi))

end

figure
bar(bf)

figure
bar(bf_log10)

xticks([1:6])
xticklabels({'rPI','rDAI','rVAI','lPI','lDAI','lVAI'})
ylabel('BF log 10')


%% dummy test

clear bf_log10 res_Bayes bf bayesFactor
% xdat = coupling_strenght_roi(iRoi,:)-median_cs_gray


xdat = randn(1,100000);
xref = mean(xdat)/std(xdat)
% xrefIroi(iRoi) = xref
[bf_log10]= my_ttest_bayes(xdat, xref);
[res_Bayes, bf] = interpret_Bayes(bf_log10)
bayesFactor = 1/(10^bf)


xdat = randn(1,100);
xdat= xdat +1;
xref = mean(xdat)/std(xdat)
% xrefIroi(iRoi) = xref
[bf_log10]= my_ttest_bayes(xdat, xref);
[res_Bayes, bf] = interpret_Bayes(bf_log10)
bayesFactor = 1/(10^bf)

%%