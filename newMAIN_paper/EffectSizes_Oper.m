
%%
LOS_all = {'PFcm','OP4','OP1','OP2-3','43','FOP1','FOP2','FOP3','FOP4','PSL','RI','Ig','PoI1','PoI2','MI','52','PBelt','LBelt','A1','MBelt','A4','A5','TA2'}

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


colors = cbrewer('div', 'RdBu', length(LOS_all));

%% plot effect sizes according to hemisphere

dat=[]
dat(:, 1) = cd_all_R
dat(:, 2) = cd_all_L
% subplot(4,7,6); % rather than a square plot, make it thinner
figure
hold on;
% if we want each bar to have a different color, loop
for b = 1:length(dat)
    bar(b*2, dat(b,1), 'FaceColor',  colors(b, : ), 'EdgeColor', 'k', 'BarWidth', 0.8,'LineStyle','--');
    bar((b*2)-0.8, dat(b,2), 'FaceColor',  colors(b, : ), 'EdgeColor', 'k', 'BarWidth', 0.8);
end

for b = 1:length(dat)
    bar(b*2, dat(b,1), 'FaceColor',  colors(b, : ), 'EdgeColor', 'k', 'BarWidth', 0.8,'LineStyle','--');
    bar((b*2)-0.8, dat(b,2), 'FaceColor',  colors(b, : ), 'EdgeColor', 'k', 'BarWidth', 0.8);
end
labelsPlot =LOS_all
set(gca,'XTick',[1.5:2:length(dat)*2+1.5])
set(gca,'XTickLabel',labelsPlot)

% set(gca,'FontSize',10)
set(gca,'XTickLabelRotation',45)

title('Effect sizes in rolandic region')
ylabel (' Cohen D')

%% plot overlap

overlapOper=[0.232695000000000,0.131356000000000;0.264967000000000,0.0976277000000000;0.184211000000000,0.0470085000000000;0.165414000000000,0.213028000000000;0.708703000000000,0.143048000000000;0.00668896000000000,0.320000000000000;0,0.0825688000000000;0,0.0244444000000000;0,0.0205198000000000;0,0.0513113000000000;0.550055000000000,0.127509000000000;0.231618000000000,0.552147000000000;0,0.0530565000000000;0,0.200734000000000;0,0.0443262000000000;0.0808679000000000,0.466321000000000;0.130201000000000,0.362725000000000;0.236887000000000,0.150985000000000;0.242500000000000,0.156146000000000;0.168333000000000,0.530837000000000;0.0794492000000000,0.120468000000000;0.0923913000000000,0.168712000000000;0.0168067000000000,0.0931677000000000]
dat=[]
dat(:, 1) = overlapOper(:,2)*100
dat(:, 2) = overlapOper(:,1)*100
% subplot(4,7,6); % rather than a square plot, make it thinner
f1=figure
hold on;
% if we want each bar to have a different color, loop
for b = 1:length(dat)
    bar(b*2, dat(b,1), 'FaceColor',  colors(b, : ), 'EdgeColor', 'k', 'BarWidth', 0.8,'LineStyle','--');
    bar((b*2)-0.8, dat(b,2), 'FaceColor',  colors(b, : ), 'EdgeColor', 'k', 'BarWidth', 0.8);
end
labelsPlot =LOS_all
set(gca,'XTick',[1.5:2:length(dat)*2+1.5])
set(gca,'XTickLabel',labelsPlot)

% set(gca,'FontSize',10)
set(gca,'XTickLabelRotation',45)

title('Overlap with gastric network in Operculum')
ylabel (' % region in gastric network')




figure(f1)
    set(gcf,'units','normalized','outerposition',[0 0 1 0.3])
        set(gca,'FontSize',14)

grid on
    print ('-dpng', '-painters', '/media/Navigastric/scripts/plots&figures/0overlapOper')

mean(dat(:,2))
mean(dat(:,1))
    
    %%
    
    
cd_all_R(isnan(cd_all_R))=0
cd_all_L(isnan(cd_all_L))=0

    dat=[]
dat(:, 1) = cd_all_R
dat(:, 2) = cd_all_L

[bR idxR] = sort(cd_all_R,'descend')
[bL idxL] = sort(cd_all_L,'descend')

% subplot(4,7,6); % rather than a square plot, make it thinner
OPER = figure
hold on;
% if we want each bar to have a different color, loop
for b = 1:length(dat)
    bar(b, bR(b), 'FaceColor',  colors(idxR(b), : ), 'EdgeColor', 'k', 'BarWidth', 0.8);
end
labelsPlotR =LOS_all(idxR)
set(gca,'XTick',[1:length(dat)])
set(gca,'XTickLabel',labelsPlotR)
ylim([0.2 0.9])
title('OPER right Effect sizes')
set(gca,'XTickLabelRotation',45)

OPEL = figure
hold on;
for b = 1:length(dat)
 bar((b), bL(b), 'FaceColor',  colors(idxL(b), : ), 'EdgeColor', 'k', 'BarWidth', 0.8);
end
labelsPlotL =LOS_all(idxL)
set(gca,'XTick',[1:length(dat)])
set(gca,'XTickLabel',labelsPlotL)
ylim([0.2 0.9])
title('OPER left Effect sizes')
set(gca,'XTickLabelRotation',45)


figure(OPEL)
    set(gcf,'units','normalized','outerposition',[0 0 1 0.3])
grid on
    set(gca,'FontSize',14)
    print ('-dpng', '-painters', '/media/Navigastric/scripts/plots&figures/EeffectSizeOPERL')

    figure(OPER)
        set(gcf,'units','normalized','outerposition',[0 0 1 0.3])
grid on
    set(gca,'FontSize',14)

    print ('-dpng', '-painters', '/media/Navigastric/scripts/plots&figures/EeffectSizeOPERR')



mean(dat(dat(:,2)>0,2))
mean(dat(:,1))



