%{
This script will 

1- Load Yeo 7 parcellation
2- Load the map of the significant gastric network regions, projected into the cortical surface
3 - Compute percentage of overlap per RSN and plot
    Here % overlap means % of gastric network in each RSN


Ignacio Rebollo
04-05-2021

%}

%% compute overlap RSNs


% First right hemisphere
% Load Yeo freesurfer native parcellation
filename = '/media/irebollo/DATAPART1/SADREST/files/rh.Yeo2011_7Networks_N1000.annot';
[verticesYeo, labelYeo, colortableYeo] = read_annotation(filename);

% load significant gastric network
[gasnet, M, mr_parms, volsz]=load_mgh('/media/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/rh.gasnetSURF.mgh');
gasnet_t = gasnet>0.01;


% number of RSNs
nROI = length(unique(labelYeo));
allROI = colortableYeo.table(:,5);

% Get vertices of each RSN
ROI_vertex = cell(nROI,1);
for iROI=1:nROI
ROI_vertex{iROI}=find(labelYeo ==allROI(iROI));
end

OverlapRH = zeros(nROI,1);
for iROI=1:nROI
    % find index of all vertex of that ROI
  OverlapRH(iROI)=  (sum(gasnet_t(ROI_vertex{iROI}))/sum(gasnet_t))*100;
end


% Then left hemisphere
filename = '/media/irebollo/DATAPART1/SADREST/files/lh.Yeo2011_7Networks_N1000.annot';
[verticesYeo, labelYeo, colortableYeo] = read_annotation(filename);

% loaf gastric network in surface
[gasnet, M, mr_parms, volsz]=load_mgh('/media/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/lh.gasnetSURF.mgh');
gasnet_t = gasnet>0.01;


nROI = length(unique(labelYeo));
% allROI = unique(labelHCP);

allROI = colortableYeo.table(:,5);

ROI_vertex = cell(nROI,1);
for iROI=1:nROI
ROI_vertex{iROI}=find(labelYeo ==allROI(iROI));
end

OverlapLH = zeros(nROI,1);
for iROI=1:nROI
    % find index of all vertex of that ROI
  OverlapLH(iROI)=  (sum(gasnet_t(ROI_vertex{iROI}))/sum(gasnet_t))*100;
  
end


meanOverlap = (OverlapLH + OverlapRH) /2;


meanOverlap = meanOverlap';
OOverlapRSNs = meanOverlap;

%% Plot 

labels = {'Visual', 'Sensory-Motor', 'Attention', 'Saliency','Limbic','Control','Default'}
colors = cbrewer('div', 'RdBu', 7);
dat=[]
dat(:, 1) = OOverlapRSNs(2:8)
[bR idxR] = sort(dat,'descend')
VISR = figure
hold on;
% if we want each bar to have a different color, loop
for b = 1:length(dat)
    bar(b*1.2, bR(b), 'FaceColor',  colors(idxR(b), : ), 'EdgeColor', 'k', 'BarWidth', 1);
end
hold on
    bar(1.2*8, OOverlapRSNs(1), 'FaceColor', 'k', 'EdgeColor', 'k', 'BarWidth', 1);
labelsPlotR =labels(idxR)
labelsPlotR{8}='Subcortical'
grid on
set(gca,'XTick',[1.2:1.2:8*1.2])
set(gca,'XTickLabel',labelsPlotR)
ylim([0 40])
title('Overlap with Resting State Networks')
ylabel('% Overlap')
set(gca,'XTickLabelRotation',45)
set(gca,'FontSize',14)
grid on
h= gcf
set(h,'PaperOrientation','landscape');
set(h,'PaperPositionMode','auto');  
set(gcf,'units','normalized','outerposition',[0 0 0.5 0.35])
print ('-dpdf', '-painters', '/media/Navigastric/scripts/plots&figures/OverlapRSNs')
    
    
    
