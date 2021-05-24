%{
This script will 

1- Load the cohen's d maps (+ bootstrapped variability) projected into the cortical surface
2- Load Yeo 7 parcellation
3 - Compute effecti size per RSN and plot


Ignacio Rebollo
03-05-2021

%}

%% Compute effect size per RSN

% Load effect sizes projected to cortical surface

[vol, M, mr_parms, volsz] = load_mgh('/media/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/rh.gasnetCD.mgh');

filename = '/media/irebollo/DATAPART1/SADREST/files/rh.Yeo2011_7Networks_N1000.annot';
[verticesYeo, labelYeo, colortableYeo] = read_annotation(filename);

[gasnet, M, mr_parms, volsz]=load_mgh('/media/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/rh.gasnetSURF.mgh');
gasnet_t = gasnet>0.01;
vol(~gasnet_t)=nan;


nROI = length(unique(labelYeo));% Get number of RSNs 


allROI = colortableYeo.table(:,5);
labelYeo(~gasnet_t)=0; % remove vertex outside the gastric network



ROI_vertex = cell(nROI,1); % find all the vertex of each rsns and store it in ROI_vertex
for iROI=1:nROI
ROI_vertex{iROI}=find(labelYeo ==allROI(iROI));
end

EffectSizeRH = zeros(nROI,1); % find the overlap 
for iROI=1:nROI
    % find index of all vertex of that ROI
  EffectSizeRH(iROI)=nanmean(vol(ROI_vertex{iROI})); 
end

%% Left hemisphere


[vol, M, mr_parms, volsz] = load_mgh('/media/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/lh.gasnetCD.mgh');

filename = '/media/irebollo/DATAPART1/SADREST/files/rh.Yeo2011_7Networks_N1000.annot';
[verticesYeo, labelYeo, colortableYeo] = read_annotation(filename);

[gasnet, M, mr_parms, volsz]=load_mgh('/media/Navigastric/scripts/bashScripts/freesurfer/warpGasnet/lh.gasnetSURF.mgh');
gasnet_t = gasnet>0.01;
vol(~gasnet_t)=nan;


nROI = length(unique(labelYeo));
% allROI = unique(labelHCP);

allROI = colortableYeo.table(:,5);

labelYeo(~gasnet_t)=0;
ROI_vertex = cell(nROI,1);
for iROI=1:nROI
ROI_vertex{iROI}=find(labelYeo ==allROI(iROI));
end

EffectSizeLH = zeros(nROI,1);
for iROI=1:nROI
    % find index of all vertex of that ROI
  EffectSizeLH(iROI)=nanmean(vol(ROI_vertex{iROI}));
  
end


%% Make figure



meanEffectSize = (EffectSizeLH + EffectSizeRH) /2;
meanEffectSize = meanEffectSize';
meanEffectSize_order = meanEffectSize(2:8);
meanEffectSize_order = meanEffectSize_order(idxR);



[bR idxR] = sort(dat,'descend')

meanEffectSize = meanEffectSize(idxR);

VISR = figure
hold on;
% if we want each bar to have a different color, loop
for b = 1:length(dat)
    bar(b*1.2, meanEffectSize_order(b), 'FaceColor',  colors(idxR(b), : ), 'EdgeColor', 'k', 'BarWidth', 1);
h = ploterr(b*1.2, meanEffectSize_order(b), [], meanOverlapSTD_order(b)/sqrt(63), 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none'); % remove marker
end
hold on
    bar(1.2*8, meanEffectSize(1), 'FaceColor', 'k', 'EdgeColor', 'k', 'BarWidth', 1);
h = ploterr(8*1.2, meanEffectSize(1), [], meanOverlapSTD(1)/sqrt(63), 'k.', 'abshhxy', 0);
labelsPlotR =labels(idxR)
labelsPlotR{8}='Subcortical'
grid on
set(gca,'XTick',[1.2:1.2:8*1.2])
set(gca,'XTickLabel',labelsPlotR)
ylim([0 40])
title('Effect sizes across Resting State Networks')
ylabel('Cohen D')
set(gca,'XTickLabelRotation',45)
set(gca,'FontSize',14)
grid on
h= gcf
set(h,'PaperOrientation','landscape');
set(h,'PaperPositionMode','auto');  
set(gcf,'units','normalized','outerposition',[0 0 0.5 0.35])
print ('-dpdf', '-painters', '/media/Navigastric/scripts/plots&figures/EffectSizesRSNs')

