% % compare heart 2 samples

subjs_phys = [8 9 10 13 15 16 17 18 19 20 21 22 23 24 25 26 29 31 32 33 34 35 36 38 39 40 41 43 44];
% subj 37 not there?
% subjs_phys = [8 9 10 13 15 16 17 18 19 20 21 22 23 24 25 26 29 31 32 33 34 35 36 37 38 39 40 41 43 44];

% ['G:\HRVts\IBIts_S',sprintf('%.2d',subj_idx)];

subjs_navi = [2 3 4 9 10 12 13 15 16 17 18 19 21 22 23 25 28 29 30 31 32 33 35 36 38];
% ['I:\dataHeartNavi\IBIts_S',sprintf('%.2d',subj_idx)];

% subjs_all = [subjs_navi,subjs_phys]';


%% get PLV

% subjs_phys = [8 9 10 13 15 16 17 18 19 20 21 22 23 24 25 26 29 31 32 33 34 35 36 37 38 39 40 41 43 44];
% subjs_navi = [2 3 4 9 10 12 13 15 16 17 18 19 21 22 23 25 28 29 30 31 32 33 35 36 38 91];


subjectsPhysiens = subjs_phys
subjectsNavigastric = subjs_navi
cfgMain = global_getcfgmain
%%l oad all PLVs
indInside = tools_getIndexBrain('inside');
indOutside = tools_getIndexBrain('outside');



empirical = zeros(length(subjectsNavigastric)+length(subjectsPhysiens),153594); % Preallocate
surrogate = empirical; % for calculating t value



for iS=1:length(subjectsNavigastric)
    subj_idx = subjectsNavigastric(iS);
        cfgMain.sample = 1

    % empirical and chance PLV filenames
    
    filenamePLV = strcat(global_filename(subj_idx,cfgMain,strcat('PLVXVoxelFilename_',cfgMain.Timeseries2Regress)),'.nii');
       filenamePLVSurrogate  = strcat(global_path2root(1),'subj',sprintf('%.2d',subj_idx),'\Timeseries','\ControlEGG_othersubjectBOTHSAMPLES\median_sPLV_s',sprintf('%.2d',subj_idx),'.nii');
    
       
           PLVGroupEmpirical{iS} = ft_read_mri(filenamePLV); % Put into cell
    PLVGroupEmpirical{iS}.Nsubject = subjectsNavigastric(iS)+40;

    
    % Load surrogate PLV and prepare structure for surrogate PLV
    
    PLVGroupSurrogate{iS} = ft_read_mri(filenamePLVSurrogate);
    PLVGroupSurrogate{iS}.Nsubject = subjectsNavigastric(iS)+40;

    
    empirical(iS,:) = PLVGroupEmpirical{iS}.anatomy(:);
    surrogate(iS,:) = PLVGroupSurrogate{iS}.anatomy(:);
    

    empirical(iS,:) = PLVGroupEmpirical{iS}.anatomy(:);
    surrogate(iS,:) = PLVGroupSurrogate{iS}.anatomy(:);
    
    
end


for iS=1:length(subjectsPhysiens)
    
    cfgMain.sample = 2
    subj_idx = subjectsPhysiens(iS);
    
    % empirical and chance PLV filenames
    
    filenamePLV = strcat(global_filename(subj_idx,cfgMain,strcat('PLVXVoxelFilename_',cfgMain.Timeseries2Regress)),'.nii');
    filenamePLVSurrogate  = strcat(global_path2root(2),'subj',sprintf('%.2d',subj_idx),'\Timeseries','\ControlEGG_othersubjectBOTHSAMPLES\median_sPLV_s',sprintf('%.2d',subj_idx),'.nii');
    
    
    
    % Load empirical PLV
    
    PLVGroupEmpirical{iS+length(subjectsNavigastric)} = ft_read_mri(filenamePLV); % Put into cell
    PLVGroupEmpirical{iS+length(subjectsNavigastric)}.Nsubject = subjectsPhysiens(iS);
    

    % Load surrogate PLV and prepare structure for surrogate PLV
    
    PLVGroupSurrogate{iS+length(subjectsNavigastric)} = ft_read_mri(filenamePLVSurrogate);
    PLVGroupSurrogate{iS+length(subjectsNavigastric)}.Nsubject = subjectsPhysiens(iS);

    
    empirical(iS+length(subjectsNavigastric),:) = PLVGroupEmpirical{iS+length(subjectsNavigastric)}.anatomy(:);
    surrogate(iS+length(subjectsNavigastric),:) = PLVGroupSurrogate{iS+length(subjectsNavigastric)}.anatomy(:);
end


gasnetmask = global_getGastricNetwork;% = ft_read_mri ('Z:\ClusterResults\All\AAAAAAAAsurrSubjectsBothSample_nR10000_CA0050_kw3_fir2_fspread_015_fOrder_5_tw_15csfr_mask.nii');
CS_bothsamples = [empirical(:,find(gasnetmask)) - surrogate(:,find(gasnetmask))] ;
CS_bothsamples_meanXSubject = mean(CS_bothsamples,2);

 cs_all =   empirical- surrogate;

%% get ibi

load([global_path2root_folder,'FrequencyHeartNaviPhys.mat'])
% 
% all_ibi_phys = zeros(length(subjs_phys),890);
% for iS = 1:length(subjs_phys)
%     subj_idx = subjs_phys(iS)
% load(['Z:\dataHeartPhys\IBIts_S' sprintf('%.2d',subj_idx) '.mat'])
% all_ibi_phys(iS,:) = ibi_int(1:890);
% end
% 
% 
% all_ibi_navi = zeros(length(subjs_navi),889);
% for iS = 1:length(subjs_navi)
%     subj_idx = subjs_navi(iS)
% load(['Z:\dataHeartNavi\IBIts_S' sprintf('%.2d',subj_idx) '.mat'])
% all_ibi_navi(iS,:) = ibi_int(1:889);
% end
% 
% 
% data = [all_ibi_phys];
% nVoxels = size(data,1);
% %load base fieldtrip data structure
% load(strcat('Z:\scripts\STOMACH_BRAIN\files\sampleFieldtripStruc.mat'))
% % Define fieldtrip structure
% channelStr=cell(nVoxels,1);
% for iVoxel = 1:nVoxels
%     channelList(iVoxel,1) = iVoxel;
%     channelStr(iVoxel) = cellstr(mat2str(iVoxel));
% end
% 
% dataStructure.hdr = EGG_downsampled.hdr;
% dataStructure.fsample = 1;
% dataStructure.time{1,1}  = [0:1:(size(data,2))-1];
% dataStructure.label = channelStr;
% dataStructure.cfg = EGG_downsampled.cfg;
% dataStructure.trial{1,1} = data;
% cfgWelch = [];
% cfgWelch.keeptrials = 'no';
% cfgWelch.lengthWindow = 120;
% cfgWelch.overlap = 6;
% len = dataStructure.fsample*cfgWelch.lengthWindow; % length of subtrials cfg.length s in samples
% dataStructure.sampleinfo=[1 max(dataStructure.time{1,1})*dataStructure.fsample];
% cfg = [];
% cfg.trl(:,1) = dataStructure.sampleinfo(1):(len/cfgWelch.overlap):dataStructure.sampleinfo(2)-len+1;%trial start in samples from begining of raw data
% cfg.trl(:,2) = dataStructure.sampleinfo(1)+len-1:(len/cfgWelch.overlap):dataStructure.sampleinfo(2);%trial ends in samples from begining of raw data
% cfg.trl(:,3) = 0; %offset of the trigger with respect to the trial
% data_trials = ft_redefinetrial(cfg,dataStructure);
% % Estimate spectrum
% cfg = [];
% cfg.method = 'mtmfft';
% cfg.taper = 'hanning';
% cfg.output = 'pow';
% cfg.pad = 1000;
% cfg.foilim = [1/cfgWelch.lengthWindow 0.4]; % 0 - 6 cpm
% cfg.keeptrials = cfgWelch.keeptrials;
% frequencyWelch_Phys = ft_freqanalysis(cfg,data_trials);
% 
% 
% data = [all_ibi_navi];
% nVoxels = size(data,1);
% %load base fieldtrip data structure
% load(strcat('Z:\scripts\STOMACH_BRAIN\files\sampleFieldtripStruc.mat'))
% % Define fieldtrip structure
% channelStr=cell(nVoxels,1);
% for iVoxel = 1:nVoxels
%     channelList(iVoxel,1) = iVoxel;
%     channelStr(iVoxel) = cellstr(mat2str(iVoxel));
% end
% 
% dataStructure.hdr = EGG_downsampled.hdr;
% dataStructure.fsample = 1;
% dataStructure.time{1,1}  = [0:1:(size(data,2))-1];
% dataStructure.label = channelStr;
% dataStructure.cfg = EGG_downsampled.cfg;
% dataStructure.trial{1,1} = data;
% cfgWelch = [];
% cfgWelch.keeptrials = 'no';
% cfgWelch.lengthWindow = 120;
% cfgWelch.overlap = 6;
% len = dataStructure.fsample*cfgWelch.lengthWindow; % length of subtrials cfg.length s in samples
% dataStructure.sampleinfo=[1 max(dataStructure.time{1,1})*dataStructure.fsample];
% cfg = [];
% cfg.trl(:,1) = dataStructure.sampleinfo(1):(len/cfgWelch.overlap):dataStructure.sampleinfo(2)-len+1;%trial start in samples from begining of raw data
% cfg.trl(:,2) = dataStructure.sampleinfo(1)+len-1:(len/cfgWelch.overlap):dataStructure.sampleinfo(2);%trial ends in samples from begining of raw data
% cfg.trl(:,3) = 0; %offset of the trigger with respect to the trial
% data_trials = ft_redefinetrial(cfg,dataStructure);
% % Estimate spectrum
% cfg = [];
% cfg.method = 'mtmfft';
% cfg.taper = 'hanning';
% cfg.output = 'pow';
% cfg.pad = 1000;
% cfg.foilim = [1/cfgWelch.lengthWindow 0.4]; % 0 - 6 cpm
% cfg.keeptrials = cfgWelch.keeptrials;
% frequencyWelch_Navi = ft_freqanalysis(cfg,data_trials);

%% get power for GLM not zscored


ind_lf = find(frequencyWelch_Phys.freq > 0.05 & frequencyWelch_Phys.freq < 0.16)
ind_hf = find(frequencyWelch_Phys.freq > 0.15 & frequencyWelch_Phys.freq < 0.41)
ind_vlf = find(frequencyWelch_Phys.freq > 0.019 & frequencyWelch_Phys.freq < 0.065)
freq4zscore = frequencyWelch_Phys.freq > 0.05 & frequencyWelch_Phys.freq < 0.41


power_all = [frequencyWelch_Navi.powspctrm',frequencyWelch_Phys.powspctrm']


%% link with cs and distributions


std_power_all = std(power_all,[],2)./sqrt(size(power_all,2));

HF_power_all = mean(power_all(ind_hf,:))';
LF_power_all = mean(power_all(ind_lf,:))';
VLF_power_all = mean(power_all(ind_vlf,:))';
Allpower_all = mean(power_all(ind_lf(1):ind_hf(end),:))';

figure
nhist(HF_power_all,'median')
title(num2str(std(HF_power_all)*2+mean(HF_power_all)))

a=[];
a.HF = HF_power_all;
a.lf= LF_power_all;
a.vlf = VLF_power_all;
a.allp = Allpower_all;

figure
nhist(a,'binfactor',10)


threshold_hf= mean(HF_power_all)+2*std(HF_power_all);
threshold_lf= mean(LF_power_all)+2*std(LF_power_all);
outliers = HF_power_all>threshold_hf | LF_power_all >threshold_lf

figure
plot(frequencyWelch_Navi.freq(freq4zscore),power_all(freq4zscore,(~outliers))','b')

figure
plot(frequencyWelch_Navi.freq(freq4zscore),mean(power_all(freq4zscore,(~outliers)),2),'b')


% 
% figure;nhist(LF_power_all(~outliers),'binfactor',10,'color',[0 0.5 1]);title('LF no outliers')
% figure;nhist(HF_power_all(~outliers),'binfactor',10,'color',[0 1 0.5]);title('HF no outliers')



figure;nhist(LF_power_all(~outliers),'binfactor',10,'color',[0 1 0.5]);title('LF')
vline(std(LF_power_all(~outliers))*2+mean(LF_power_all(~outliers)),'k-','threshold outliers')

figure;nhist(HF_power_all(~outliers),'binfactor',10,'color',[0 0.5 1]);title('HF no outliers')
vline(std(HF_power_all(~outliers))*2+mean(HF_power_all(~outliers)),'k-','threshold outliers')



power_all_thresh =power_all(:,~outliers);
HFpowerallt = mean(power_all(ind_hf,~outliers));
LFpowerallt = mean(power_all(ind_lf,~outliers));
Allpower_all_t = mean(power_all(ind_lf(1):ind_hf(end),~outliers));


figure
plot(CS_bothsamples_meanXSubject(~outliers),HFpowerallt,'ok')
[r,p]=corrcoef(CS_bothsamples_meanXSubject(~outliers),HFpowerallt)
lsline
title(['wo hf r ',num2str(r(3)),' p ',num2str(p(3))])


figure
plot(LFpowerallt,CS_bothsamples_meanXSubject(~outliers),'ok')
[r,p]=corrcoef(LFpowerallt,CS_bothsamples_meanXSubject(~outliers))
lsline
title(['wo Lf r ',num2str(r(3)),' p ',num2str(p(3))])


figure
nhist(HF_power_all,'median','binfactor',10,'color',[0 0.5 1])
title('HF power')
vline(std(HF_power_all)*2+mean(HF_power_all),'k-','threshold outliers')

figure
nhist(LF_power_all,'median','binfactor',10,'color',[0 1 0.5])
title('LF power')
vline(std(LF_power_all)*2+mean(LF_power_all),'k-','threshold outliers')

figure
plot(HF_power_all,CS_bothsamples_meanXSubject,'ok')
[r,p]=corrcoef(CS_bothsamples_meanXSubject(~outliers),HFpowerallt)
lsline
title(['with o hf r ',num2str(r(3)),' p ',num2str(p(3))])

figure
plot(LF_power_all,CS_bothsamples_meanXSubject,'ok')
[r,p]=corrcoef(CS_bothsamples_meanXSubject(~outliers),HFpowerallt)
lsline
title(['with outlier  lf r ',num2str(r(3)),' p ',num2str(p(3))])



%% PLOT POWER


figure
plot(frequencyWelch_Navi.freq(freq4zscore),power_all(freq4zscore,:)','k');title('Spectrum all')
a=gca
set(a,'fontsize',12)

figure
plot(frequencyWelch_Navi.freq(freq4zscore),power_all(freq4zscore,~outliers)','k');title('Spectrum no outliers')
hold on
plot(frequencyWelch_Navi.freq(freq4zscore),mean(power_all(freq4zscore,~outliers),2),'r','linewidth',4);title('Spectrum mean no outliers')
a=gca
set(a,'fontsize',12)


%% Plot correlation of HRV with CS ins GASNET in subjects without outliers

figure
plot(CS_bothsamples_meanXSubject(~outliers),HFpowerallt,'ob','Linewidth',3)
[r,p]=corrcoef(CS_bothsamples_meanXSubject(~outliers),HFpowerallt)
lsline
title(['wo hf r ',num2str(r(3)),' p ',num2str(p(3))])
figure
plot(LFpowerallt,CS_bothsamples_meanXSubject(~outliers),'og','Linewidth',3)
[r,p]=corrcoef(LFpowerallt,CS_bothsamples_meanXSubject(~outliers))
lsline
title(['wo lf r ',num2str(r(3)),' p ',num2str(p(3))])

% plot 2 without more outliers
cs_thresh4outliers1 = CS_bothsamples_meanXSubject(~outliers);
% new_outliers = LF_power_all_t<2.5e-05 & HF_power_all_t<1.25e-05
% new_outliers = LF_power_all_t<4e-05 & HF_power_all_t<2e-05
% new_outliers = new_outliers'
new_outliers = LFpowerallt<std(LFpowerallt)*2+mean(LFpowerallt) & HFpowerallt<std(HFpowerallt)*2+mean(HFpowerallt)
% new_outliers = LF_power_all_t<4e-05 & HF_power_all_t<2e-05
new_outliers = new_outliers'
% vline(std(LF_power_all(~outliers))*2+mean(LF_power_all(~outliers)),'k-','threshold outliers')

new_outliers_HF =HFpowerallt<std(HFpowerallt)*2+mean(HFpowerallt)
new_outliers_LF =LFpowerallt<std(LFpowerallt)*2+mean(LFpowerallt)


figure
plot(HFpowerallt(new_outliers_HF),cs_thresh4outliers1(new_outliers_HF),'og')
[r,p]=corrcoef(HFpowerallt(new_outliers_HF),cs_thresh4outliers1(new_outliers_HF))
lsline
hold on
plot(HFpowerallt(new_outliers_HF),cs_thresh4outliers1(new_outliers_HF),'ok')
title(['no HF r ',num2str(r(3)),' p ',num2str(p(3))])
ylabel('cs')
xlabel('HF power')
set(gca,'fontsize',14)


figure
plot(LFpowerallt(new_outliers_LF),cs_thresh4outliers1(new_outliers_LF),'ob')
hold on
[r,p]=corrcoef(LFpowerallt(new_outliers_LF),cs_thresh4outliers1(new_outliers_LF))
lsline
plot(LFpowerallt(new_outliers_LF),cs_thresh4outliers1(new_outliers_LF),'ok')

title(['no lf r ',num2str(r(3)),' p ',num2str(p(3))])
ylabel('cs')
xlabel('LF power')
set(gca,'fontsize',14)




%corr btw h and l f

figure
plot(LFpowerallt(new_outliers),HFpowerallt(new_outliers),'or')
[r,p]=corrcoef(LFpowerallt(new_outliers),HFpowerallt(new_outliers))
lsline
title(['no hf lf r ',num2str(r(3)),' p ',num2str(p(3))])
ylabel('HF power')
xlabel('LF power')
set(gca,'fontsize',14)


figure
plot(Allpower_all_t(new_outliers),cs_thresh4outliers1(new_outliers),'ob')
[r,p]=corrcoef(Allpower_all_t(new_outliers),cs_thresh4outliers1(new_outliers))
lsline
title(['no allpower r ',num2str(r(3)),' p ',num2str(p(3))])
ylabel('CS')
xlabel('All power')
set(gca,'fontsize',14)

% corr with ratio
figure
plot(LFpowerallt(new_outliers)./HFpowerallt(new_outliers),cs_thresh4outliers1(new_outliers),'ok')
[r,p]=corrcoef(LFpowerallt(new_outliers)./HFpowerallt(new_outliers),cs_thresh4outliers1(new_outliers))
lsline
title(['no ratio r ',num2str(r(3)),' p ',num2str(p(3))])
ylabel('CS')
xlabel('Ratio LF to HF')

set(gca,'fontsize',14)





HF_power_all_c =HFpowerallt(new_outliers)'-mean(HFpowerallt(new_outliers))
LF_power_all_c = LFpowerallt(new_outliers)'-mean(LFpowerallt(new_outliers))'


HF_power_all_c =zscore(HFpowerallt(new_outliers))'
LF_power_all_c = zscore(LFpowerallt(new_outliers))'

RATIO_c = LFpowerallt(new_outliers)./HFpowerallt(new_outliers) 
figure
nhist(RATIO_c)

% RATIO_c(13)=[]
RATIO_c =RATIO_c - nanmean(RATIO_c)
RATIO_c = RATIO_c'

figure
nhist(RATIO_c)
% 
% 
% hf_map = ft_read_mri('F:\navigastric\ClusterResults\demographics\models\HRV_HFonly_noratio_noZ_centered_no\spmT_0011.nii')
% lf_map = ft_read_mri('F:\navigastric\ClusterResults\demographics\models\HRV_LFonly_noratio_noZ_centered_no\spmT_0011.nii')
% hf_map_t = hf_map.anatomy>3.3;
% lf_map_t = lf_map.anatomy>3.3;

cs_all_no= cs_all(~outliers,:);
% 
% figure
% plot(HF_power_all_c,mean(cs_all_no(new_outliers,hf_map_t),2),'ob')
% [r,p]=corrcoef(HF_power_all_c,mean(cs_all_no(new_outliers,hf_map_t),2))
% lsline
% title(['no allpower r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% figure
% plot(LF_power_all_c,mean(cs_all_no(new_outliers,hf_map_t),2),'ob')
% [r,p]=corrcoef(LF_power_all_c,mean(cs_all_no(new_outliers,hf_map_t),2))
% lsline
% title(['no allpower r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% figure
% plot(LF_power_all_c,mean(cs_all_no(new_outliers,lf_map_t),2),'ob')
% [r,p]=corrcoef(LF_power_all_c,mean(cs_all_no(new_outliers,lf_map_t),2))
% lsline
% title(['no allpower r ',num2str(r(3)),' p ',num2str(p(3))])
% 



% figure
% plot(HF_power_all_t,CS_bothsamples_meanXSubject(~outliers),'ob','Linewidth',3)
% hold on
% [r,p]=corrcoef(HF_power_all_t,CS_bothsamples_meanXSubject(~outliers))
% lsline
% title(['wo hf r ',num2str(r(3)),' p ',num2str(p(3))])

%% plots for paper

figure
scatterhist(LFpowerallt,CS_bothsamples_meanXSubject(~outliers),'Kernel','on','Location','NorthEast',...
    'Direction','out','Color','k','LineStyle',{'-'},...
    'LineWidth',[2],'Marker','+od','MarkerSize',[4]);
hold on
plot(LFpowerallt(new_outliers_LF),cs_thresh4outliers1(new_outliers_LF),'o','color',[0.5 0.5 0.5])
[r,p]=corrcoef(CS_bothsamples_meanXSubject(~outliers),LFpowerallt)
lsline
[r2,p2]=corrcoef(cs_thresh4outliers1(new_outliers_LF),LFpowerallt(new_outliers_LF))
title(['wo lf r ',num2str(r(3)),' p ',num2str(p(3)),' ','no lf r ',num2str(r2(3)),' p ',num2str(p2(3))])
ylabel('cs')
xlabel('LF power')
set(gca,'fontsize',14)

[RHO,PVAL] = corr(CS_bothsamples_meanXSubject(~outliers),LFpowerallt','Type','Spearman');



figure
scatterhist(HFpowerallt,CS_bothsamples_meanXSubject(~outliers),'Kernel','on','Location','NorthEast',...
    'Direction','out','Color','k','LineStyle',{'-'},...
    'LineWidth',[2],'Marker','+od','MarkerSize',[4]);
hold on
plot(HFpowerallt(new_outliers_HF),cs_thresh4outliers1(new_outliers_HF),'o','color',[0.5 0.5 0.5])
[r,p]=corrcoef(CS_bothsamples_meanXSubject(~outliers),HFpowerallt)
lsline
[r2,p2]=corrcoef(cs_thresh4outliers1(new_outliers_HF),HFpowerallt(new_outliers_HF))
title(['wo Hf r ',num2str(r(3)),' p ',num2str(p(3)),' ','no lf r ',num2str(r2(3)),' p ',num2str(p2(3))])
ylabel('cs')
xlabel('HF power')
set(gca,'fontsize',14)

[RHO,PVAL] = corr(CS_bothsamples_meanXSubject(~outliers),HFpowerallt','Type','Spearman');
[RHO,PVAL] = corr(cs_thresh4outliers1(new_outliers_HF),HFpowerallt(new_outliers_HF)','Type','Spearman');


ratio_all = LFpowerallt./HFpowerallt;
outliers_ratio = ratio_all<std(ratio_all)*2+mean(ratio_all)
% corr with ratio
figure
scatterhist(ratio_all,CS_bothsamples_meanXSubject(~outliers),'Kernel','on','Location','NorthEast',...
    'Direction','out','Color','k','LineStyle',{'-'},...
    'LineWidth',[4],'Marker','+od','MarkerSize',[4]);hold on
plot(ratio_all(outliers_ratio),cs_thresh4outliers1(outliers_ratio),'o','color',[0.5 0.5 0.5])
[r,p]=corrcoef(CS_bothsamples_meanXSubject(~outliers),ratio_all)
[r2,p2]=corrcoef(ratio_all(outliers_ratio),cs_thresh4outliers1(outliers_ratio))
lsline
title(['ratio O r ',num2str(r(3)),' p ',num2str(p(3)),' ','no ratio r ',num2str(r2(3)),' p ',num2str(p2(3))])
ylabel('CS')
xlabel('Ratio LF to HF')


[RHO,PVAL] = corr(CS_bothsamples_meanXSubject(~outliers),ratio_all','Type','Spearman');
[RHO,PVAL] = corr(ratio_all(outliers_ratio)',cs_thresh4outliers1(outliers_ratio),'Type','Spearman');




% corr with all power HF+LF

Allpower_all_t
% corr with ratio
figure
scatterhist(Allpower_all_t,CS_bothsamples_meanXSubject(~outliers),'Kernel','on','Location','NorthEast',...
    'Direction','out','Color','k','LineStyle',{'-'},...
    'LineWidth',[4],'Marker','+od','MarkerSize',[4]);hold on
plot(Allpower_all_t,cs_thresh4outliers1,'o','color',[0.5 0.5 0.5])
[r,p]=corrcoef(CS_bothsamples_meanXSubject(~outliers),Allpower_all_t)

% [r2,p2]=corrcoef(ratio_all(outliers_ratio),cs_thresh4outliers1(outliers_ratio))
lsline
title(['Allpower r ',num2str(r(3)),' p ',num2str(p(3))])
xlabel('CS')
ylabel('Power LF and HF')
[RHO,PVAL] = corr(CS_bothsamples_meanXSubject(~outliers),Allpower_all_t','Type','Spearman');





% peak_frequencies - mean_amplitudes
[b,stats] = robustfit([LFpowerallt;HFpowerallt;ratio_all]',CS_bothsamples_meanXSubject(~outliers));

figure;
scatter(LFpowerallt,CS_bothsamples_meanXSubject(~outliers),'filled'); grid on; hold on
plot(LFpowerallt,b(1)+b(2)*LFpowerallt,'g','LineWidth',2);
ax = gca;
ax.FontSize = 16;
xlabel('LF power', 'FontSize', 18);
ylabel('Coupling strenght','FontSize', 18);
title(['LF Robustfit, p=' num2str(stats.p(2))], 'FontSize', 18);
lsline

figure;
scatter(HFpowerallt,CS_bothsamples_meanXSubject(~outliers),'filled'); grid on; hold on
plot(HFpowerallt,b(1)+b(3)*HFpowerallt,'g','LineWidth',2);
ax = gca;
ax.FontSize = 16;
xlabel('HF power', 'FontSize', 18);
ylabel('Coupling strenght','FontSize', 18);
title(['HF Robustfit, p=' num2str(stats.p(3))], 'FontSize', 18);
lsline

figure;
scatter(ratio_all,CS_bothsamples_meanXSubject(~outliers),'filled'); grid on; hold on
plot(ratio_all,b(1)+b(4)*ratio_all,'g','LineWidth',2);
ax = gca;
ax.FontSize = 16;
xlabel('LF/HF ratio', 'FontSize', 18);
ylabel('Coupling strenght','FontSize', 18);
title(['Ratio Robustfit, p=' num2str(stats.p(4))], 'FontSize', 18);
lsline


[b,stats] = robustfit([Allpower_all_t;ratio_all]',CS_bothsamples_meanXSubject(~outliers));
% [b,stats] = robustfit(Allpower_all_t,cs_thresh4outliers1);

figure
scatter(Allpower_all_t,cs_thresh4outliers1,'filled'); grid on; hold on
plot(Allpower_all_t,b(1)+b(2)*Allpower_all_t,'g','LineWidth',2);
lsline
title(['allpower p ',num2str(stats.p(2))])
ylabel('CS')
xlabel('All power')
set(gca,'fontsize',14)

figure
scatter(ratio_all,cs_thresh4outliers1,'filled'); grid on; hold on
plot(ratio_all,b(1)+b(3)*ratio_all,'g','LineWidth',2);
lsline
title(['ratio p ',num2str(stats.p(3))])
ylabel('CS')
xlabel('ratio')
set(gca,'fontsize',14)



[b,stats] = robustfit(HFpowerallt(new_outliers_HF),cs_thresh4outliers1(new_outliers_HF));
figure;
scatter(HFpowerallt(new_outliers_HF),cs_thresh4outliers1(new_outliers_HF),'filled'); grid on; hold on
plot(HFpowerallt(new_outliers_HF),b(1)+b(2)*HFpowerallt(new_outliers_HF),'g','LineWidth',2);
ax = gca;
ax.FontSize = 16;
xlabel('EGG peak frequency (Hz)', 'FontSize', 18);
ylabel('mean ampl. (microvolt)', 'FontSize', 18);
title(['Robustfit, p=' num2str(stats.p(2))], 'FontSize', 18);


[b,stats] = robustfit(HFpowerallt,cs_thresh4outliers1);
figure;
scatter(HFpowerallt,cs_thresh4outliers1,'filled'); grid on; hold on
plot(HFpowerallt,b(1)+b(2)*HFpowerallt,'g','LineWidth',2);
ax = gca;
ax.FontSize = 16;
xlabel('EGG peak frequency (Hz)', 'FontSize', 18);
ylabel('mean ampl. (microvolt)', 'FontSize', 18);
title(['Robustfit, p=' num2str(stats.p(2))], 'FontSize', 18);


% log normal test
figure
scatterhist(CS_bothsamples_meanXSubject(~outliers),log(Allpower_all_t),'Kernel','on','Location','NorthEast',...
    'Direction','out','Color','k','LineStyle',{'-'},...
    'LineWidth',[4],'Marker','+od','MarkerSize',[4]);hold on
plot(CS_bothsamples_meanXSubject(~outliers),log(Allpower_all_t),'o','color',[0.5 0.5 0.5])
[r,p]=corrcoef(CS_bothsamples_meanXSubject(~outliers),log(Allpower_all_t))

% [r2,p2]=corrcoef(ratio_all(outliers_ratio),cs_thresh4outliers1(outliers_ratio))
lsline
title(['Allpower r ',num2str(r(3)),' p ',num2str(p(3))])
xlabel('CS')
ylabel('Power LF and HF')

[RHO,PVAL] = corr(CS_bothsamples_meanXSubject(~outliers),Allpower_all_t','Type','Spearman');



[b,stats] = robustfit(ratio_all,CS_bothsamples_meanXSubject(~outliers));
figure
scatter(ratio_all,cs_thresh4outliers1,'filled'); grid on; hold on
plot(ratio_all,b(1)+b(2)*ratio_all,'g','LineWidth',2);
lsline
title(['ratio p ',num2str(stats.p(2))])
ylabel('CS')
xlabel('ratio')
set(gca,'fontsize',14)

%% perform bayes stats

[r,p]=corrcoef(CS_bothsamples_meanXSubject(~outliers),HFpowerallt)
r_HF = r(3)
bf_HF= corrbf(r_HF,52)

[r,p]=corrcoef(LFpowerallt,CS_bothsamples_meanXSubject(~outliers))
r_LF = r(3)
bf_LF= corrbf(r_LF,52)

[r,p]=corrcoef(LFpowerallt,CS_bothsamples_meanXSubject(~outliers))
r_LF = r(3)
bf_LF= corrbf(r_LF,52)





set(gca,'fontsize',14)
r_RA = r(3)
bf_RA= corrbf(r_RA,52)

%% export table to excell
subj_all = [subjs_navi,subjs_phys+100]
LF_power_all
HF_power_all
ratio_allexport =  LF_power_all./HF_power_all;
table4export = [subj_all',LF_power_all,HF_power_all,ratio_allexport,CS_bothsamples_meanXSubject]
table4export(:,2) = table4export(:,2)*10000
table4export(:,3) = table4export(:,3)*10000

%% Correlation HRV demographics
% 
% variables_demos = []
% 
% figure
% subplot(3,1,1)
% plot(variables_demos(~outliers,1),variables_demos(~outliers,4),'ok')
% lsline
% subplot(3,1,2)
% plot(variables_demos(~outliers,1),variables_demos(~outliers,5),'ok')
% lsline
% subplot(3,1,3)
% plot(variables_demos(~outliers,1),variables_demos(~outliers,6),'ok')
% lsline
% 
% figure
% % subplot(3,1,1)
% scatterhist(variables_demos(~outliers,1),variables_demos(~outliers,4),'Kernel','on','Location','NorthEast',...
%     'Direction','out','Color','k','LineStyle',{'-'},...
%     'LineWidth',[4],'Marker','+od','MarkerSize',[4]);
% lsline
% [r,p]=corrcoef(variables_demos(~outliers,1),variables_demos(~outliers,4))
% title(['LFHRV EGG_FREQ r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% figure
% scatterhist(variables_demos(~outliers,1),variables_demos(~outliers,5),'Kernel','on','Location','NorthEast',...
%     'Direction','out','Color','k','LineStyle',{'-'},...
%     'LineWidth',[4],'Marker','+od','MarkerSize',[4]);
% lsline
% [r,p]=corrcoef(variables_demos(~outliers,1),variables_demos(~outliers,5))
% title(['HFHRV EGG_FREQ r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% figure
% scatterhist(variables_demos(~outliers,1),variables_demos(~outliers,6),'Kernel','on','Location','NorthEast',...
%     'Direction','out','Color','k','LineStyle',{'-'},...
%     'LineWidth',[4],'Marker','+od','MarkerSize',[4]);
% lsline
% [r,p]=corrcoef(variables_demos(~outliers,1),variables_demos(~outliers,6))
% title(['RATIOHRV EGG_FREQ r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% 
% 
% % correlation with power
% outliers_power = variables_demos(~outliers,2)> 1000
% powernoout = variables_demos(~outliers,2)
% powernoout(outliers_power)=[];
% 
% heartdemos_noout = variables_demos(~outliers,[4:6])
% heartdemos_noout(outliers_power,:)=[];
% 
% figure
% % subplot(3,1,1)
% scatterhist(powernoout,heartdemos_noout(:,1),'Kernel','on','Location','NorthEast',...
%     'Direction','out','Color','k','LineStyle',{'-'},...
%     'LineWidth',[4],'Marker','+od','MarkerSize',[4]);
% lsline
% [r,p]=corrcoef(powernoout,heartdemos_noout(:,1))
% title(['LFHRV EGG_POW r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% figure
% % subplot(3,1,1)
% scatterhist(powernoout,heartdemos_noout(:,2),'Kernel','on','Location','NorthEast',...
%     'Direction','out','Color','k','LineStyle',{'-'},...
%     'LineWidth',[4],'Marker','+od','MarkerSize',[4]);
% lsline
% [r,p]=corrcoef(powernoout,heartdemos_noout(:,2))
% title(['HFHRV EGG_POW r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% figure
% scatterhist(powernoout,heartdemos_noout(:,3),'Kernel','on','Location','NorthEast',...
%     'Direction','out','Color','k','LineStyle',{'-'},...
%     'LineWidth',[4],'Marker','+od','MarkerSize',[4]);
% lsline
% [r,p]=corrcoef(powernoout,heartdemos_noout(:,3))
% title(['RATIOHRV EGG_POW r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% 
% figure
% % subplot(3,1,1)
% scatterhist(variables_demos(~outliers,3),variables_demos(~outliers,4),'Kernel','on','Location','NorthEast',...
%     'Direction','out','Color','k','LineStyle',{'-'},...
%     'LineWidth',[4],'Marker','+od','MarkerSize',[4]);
% lsline
% [r,p]=corrcoef(variables_demos(~outliers,3),variables_demos(~outliers,4))
% title(['LFHRV EGG_FREQ r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% figure
% scatterhist(variables_demos(~outliers,3),variables_demos(~outliers,5),'Kernel','on','Location','NorthEast',...
%     'Direction','out','Color','k','LineStyle',{'-'},...
%     'LineWidth',[4],'Marker','+od','MarkerSize',[4]);
% lsline
% [r,p]=corrcoef(variables_demos(~outliers,3),variables_demos(~outliers,5))
% title(['HFHRV EGG_FREQ r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% figure
% scatterhist(variables_demos(~outliers,3),variables_demos(~outliers,6),'Kernel','on','Location','NorthEast',...
%     'Direction','out','Color','k','LineStyle',{'-'},...
%     'LineWidth',[4],'Marker','+od','MarkerSize',[4]);
% lsline
% [r,p]=corrcoef(variables_demos(~outliers,3),variables_demos(~outliers,6))
% title(['RATIOHRV EGG_FREQ r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% 
% % fwd = []
% figure
% % subplot(3,1,1)
% scatterhist(fwd(~outliers),variables_demos(~outliers,4),'Kernel','on','Location','NorthEast',...
%     'Direction','out','Color','k','LineStyle',{'-'},...
%     'LineWidth',[4],'Marker','+od','MarkerSize',[4]);
% lsline
% [r,p]=corrcoef(fwd(~outliers),variables_demos(~outliers,4))
% title(['LFHRV fWD r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% figure
% scatterhist(fwd(~outliers),variables_demos(~outliers,5),'Kernel','on','Location','NorthEast',...
%     'Direction','out','Color','k','LineStyle',{'-'},...
%     'LineWidth',[4],'Marker','+od','MarkerSize',[4]);
% lsline
% [r,p]=corrcoef(fwd(~outliers),variables_demos(~outliers,5))
% title(['HFHRV fWD r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% figure
% scatterhist(fwd(~outliers),variables_demos(~outliers,6),'Kernel','on','Location','NorthEast',...
%     'Direction','out','Color','k','LineStyle',{'-'},...
%     'LineWidth',[4],'Marker','+od','MarkerSize',[4]);
% lsline
% [r,p]=corrcoef(fwd(~outliers),variables_demos(~outliers,6))
% title(['RATIOHRV FWD r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% %% zscored
% 
% freq4zscore = frequencyWelch_Phys.freq > 0.05 & frequencyWelch_Phys.freq < 0.41
% power_phys_z =zscore(frequencyWelch_Phys.powspctrm(:,freq4zscore),[],2);
% power_navi_z=zscore(frequencyWelch_Navi.powspctrm(:,freq4zscore),[],2);
% 
% 
% std_power_navi = std(frequencyWelch_Navi.powspctrm(:,freq4zscore),[],2)./sqrt(size(frequencyWelch_Navi.powspctrm(:,freq4zscore),2));
% std_power_phys = std(frequencyWelch_Phys.powspctrm(:,freq4zscore),[],2)./sqrt(size(frequencyWelch_Phys.powspctrm(:,freq4zscore),2));
% 
% power_phys_unit =frequencyWelch_Phys.powspctrm(:,freq4zscore)./std_power_phys
% power_navi_unit=frequencyWelch_Navi.powspctrm(:,freq4zscore)./std_power_navi
% 
% 
% power_all_unit = [power_phys_unit',power_navi_unit'];
% 
% power_all_unit_thresh =power_all_unit(:,~outliers);
% power_all_unit_thresh_2 = power_all_unit_thresh(:,new_outliers);
% 
% figure
% plot(frequencyWelch_Navi.freq(freq4zscore),mean(power_all_unit_thresh_2,2),'b')
% hold on
% title('HRV power')
% xlabel('freq')
% ylabel('STD unit power')
% set(gca,'fontsize',16)
% 
% 
% 
% 
% ratio_lf_hf_all_unit = mean(power_all_unit_thresh_2(1:length(ind_lf),:),1)./mean(power_all_unit_thresh_2(length(ind_lf):end,:),1)
% figure
% nhist(ratio_lf_hf_all_unit)
% 
% 
% % corr with ratio
% figure
% plot(ratio_lf_hf_all_unit,cs_thresh4outliers1(new_outliers),'ok')
% [r,p]=corrcoef(ratio_lf_hf_all_unit,cs_thresh4outliers1(new_outliers))
% lsline
% title(['no ratio r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% % corr with hf
% figure
% plot(mean(power_all_unit(1:length(ind_lf),:),1),CS_bothsamples_meanXSubject,'ok')
% % [r,p]=corrcoef(ratio_lf_hf_all_unit(ratio_lf_hf_all_unit<8),CS_bothsamples_meanXSubject(ratio_lf_hf_all_unit<8))
% lsline
% title(['no ratio r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% % corr with hf
% figure
% plot(mean(power_all_unit(length(ind_lf):end,:),1),CS_bothsamples_meanXSubject,'ok')
% % [r,p]=corrcoef(ratio_lf_hf_all_unit(ratio_lf_hf_all_unit<8),CS_bothsamples_meanXSubject(ratio_lf_hf_all_unit<8))
% lsline
% % title(['no ratio r ',num2str(r(3)),' p ',num2str(p(3))])
% 
% power_lf_unit = mean(power_all_unit_thresh_2(1:length(ind_lf),:),1)';
% power_hf_unit= mean(power_all_unit_thresh_2(length(ind_lf):end,:),1)';
% ratio_lf_hf_all_unit= ratio_lf_hf_all_unit';
% 
% figure
% plot(power_lf_unit,power_hf_unit,'ok')
% 
% 
% figure
% nhist(power_lf_unit)
% figure
% nhist(power_hf_unit)
% figure
% nhist(ratio_lf_hf_all_unit)
% 
% outliers3 = power_lf_unit<mean(power_lf_unit)+2*std(power_lf_unit) & power_hf_unit<mean(power_hf_unit)+2*std(power_hf_unit) & ratio_lf_hf_all_unit<mean(ratio_lf_hf_all_unit)+2*std(ratio_lf_hf_all_unit)
% 
% 
% 
% power2export =[power_hf_unit(outliers3)-mean(power_hf_unit(outliers3)),power_lf_unit(outliers3)-mean(power_lf_unit(outliers3)),ratio_lf_hf_all_unit(outliers3)-mean(ratio_lf_hf_all_unit(outliers3))]
% %% Get all parameters both samples
% power_all = [frequencyWelch_Phys.powspctrm',frequencyWelch_Navi.powspctrm']
% std_power_all = std(power_all,[],2)./sqrt(size(power_all,2));
% 
% HF_power_all = mean(power_all(ind_hf,:))';
% LF_power_all = mean(power_all(ind_lf,:))';
% VLF_power_all = mean(power_all(ind_vlf,:))';
% 
% figure
% plot(frequencyWelch_Navi.freq,mean(power_all,2))
% hold on
% % errorbar(frequencyWelch_Navi.freq,mean(power_all,2),std_power_all,'color',[0.7,0.25,0])
% plot(frequencyWelch_Navi.freq,mean(power_all,2),'r','Linewidth',10)
% title('HRV power')
% xlabel('freq')
% ylabel('power (squared ms)')
% set(gca,'fontsize',16)
% shg
% 
% 
% %% threshold
% % figure
% power_all = [frequencyWelch_Phys.powspctrm',frequencyWelch_Navi.powspctrm'];
% outliers=zeros(1,size(power_all,2));
% threshold=(mean(mean(power_all))+2*std(mean(power_all)));
% outliers=mean(power_all)>threshold;
% power_all_thresh =power_all(:,~outliers);
% threshold=(mean(mean(power_all_thresh))+2*std(mean(power_all_thresh)));
% outliers=mean(power_all)>threshold;
% power_all_thresh =power_all(:,~outliers);
% threshold=(mean(mean(power_all_thresh))+2*std(mean(power_all_thresh)));
% outliers=mean(power_all)>threshold;
% power_all_thresh =power_all(:,~outliers);
% 
% HFpowerallt = mean(power_all(ind_hf,~outliers));
% LFpowerallt = mean(power_all(ind_lf,~outliers));
% VLF_power_all_t = mean(power_all(ind_vlf,~outliers));
% % pLOT 
% 
% figure
% plot(CS_bothsamples_meanXSubject(~outliers),HFpowerallt,'ok')
% [r,p]=corrcoef(CS_bothsamples_meanXSubject(~outliers),HFpowerallt)
% lsline
% title(num2str(r(3)))
% figure
% plot(LFpowerallt,CS_bothsamples_meanXSubject(~outliers),'ob')
% [r,p]=corrcoef(LFpowerallt,CS_bothsamples_meanXSubject(~outliers))
% lsline
% title(num2str(r(3)))
% figure
% plot(VLF_power_all_t,CS_bothsamples_meanXSubject(~outliers),'or')
% [r,p]=corrcoef(VLF_power_all_t,CS_bothsamples_meanXSubject(~outliers))
% lsline
% title(num2str(r(3)))
% 
% figure
% plot(LFpowerallt./CS_bothsamples_meanXSubject(~outliers),HFpowerallt,'ob')
% [r,p]=corrcoef(LFpowerallt./CS_bothsamples_meanXSubject(~outliers),HFpowerallt)
% lsline
% title(num2str(r(3)))
% 
% 
% 
% %% power zscored
% power_phys_z =zscore(frequencyWelch_Phys.powspctrm,[],2);
% power_navi_z=zscore(frequencyWelch_Navi.powspctrm,[],2);
% % figure
% 
% power_all = [power_phys_z',power_navi_z'];
% outliers=zeros(1,size(power_all,2));
% threshold=(mean(mean(power_all))+2*std(mean(power_all)));
% outliers=mean(power_all)>threshold;
% % power_all_thresh =power_all(:,~outliers);
% % threshold=(mean(mean(power_all_thresh))+2*std(mean(power_all_thresh)));
% % outliers=mean(power_all)>threshold;
% % power_all_thresh =power_all(:,~outliers);
% % threshold=(mean(mean(power_all_thresh))+2*std(mean(power_all_thresh)));
% % outliers=mean(power_all)>threshold;
% % power_all_thresh =power_all(:,~outliers);
% 
% HFpowerallt = mean(power_all(ind_hf,~outliers));
% LFpowerallt = mean(power_all(ind_lf,~outliers));
% VLF_power_all_t = mean(power_all(ind_vlf,~outliers));
% % pLOT 
% 
% figure
% plot(CS_bothsamples_meanXSubject(~outliers),HFpowerallt,'ok')
% [r,p]=corrcoef(CS_bothsamples_meanXSubject(~outliers),HFpowerallt)
% lsline
% title(['z HF ',num2str(r(3)),' p ',num2str(p(3)) ])
% figure
% plot(LFpowerallt,CS_bothsamples_meanXSubject(~outliers),'ob')
% [r,p]=corrcoef(LFpowerallt,CS_bothsamples_meanXSubject(~outliers))
% lsline
% title(['z LF ',num2str(r(3)),' p ',num2str(p(3)) ])
% figure
% plot(VLF_power_all_t,CS_bothsamples_meanXSubject(~outliers),'or')
% [r,p]=corrcoef(VLF_power_all_t,CS_bothsamples_meanXSubject(~outliers))
% lsline
% title(['z VLF ',num2str(r(3)),' p ',num2str(p(3)) ])
% 
% figure
% plot(LFpowerallt./CS_bothsamples_meanXSubject(~outliers),HFpowerallt,'ob')
% [r,p]=corrcoef(LFpowerallt./CS_bothsamples_meanXSubject(~outliers),HFpowerallt)
% lsline
% title(['z ratio ',num2str(r(3)),' p ',num2str(p(3)) ])
% 
% 
% %% for spm
% 
% HF_power_all_c =HFpowerallt'-mean(HFpowerallt)
% LF_power_all_c = LFpowerallt'-mean(LFpowerallt)'
% VLF_power_all_c=VLF_power_all_t'-mean(VLF_power_all_t)'
% RATIO_c = HFpowerallt./LFpowerallt - nanmean(HFpowerallt./LFpowerallt)
% RATIO_c = RATIO_c'
% 
% figure
% nhist(RATIO_c)
% 
% 
% figure
% plot(RATIO_c,CS_bothsamples_meanXSubject,'ob')
% [r,p]=corrcoef(RATIO_c,CS_bothsamples_meanXSubject)
% lsline
% title(['z ratio ',num2str(r(3)),' p ',num2str(p(3)) ])
% 
% 
% figure
% plot(HF_power_all_c,RATIO_c,'ok')
% [r,p]=corrcoef(RATIO_c,HF_power_all_c)
% ylabel('Ratio')
% xlabel('HFpower')
% 
% title(['HF powerc,RATIO' ,num2str(r(3))])
% lsline
% figure
% plot(LF_power_all_c,RATIO_c,'ok')
% [r,p]=corrcoef(RATIO_c,LF_power_all_c)
% ylabel('Ratio')
% xlabel('LFpower')
% 
% title(['LF power,RATIO' ,num2str(r(3))])
% lsline
% 
% figure
% plot(LF_power_all_c,HF_power_all_c,'ok')
% [r,p]=corrcoef(LF_power_all_c,HF_power_all_c)
% xlabel('LFpower')
% ylabel('HFpower')
% 
% title(['LF power,HF power' ,num2str(r(3))])
% lsline
% 
% %% power no
% 
% %% check effect of regressing ratio
% 
% toBeExplained = [LF_power_all_c HF_power_all_c,CS_bothsamples_meanXSubject]; % BOLD timeseries will be the variable to be explained out in the GLM
% 
% % csf_timeseries = nanmean(BOLD_filtered_zscored(logical(csf),:)); % Average the timeseries in the csf compartment
% explainingVariables = [RATIO_c]; % Variables used to explain the BOLD data
% betas_regression = tools_EfficientGLM(toBeExplained,explainingVariables); % Obtain the betas indicating how much the predicting variable predicts the data
% predictedBOLD = explainingVariables(:,1) *betas_regression(1,:); % What the BOLD timeseries should look like if CSF predicted at 100% accuracy the data
% error_hrv = toBeExplained - predictedBOLD; % The error is the portion of the data not predicted by the CSF signal
% 
% figure
% plot(error_hrv(:,1),error_hrv(:,3),'ok');lsline
% [r,p]=corrcoef(error_hrv(:,1),error_hrv(:,3))
% title(['residuals LF power,CS' ,num2str(r(3))])
% 
% figure
% plot(error_hrv(:,1),error_hrv(:,2),'ok');lsline
% [r,p]=corrcoef(error_hrv(:,1),error_hrv(:,2))
% title(['residuals LF power,HF power' ,num2str(r(3))])
% 
% figure
% plot(error_hrv(:,2),error_hrv(:,3),'ok');lsline
% [r,p]=corrcoef(error_hrv(:,2),error_hrv(:,3))
% title(['residuals HF power,CS' ,num2str(r(3))])
% 
% 
% 
% 
% toBeExplained = [CS_bothsamples_meanXSubject]; % BOLD timeseries will be the variable to be explained out in the GLM
% explainingVariables = [HF_power_all_c ,LF_power_all_c,VLF_power_all_c,RATIO_c]; % Variables used to explain the BOLD data
% betas_regression = tools_EfficientGLM(toBeExplained,explainingVariables); % Obtain the betas indicating how much the predicting variable predicts the data
% predictedBOLD = explainingVariables(:,1)*betas_regression(1)+explainingVariables(:,2) *betas_regression(2)+explainingVariables(:,3) *betas_regression(3)+explainingVariables(:,4) *betas_regression(4); % What the BOLD timeseries should look like if CSF predicted at 100% accuracy the data
% error_cs= toBeExplained - predictedBOLD; % The error is the portion of the data not predicted by the CSF signal
% 
% figure
% plot(explainingVariables(:,1)*betas_regression(1),CS_bothsamples_meanXSubject,'ok')
% 
% figure
% plot(error_hrv(:,1),error_hrv(:,3),'ok');lsline
% [r,p]=corrcoef(error_hrv(:,1),error_hrv(:,3))
% title(['residuals LF power,CS' ,num2str(r(3))])
% 
% figure
% plot(error_hrv(:,1),error_hrv(:,2),'ok');lsline
% [r,p]=corrcoef(error_hrv(:,1),error_hrv(:,2))
% title(['residuals LF power,HF power' ,num2str(r(3))])
% 
% figure
% plot(error_hrv(:,2),error_hrv(:,3),'ok');lsline
% [r,p]=corrcoef(error_hrv(:,2),error_hrv(:,3))
% title(['residuals HF power,CS' ,num2str(r(3))])
% 
% 
% %% OLD
% %% plot figure poster
% 
% figure
% plot(frequencyWelch_Navi.freq,mean(zscore(frequencyWelch_Phys.powspctrm,[],2)),'r')
% [h p ci stats] = ttest2(zscore(frequencyWelch_Navi.powspctrm,[],2),zscore(frequencyWelch_Phys.powspctrm,[],2))
% errorbar(frequencyWelch_Navi.freq,mean(power_phys_z),std_power_phys,'r')
% title('HRV power')
% xlabel('freq')
% ylabel('Z power')
% set(gca,'fontsize',16)
% 
% 
% 
% figure
% plot(frequencyWelch_Navi.freq,mean(zscore(frequencyWelch_Phys.powspctrm,[],2)),'r')
% [h p ci stats] = ttest2(zscore(frequencyWelch_Navi.powspctrm,[],2),zscore(frequencyWelch_Phys.powspctrm,[],2))
% errorbar(frequencyWelch_Navi.freq,mean(power_phys_z),std_power_phys,'r')
% title('HRV power')
% xlabel('freq')
% ylabel('Z power')
% set(gca,'fontsize',16)
% %% median split high and low and ttest
% 
% 
% 
% %% low high power
% % 
% % subjs_phys = [8 9 10 13 15 16 17 18 19 20 21 22 23 24 25 26 29 31 32 33 34 35 36 37 38 39 40 41 43 44];
% % subjs_navi = [2 3 4 9 10 12 13 15 16 17 18 19 21 22 23 25 28 29 30 31 32 33 35 36 38 91];
% 
% power_vlf_navi = mean(power_navi_z(:,ind_vlf),2)
% power_vlf_phys = mean(power_phys_z(:,ind_vlf),2)
% 
% 
% powervlf_all = [power_vlf_navi',power_vlf_phys']
% 
% %high power
% subjs_navi_hr = power_vlf_navi > median(powervlf_all)
% subjs_phys_hr = power_vlf_phys > median(powervlf_all)
% 
% subjNavi2Cluster_hr = subjs_navi(subjs_navi_hr)
% subjPhys2cluster_hr =subjs_phys(subjs_phys_hr)
% 
% %low power
% subjs_navi_lr = power_vlf_navi <= median(powervlf_all)
% subjs_phys_lr = power_vlf_phys <= median(powervlf_all)
% 
% subjNavi2Cluster_lr = subjs_navi(subjs_navi_lr)
% subjPhys2cluster_lr =subjs_phys(subjs_phys_lr)
% 
% cfgMain.clusterAlpha = 0.025
% 
% 
% cfgMain.numberofrandomizations = 1000
% timeseries_statsCluster_PLV2samplettest(subjs_navi(subjs_navi_hr),subjs_phys(subjs_phys_hr),subjs_navi(subjs_navi_lr),subjs_phys(subjs_phys_lr),cfgMain)
% 
% 
% 
% %% freq statistics
% % something is wrong with the design matrix
% 
% frequencyWelch_Navi_z = frequencyWelch_Navi;
% frequencyWelch_Navi_z.powspctrm = power_navi_z;
% frequencyWelch_Phys_z = frequencyWelch_Phys;
% frequencyWelch_Phys_z.powspctrm = power_phys_z;
% 
% cfg = [];
% cfg.statistic        = 'ft_statfun_indepsamplesT';
% cfg.method = 'montecarlo'
% cfg.tail             = 0;
% cfg.clustertail      = 0;
% cfg.alpha            = 0.025;
% cfg.numrandomization = 1000;
% cfgStats.clusteralpha = 0.025;
% 
% 
% 
% cfg.design(1,:) = [1:56];
% % cfgStats.design(1,:) = [1:length(subjects_highStai) 1:length(subjects_lowStai)];
% cfg.design(2,:) = [ones(1,size(frequencyWelch_Navi_z.powspctrm,1)) ones(1,size(frequencyWelch_Phys_z.powspctrm,1))*2];
% % cfgStats.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
% cfg.ivar        = 2; % row of design matrix that contains independent variable (the conditions)
% 
% % design = zeros(1,size(frequencyWelch_Navi_z.powspctrm,1) + size(frequencyWelch_Phys_z.powspctrm,1));
% % design(1,1:size(frequencyWelch_Navi_z.powspctrm,1)) = 1;
% % design(1,(size(frequencyWelch_Navi_z.powspctrm,1)+1):(size(frequencyWelch_Navi_z.powspctrm,1)+...
% % size(frequencyWelch_Phys_z.powspctrm,1))) = 2;
% % 
% % cfg.design           = design;
% % cfg.ivar             = 1;
% 
% 
% 
% stats = ft_freqstatistics(cfg,frequencyWelch_Navi_z,frequencyWelch_Phys_z)
% %% 
% ratio_all= [ratio_lf_hf_navi',ratio_lf_hf_phys']
% 
% 
% subjectsTableheart = [ratio_all;powervlf_all;subjs_navi,subjs_phys]
% subjectsTableheart = subjectsTableheart'
% % subjs_phys = [8 9 10 13 15 16 17 18 19 20 21 22 23 24 25 26 29 31 32 33 34 35 36 37 38 39 40 41 43 44];
% % subjs_navi = [2 3 4 9 10 12 13 15 16 17 18 19 21 22 23 25 28 29 30 31 32 33 35 36 38 91];
% 
% %high ratio
% subjs_navi_hr = ratio_lf_hf_navi > median(ratio_all)
% subjs_phys_hr = ratio_lf_hf_phys > median(ratio_all)
% 
% subjNavi2Cluster = subjs_navi(subjs_navi_hr)
% subjPhys2cluster =subjs_phys(subjs_phys_hr)
% 
% cfgMain.clusterAlpha = 0.025
% % timeseries_statsCluster_Regression_surrSubj_2s(subjNavi2Cluster,subjPhys2cluster,cfgMain)
% 
% %low ratio
% subjs_navi_hr = ratio_lf_hf_navi <= median(ratio_all)
% subjs_phys_hr = ratio_lf_hf_phys <= median(ratio_all)
% 
% subjNavi2Cluster = subjs_navi(subjs_navi_hr)
% subjPhys2cluster =subjs_phys(subjs_phys_hr)
% 
% cfgMain.clusterAlpha = 0.025
% % timeseries_statsCluster_Regression_surrSubj_2s(subjNavi2Cluster,subjPhys2cluster,cfgMain)
% 
% 
% %% nonz
% % 
% % [h p ci stats] = ttest2(frequencyWelch_Phys.powspctrm,frequencyWelch_Navi.powspctrm)
% % 
% % figure
% % plot(frequencyWelch_Navi.freq,mean(frequencyWelch_Navi.powspctrm),'b')
% % hold on
% % plot(frequencyWelch_Navi.freq,mean(frequencyWelch_Phys.powspctrm),'r')
% % plot(frequencyWelch_Navi.freq,h*0.1,'k')
% % 
% % 
% % 
% % ind_lf = find(frequencyWelch_Phys.freq > 0.05 & frequencyWelch_Phys.freq < 0.16)
% % ind_hf = find(frequencyWelch_Phys.freq > 0.15 & frequencyWelch_Phys.freq < 0.41)
% % 
% % 
% % ratio_lf_hf_phys = mean(frequencyWelch_Phys.powspctrm(:,ind_lf),2)./mean(frequencyWelch_Phys.powspctrm(:,ind_hf),2)
% % ratio_lf_hf_navi = mean(frequencyWelch_Navi.powspctrm(:,ind_lf),2)./mean(frequencyWelch_Navi.powspctrm(:,ind_hf),2)
% % 
% % 
% % [h p ci stats] = ttest2(ratio_lf_hf_phys,ratio_lf_hf_navi)
% % 
% % 
% % figure
% % bar([mean(ratio_lf_hf_phys);mean(ratio_lf_hf_navi)])
% % frequencyWelch_Phys.powspctrm
% 
% 
% 
% 
% lsline
% [r,p]=corrcoef(ratio_all,CS_bothsamples_meanXSubject)
% [r,p]=corrcoef(ratio_all(1:26),CS_bothsamples_meanXSubject(1:26))
% [r,p]=corrcoef(ratio_all(27:55),CS_bothsamples_meanXSubject(27:55))
% 
% title('r = -0.27 p = 0.042')
% xlabel('LF to HF ratio')
% ylabel('PLV in gasnet')
% hold on
% line([median(ratio_all),median(ratio_all)],[0.15 0.4])
% 
% notoutliers_ratio = find(ratio_all >-6 & ratio_all<6)
% ratio_all_nooutlier = ratio_all(notoutliers_ratio)
% figure
% plot(ratio_all_nooutlier,CS_bothsamples_meanXSubject(notoutliers_ratio),'ok')
% [r,p]=corrcoef(ratio_all_nooutlier,CS_bothsamples_meanXSubject(notoutliers_ratio))
% 
% 
% figure
% plot(powervlf_all,CS_bothsamples_meanXSubject,'ok')
% lsline
% hold on
% plot(powervlf_all(1:26),CS_bothsamples_meanXSubject(1:26),'ob')
% plot(powervlf_all(27:55),CS_bothsamples_meanXSubject(27:55),'or')
% % [r,p]=corrcoef(powervlf_all,PLV_bothsamples_meanXSubject)
% % title('r = -0.27 p = 0.042')
% xlabel('VLF power')
% ylabel('PLV in gasnet')
% hold on
% line([median(powervlf_all),median(powervlf_all)],[0.15 0.4])
% 
% 
% 
% 
%  [h p ci stats] = ttest2(CS_bothsamples_meanXSubject(ratio_all<median(ratio_all),:),CS_bothsamples_meanXSubject(ratio_all>median(ratio_all),:));
%  
% %% Distribution across participants
% figure
% subplot(2,1,1)
% violin(ratio_all')
% title('ratio_all')
% subplot(2,1,2)
% plot(ratio_all','ok')
% hold on
% line([26.5 26.5],[-10 10])
% 
% 
% figure
% subplot(2,1,1)
% violin({powervlf_all(1:26)',powervlf_all(27:56)'})
% title('powervlf_all')
% subplot(2,1,2)
% plot(powervlf_all','ok')
% hold on
% line([26.5 26.5],[-2 2])
% 
% 
% %% variables per voxel lf/hf ratio
% a = []
% a.highRatio =  mean(CS_bothsamples(ratio_all>median(ratio_all),:));
% a.lowRatio =  mean(CS_bothsamples(ratio_all<median(ratio_all),:));
% a.all =   mean(CS_bothsamples);
% [h p ci stats]= ttest2(a.highRatio,a.lowRatio);
% figure
% nhist(a,'median')
% title(['Mean CS  per voxel across participants p',num2str(p),' t ', num2str(stats.tstat)])
% 
% a = []
% a.highRatio =  var(CS_bothsamples(ratio_all>median(ratio_all),:));
% a.lowRatio =  var(CS_bothsamples(ratio_all<median(ratio_all),:));
% a.all =   var(CS_bothsamples);
% 
% figure
% nhist(a,'median')
% [h p ci stats]= ttest2(a.highRatio,a.lowRatio);
% title(['Variance CS per voxel across participants p',num2str(p),' t ', num2str(stats.tstat)])
% 
% 
% a = []
% a.highRatio =  mean(CS_bothsamples(ratio_all>median(ratio_all),:),2);
% a.lowRatio =  mean(CS_bothsamples(ratio_all<median(ratio_all),:),2);
% a.all =   mean(CS_bothsamples,2);
% 
% figure
% nhist(a,'median')
% [h p ci stats]= ttest2(a.highRatio,a.lowRatio);
% title(['Mean CS  per participants across voxels  p',num2str(p),' t ', num2str(stats.tstat)])
% 
% a = []
% a.highRatio =  var(CS_bothsamples(ratio_all>median(ratio_all),:),[],2);
% a.lowRatio =  var(CS_bothsamples(ratio_all<median(ratio_all),:),[],2);
% a.all=   var(CS_bothsamples,[],2);
% figure
% nhist(a,'median')
% [h p ci stats]= ttest2(a.highRatio,a.lowRatio);
% title(['Variance CS per participants across voxels  p',num2str(p),' t ', num2str(stats.tstat)])
% 
% %%
% a = []
% a.highPower =  mean(CS_bothsamples(powervlf_all>median(powervlf_all),:));
% a.lowPower =  mean(CS_bothsamples(powervlf_all<median(powervlf_all),:));
% a.all =   mean(CS_bothsamples);
% 
% [h p ci stats]= ttest2(a.highPower,a.lowPower);
% figure
% nhist(a,'median')
% title(['Mean CS  per voxel across participants p',num2str(p),' t ', num2str(stats.tstat)])
% 
% a = []
% a.highPower =  var(CS_bothsamples(powervlf_all>median(powervlf_all),:));
% a.lowPower =  var(CS_bothsamples(powervlf_all<median(powervlf_all),:));
% a.all =   var(CS_bothsamples);
% 
% figure
% nhist(a,'median')
% [h p ci stats]= ttest2(a.highPower,a.lowPower);
% title(['Variance CS per voxel across participants p',num2str(p),' t ', num2str(stats.tstat)])
% 
% 
% a = []
% a.highPower =  mean(CS_bothsamples(powervlf_all>median(powervlf_all),:),2);
% a.lowPower =  mean(CS_bothsamples(powervlf_all<median(powervlf_all),:),2);
% a.all =   mean(CS_bothsamples,2);
% 
% figure
% nhist(a,'median')
% [h p ci stats]= ttest2(a.highPower,a.lowPower);
% title(['Mean CS  per participants across voxels  p',num2str(p),' t ', num2str(stats.tstat)])
% 
% a = []
% a.highPower =  var(CS_bothsamples(powervlf_all>median(powervlf_all),:),[],2);
% a.lowPower =  var(CS_bothsamples(powervlf_all<median(powervlf_all),:),[],2);
% a.all=   var(CS_bothsamples,[],2);
% 
% figure
% nhist(a,'median')
% [h p ci stats]= ttest2(a.highPower,a.lowPower);
% title(['Variance CS per participants across voxels  p',num2str(p),' t ', num2str(stats.tstat)])
% 
% 
% %% perform GLM on ratio
% 
% 
% clusterMap = ft_read_mri ('Z:\ClusterResults\All\surrSubjectsBothSample_nR1000_CA0010_kw3_fir2_fspread_015_fOrder_5_tw_15csfr_ClusterMap.nii');
% clusterMap = clusterMap.anatomy(:);
% nclusters = max(clusterMap);
% 
% 
% nSubjects = size(coupling_strenght,1);
% CS_cluster_subject = zeros(nSubjects,nclusters);
% for iCluster = 1: nclusters
%     indCluster = find(clusterMap==iCluster);
%     for iSubj =1:nSubjects
%    CS_cluster_subject(iSubj,iCluster)= mean(coupling_strenght(iSubj,indCluster),2);
%     end
% end
% 
% r_ratio_cluster = zeros(nclusters,1);
% p_ratio_cluster = zeros(nclusters,1);
% 
% for iCluster = 1: nclusters
% [r,p]= corrcoef(CS_cluster_subject(:,iCluster),ratio_all);
% r_ratio_cluster(iCluster) = r(3)
% p_ratio_cluster(iCluster) = p(3)
% end
% r_ratio_cluster(find(p_ratio_cluster<0.05))
% p_ratio_cluster(find(p_ratio_cluster<0.05))
% 
% find(p_ratio_cluster<0.05)
% 
% 
% % vlf power
% r_vlf_cluster = zeros(nclusters,1);
% p_vlf_cluster = zeros(nclusters,1);
% 
% for iCluster = 1: nclusters
% [r,p]= corrcoef(CS_cluster_subject(:,iCluster),powervlf_all);
% r_vlf_cluster(iCluster) = r(3);
% p_vlf_cluster(iCluster) = p(3);
% end
% r_vlf_cluster(find(p_vlf_cluster<0.05))
% p_vlf_cluster(find(p_vlf_cluster<0.05))
% 
% find(p_vlf_cluster<0.05)
% 
% 
% %% export table for demographics
% 
% exportTable = [subjs_navi,subjs_phys;powervlf_all;ratio_all]'
% % 
% % subjs_phys = [8 9 10 13 15 16 17 18 19 20 21 22 23 24 25 26 29 31 32 33 34 35 36 37 38 39 40 41 43 44];
% % subjs_navi = [2 3 4 9 10 12 13 15 16 17 18 19 21 22 23 25 28 29 30 31 32 33 35 36 38 91];
% 
% 
