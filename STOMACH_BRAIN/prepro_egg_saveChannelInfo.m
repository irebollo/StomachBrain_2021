function prepro_egg_saveChannelInfo(cfgMain)

%{
Create a file in the folder scripts/files called EGG_peaks_info based on the log of
the preprocessing of each subject, it retrieves the mostpowerfull
frequency and channel of each subject and put it on a table. The original
selection of channels and frequency peaks is done on the first call to
prepro_egg, with the paramter cfgMainMain.automaticChannelSelection == 0

Input preprocessing log of each subject
Y:\Subjects\Subject13\Timeseries\EGGtimeseries\PhaseXvolume_S_13_fir2_fspread_015_ord_5_tw_15

Output
scripts/files/EGG_peaks_info

%}
subjects = global_subjectList

EGGpeaks = zeros(5,length(subjects))';

for iSubj = 1 : length(subjects)
    %load
    subj_idx = subjects(iSubj);
    dataDir = strcat(global_path2subject(subj_idx),'Timeseries',filesep,'data',filesep);
    filenameEGG = global_filename(subj_idx,cfgMain,'EGGPhaseXVolumeFilename');
    load(filenameEGG)
    EGGpeaks(iSubj,1)= subj_idx;
    EGGpeaks(iSubj,2)= logEGGpreprocessing.bestChannel;
    EGGpeaks(iSubj,3)= logEGGpreprocessing.mostPowerfullFrequency;
    EGGpeaks(iSubj,4)= logEGGpreprocessing.bestChannel == logEGGpreprocessing.mostPowerfullChannel; % to keep track which subjects the channel was chosen manually
%     EGGpeaks(iSubj,5)= logEGGpreprocessing.powerinChosenChannel; % to keep track which subjects the channel was chosen manually

end


rootFolder = global_path2root_folder();



save(strcat(rootFolder,'/scripts_4_github/files/','EGG_peaks_info'),'EGGpeaks')

% save(strcat('D:\NAVIGASTRIC\stomachBrainNavigastric\files\','EGG_peaks_info'),'EGGpeaks')
end