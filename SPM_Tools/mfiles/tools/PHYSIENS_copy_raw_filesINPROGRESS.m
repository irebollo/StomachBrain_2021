function    PHYSIENS_copy_raw_files(subj_idx)

rawFolder = strcat('G:\PHYSIENS\',sprintf('%.2d',subj_idx));
targetFolder = strcat('C:\PHYSIENS\Physiens\Subjects\Subject',sprintf('%.2d',subj_idx),'\');

brainampDir = strcat(subjectFolder,'\brainamp\with MRI\');
    files= dir( fullfile( brainampDir,'*.vhdr')); %# list all *.vhdr files
    filename = {files.name}';%'# file names
    output = char(strcat(brainampDir,filename));

end