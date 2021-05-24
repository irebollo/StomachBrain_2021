

rootDir= 'C:\PHYSIENS\Physiens\Subjects\';


subjectlist= [9 10 13 15 16 17 18 19 20 21 22 23 24 25 26 28 29 30 31 32 33];

for iSubject = 1 : length(subjectlist)
    
    subj_idx = subjectlist(iSubject)
    
    
subjectDir = strcat(rootDir,'Subject',sprintf('%.2d',subj_idx)) ; 
inputDir = strcat(subjectDir,'\fMRI\acquisition1\RestingState');
outputDir = strcat(subjectDir,'\fMRI\acquisition1\RestingStateSmooth3');
%mkdir(outputDir)
%copyfile(strcat(inputDir,'\f*'),outputDir)

% found a more effiecient way to save files in same folder, will delete
% files now
rmdir(outputDir,'s')

end