function copy_anat_disque_dur(confirm)

%----------------------------------------------------------------------
% GENERAL PARAMETERS
%----------------------------------------------------------------------

addpath('/neurospin/unicog/protocols/IRMf/MainDatabaseLocalizers_PinelMoreno_2008/Tools/ContrastsSummary/mfiles');

data_path = '/neurospin/unicog/protocols/IRMf/MainDatabaseLocalizers_PinelMoreno_2008/Subjects';

output_dir = '/media/IOMEGA HDD/MainDatabaseLocalizers_PinelMoreno_2008/Subjects';
output_dir_string = '/media/IOMEGA\ HDD/MainDatabaseLocalizers_PinelMoreno_2008/Subjects';

%----------------------------------------------------------------------
% INITIALIZATIONS
%----------------------------------------------------------------------

% All the subjects are selected
list_images = dir(data_path);
for i = 3:length(list_images) % To avoid "." and ".."
    list_subjects{i-2}=list_images(i).name;
end

anat_dir = 't1mri';

%----------------------------------------------------------------------
% COPY OF ALL THE ORIGINAL ANATOMIES IN THE EXTERNAL HARD DISK
%----------------------------------------------------------------------

for k=1:length(list_subjects)

    mkdir(output_dir,list_subjects{k});
    mkdir(fullfile(output_dir,list_subjects{k}),anat_dir);
    
    % The number of acquisitions is detected automatically for each subject
    cd(fullfile(data_path,list_subjects{k},anat_dir));
    list_dir = dir(pwd);
    nb_acquisition = 0;
    for i = 3:length(list_dir)
        dirname = list_dir(i).name;
        if strcmp(dirname(1:4),'acqu')
            nb_acquisition = nb_acquisition+1;
        end
    end
    
    % FOR EACH ACQUISITION
    for mm = 1:nb_acquisition

        cd(fullfile(data_path,list_subjects{k},anat_dir,sprintf('acquisition%d',mm)));
        
        mkdir(fullfile(output_dir,list_subjects{k},anat_dir),sprintf('acquisition%d',mm));
%        system(sprintf('cp anat_* %s',fullfile(output_dir_string,list_subjects{k},anat_dir,sprintf('acquisition%d',mm))));
        system(sprintf('rm %s/*.mat',fullfile(output_dir_string,list_subjects{k},anat_dir,sprintf('acquisition%d',mm))));

    end

end



%end

