function path2root = global_path2root(sample)
% path to Physiens root folder

% path2root = ['D:\NAVIGASTRIC\Navigastric\'];
% path2root = ['D:\NAVIGASTRIC\test2pipelines\'];
% path2root = ['Z:\dataPhysiens\'];

rootFolder = global_path2root_folder();

% rootFolder = '/media/Navigastric/';
switch sample
    case 1
        path2root = [rootFolder 'subjects' filesep];
    case 2
        path2root = [rootFolder 'dataPhysiens' filesep];

end

end