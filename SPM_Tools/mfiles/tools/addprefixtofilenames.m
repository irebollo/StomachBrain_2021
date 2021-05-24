function newlist = addprefixtofilenames(filelist,prefix)
% filelist must be a cell array where each element is a char-array
% containing a list of file names.

newlist=cell(1,length(filelist));

for i=1:length(filelist)
    wfiles=[];
    for img=1:size(filelist{i},1)
        [pth,nm,xt] = fileparts(deblank(filelist{i}(img,:)));
        wfiles  = strvcat(wfiles, fullfile(pth,[prefix nm xt]));
    end
    newlist{i}=wfiles;
end
