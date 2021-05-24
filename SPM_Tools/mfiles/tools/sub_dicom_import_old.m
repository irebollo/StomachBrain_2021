function sub_dicom_import(rootdir,remove)

d = list_dirs(rootdir);

for i=1:length(d)
	sub_dicom_import(fullfile(rootdir,d{i}),remove);
end

f = spm_select('List',rootdir,'^*.MRDC');
if ~isempty(f)
	fprintf('DICOM conversion in %s\n',rootdir);
	jobs{1}.util{1}.dicom.data   = cellstr(f);
	jobs{1}.util{1}.dicom.outdir = cellstr(rootdir);
	jobs{1}.util{1}.dicom.convopts.format = 'img';
	spm_jobman('run',jobs);
	mkdir(rootdir,'dicom');
	system(sprintf('mv %s/*.MRDC %s',rootdir,fullfile(rootdir,'dicom')));
    if remove
        system(sprintf('rm -r %s',fullfile(rootdir,'dicom')));
    end
end

%=======================================================================
function d = list_dirs(rootdir)

d = dir(rootdir); 
d = {d([d.isdir]).name};
d = {d{cellfun('isempty',regexp(d,'^\.'))}};
