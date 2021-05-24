function d = subdirectories(wd, filt)
%Extract subdirectories of a given directory

error(nargchk(1,2,nargin));
if nargin == 1, filt = '.*'; end

[dummy, d] = spm_select('List',wd,filt);
d(ismember(d,{'.' '..'}),:) = '';
d(cellfun('isempty',regexp(cellstr(d),filt)),:) = '';