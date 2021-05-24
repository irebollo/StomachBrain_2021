function NS_realign_unwarp(fm_dir,epi_dir)
%----------------------------------------------------------------------
% Load epi data from data directory for each session if necessary
%----------------------------------------------------------------------

if iscell(epi_dir)
    nsessions = size(epi_dir,2);
    for sessnum=1:nsessions
        epi_all = spm_select('List',epi_dir{sessnum},'^.*\.img$');
        nvol = size(epi_all,1);
        for ivol=1:nvol
            epi_img{sessnum, ivol}=fullfile(epi_dir{sessnum},epi_all(ivol,:));
        end
    end
elseif isstr(epi_dir)
    nsessions=1;
    sessnum=1;
    epi_all = spm_select('List',epi_dir,'^.*\.img$');
    nvol = size(epi_all,1);
    for ivol=1:nvol
        epi_img{sessnum, ivol}=fullfile(epi_dir,epi_all(ivol,:));
    end
end

opts_unwarp.M = spm_get_space(epi_img{1,1});

for sessnum=1:nsessions
    nvol = size(epi_img, 2);
    EndOfVDMName = sprintf('session%d.img', sessnum);
    % Corrected by A. Moreno - 2009-12-02: the vdm5_ files do not have the
    % name of the session at the end of their name when there is only 1 session !!
    if nsessions==1
        VDMName = spm_select('List',fm_dir, 'vdm5_.*\.img$');
    else
        VDMName = spm_select('List',fm_dir, ['vdm5_.*' EndOfVDMName]);
    end
    opts_unwarp.sfP = fullfile(fm_dir, VDMName);
    for ivol=2:nvol
        if ~isempty(epi_img{sessnum, ivol})
            s = spm_uw_estimate(epi_img{sessnum, ivol},opts_unwarp);
            spm_uw_apply(s);
        end
    end
end
