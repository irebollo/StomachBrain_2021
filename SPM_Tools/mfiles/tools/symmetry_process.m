%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% am - 01/08/2007
% am - update 05/03/2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function symmetry_process(data_path, subject, list_sessions, anat_path, orig_anat_name, contrasts_path)

load(fullfile(data_path,subject,contrasts_path,'SPM.mat'));

% Selection of the contrasts for which a T test is computed
ncons=[];
numberOfFtests = 0;
for i=1:size(SPM.xCon,2)
    if (SPM.xCon(i).STAT == 'T')    
        ncons(i-numberOfFtests) = i;
    elseif (SPM.xCon(i).STAT == 'F')
        numberOfFtests = numberOfFtests + 1; 
    end
end

contrasts_name_orig = 'con_';
contrasts_name = 'con_';

anat_name = sprintf('wm%s', char(orig_anat_name)); % We take the normalized image
gm_name = sprintf('wGM_%s', char(orig_anat_name));
wm_name = sprintf('wWM_%s', char(orig_anat_name));

sgm_name = sprintf('swGM_%s', char(orig_anat_name));
swm_name = sprintf('swWM_%s', char(orig_anat_name));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copy all the functional contrasts into an output directory (from "batch_copy.m")

disp('Copying all the functional contrasts into an output directory...');

for k=1:size(ncons,2)
    ncon=ncons(k);
    clear swd;
    contrastcompletename=[];
    swd = fullfile(data_path,subject,sprintf('%s%s%04d.img',contrasts_path,contrasts_name_orig,ncon));
    contrastcompletename=[contrastcompletename;swd];

    % Name and data of the images to be copied
    V=spm_vol(swd);
    D=spm_read_vols(V);

    % Name and path of the output (copied) images
    outDir=fullfile(data_path,subject,contrasts_path);
    name=fullfile(outDir,sprintf('%s%04d.img',contrasts_name,ncon));

    Vout=V;
    Vout.fname=name;

    Voutb=spm_write_vol(Vout, D);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flip all the anatomical images and write them into an output directory (from "batch_flipAnat.m")

disp('Flipping all the anatomical images and write them into an output directory...');

% FIRST FOR THE NORMALIZED ANATOMY, THEN FOR THE GREY MATTER AND THE WHITE
% MATTER, AND FINALLY FOR THE SMOOTHED GREY AND WHITE MATTERS

if ~isempty(spm_select('List',fullfile(data_path,subject,anat_path),'^wGM.*\.img')) & ~isempty(spm_select('List',fullfile(data_path,subject,anat_path),'^wWM.*\.img')) % in case of grey / white matter already processed !!

    for gmwmcount = 1:5

        swd = [];
        if (gmwmcount == 1)
            name = fullfile(data_path,subject,anat_path,sprintf('%s.img',anat_name));
        elseif (gmwmcount == 2)
            name = fullfile(data_path,subject,anat_path,sprintf('%s.img',gm_name));
        elseif (gmwmcount == 3)
            name = fullfile(data_path,subject,anat_path,sprintf('%s.img',wm_name));
        elseif (gmwmcount == 4)
            name = fullfile(data_path,subject,anat_path,sprintf('%s.img',sgm_name));
        elseif (gmwmcount == 5)
            name = fullfile(data_path,subject,anat_path,sprintf('%s.img',swm_name));
        end
        swd=[swd;name];

        flip_name = swd ;

        for i = 1:size(flip_name,1)

            if ( isempty(deblank(flip_name(i,:))) ~= 1 )

                [pathDir,pathName,niu,niu] = fileparts(deblank(flip_name(i,:)));
                outDir=fullfile(data_path,subject,anat_path);
                flip_img = fullfile(outDir, strcat('flip_',pathName,'.img'));

                V_P = spm_vol(deblank(flip_name(i,:)));
                dim_x = V_P.dim(1);
                dim_y = V_P.dim(2);
                dim_z = V_P.dim(3);
                M = V_P.mat;

                % Initialization of matrix "mat_image" (reading slice by slice the initial image)
                mat_image = ones(dim_x,dim_y);
                % Initialization of the resulting flipped matrix
                flip_image = ones(dim_x,dim_y,dim_z) ;
                % Image reading parameter (cf. spm_slice_vol.m)
                hold = 0 ;

                % Reading the image slice by slice (in z) and flipping
                for k = 1:1:dim_z
                    % Reading slice z=k
                    Ma = spm_matrix( [0 0 k] );
                    mat_image(:,:) = spm_slice_vol( V_P, Ma, V_P.dim(1:2), hold );
                    % Flip in x
                    for i = 1:1:dim_x
                        flip_image(i,:,k) = mat_image(dim_x+1-i,:);
                    end
                end

                % Writing the resulting flipped matrix
                descrip = 'Flipped image';
                Vo = struct('fname', flip_img,...
                    'dim', [dim_x, dim_y, dim_z],...
                    'dt', V_P.dt,...
                    'mat', M,...
                    'descrip', descrip);
                Vo = spm_write_vol(Vo,flip_image);

                clear('V_P','Vo');

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flip all the functional contrast images and write them into an output
% directory (from "batch_flipCon.m")

disp('Flipping all the functional contrast images and write them into an output directory...');

% For each functional image (contrast)
for k=1:size(ncons,2)   % CHANGER ORDRE DES BOUCLES
    ncon=ncons(k);

    swd = [];
    name = fullfile(data_path,subject,contrasts_path,sprintf('%s%04d.img',contrasts_name,ncon));
    swd=[swd;name];

    flip_name = swd;

    for i = 1:size(flip_name,1)

        if ( isempty(deblank(flip_name(i,:))) ~= 1 )

            [pathDir,pathName,niu,niu] = fileparts(deblank(flip_name(i,:)));
            outDir=fullfile(data_path,subject,contrasts_path);
            flip_img = fullfile(outDir,strcat('flip_',pathName,'.img'));

            V_P = spm_vol(deblank(flip_name(i,:)));
            dim_x = V_P.dim(1);
            dim_y = V_P.dim(2);
            dim_z = V_P.dim(3);
            M = V_P.mat;

            % Initialization of matrix "mat_image" (reading slice by slice the initial image)
            mat_image = ones(dim_x,dim_y);
            % Initialization of the resulting flipped matrix
            flip_image = ones(dim_x,dim_y,dim_z);
            % Image reading parameter (cf. spm_slice_vol.m)
            hold = 0;

            % Reading the image slice by slice (in z) and flipping
            for k = 1:1:dim_z
                % Reading slice z=k
                Ma = spm_matrix( [0 0 k] );
                mat_image(:,:) = spm_slice_vol( V_P, Ma, V_P.dim(1:2), hold );
                % Flip in x
                for i = 1:1:dim_x
                    flip_image(i,:,k) = mat_image(dim_x+1-i,:);
                end
            end

            % Reading the resulting flipped matrix
            descrip = 'Flipped image';
            Vo = struct('fname', flip_img,...
                'dim', [dim_x, dim_y, dim_z],...
                'dt', V_P.dt,...
                'mat', M,...
                'descrip', descrip);
            Vo = spm_write_vol(Vo,flip_image);

            clear('V_P','Vo');

        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalize individually each flipped anatomy towards the original anatomy and
% apply individually each resulting normalization matrix to each flipped contrast image
% (from "batch_normalizedflipped.m" and "batch_ApplyIndividualNormalization.m")

disp('Normalizing individually each flipped anatomy towards the original anatomy and applying individually each resulting normalization matrix to each flipped contrast image...');

clear swd;
original_anat = fullfile(data_path,subject,anat_path,sprintf('%s.img',anat_name));
flipped_anat = fullfile(data_path,subject,anat_path,sprintf('flip_%s.img',anat_name));

disp('Processing subject:');
subject

% Normalize "flipped_anat" onto "original_anat"

jobs{1}.spatial{1}.normalise{1}.estwrite.ref = cellstr(flipped_anat);
jobs{1}.spatial{1}.normalise{1}.estwrite.source = cellstr(original_anat);

% Normalization of constrats

for k=1:size(ncons,2)
    ncon=ncons(k);
    clear original_flippedimage;
    original_flippedimage = fullfile(data_path,subject,contrasts_path,sprintf('flip_%s%04d.img',contrasts_name,ncon));
    jobs{1}.spatial{1}.normalise{1}.estwrite.other(k) = cellstr(original_flippedimage);
end

if ~isempty(spm_select('List',fullfile(data_path,subject,anat_path),'^wGM.*\.img')) & ~isempty(spm_select('List',fullfile(data_path,subject,anat_path),'^wWM.*\.img')) % in case of grey / white matter already processed !!
    % Normalization of the flipped grey and white matters
    kmax = size(ncons,2);
    flipped_gm = fullfile(data_path,subject,anat_path,sprintf('flip_%s.img',gm_name));
    flipped_wm = fullfile(data_path,subject,anat_path,sprintf('flip_%s.img',wm_name));
    jobs{1}.spatial{1}.normalise{1}.estwrite.other(kmax+1) = cellstr(flipped_gm);
    jobs{1}.spatial{1}.normalise{1}.estwrite.other(kmax+2) = cellstr(flipped_wm);

    % Normalization of the flipped smoothed grey and white matters
    flipped_sgm = fullfile(data_path,subject,anat_path,sprintf('flip_%s.img',sgm_name));
    flipped_swm = fullfile(data_path,subject,anat_path,sprintf('flip_%s.img',swm_name));
    jobs{1}.spatial{1}.normalise{1}.estwrite.other(kmax+3) = cellstr(flipped_sgm);
    jobs{1}.spatial{1}.normalise{1}.estwrite.other(kmax+4) = cellstr(flipped_swm);
end

spm_jobman('run',jobs{1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply individually the substraction "scon - nflipped scon" (original contrast minus
% normalized flipped contrast (from "batch_symmetry.m")

disp('Apply individually the substraction "scon - nflipped scon" (original contrast minus normalized flipped contrast...');

    for k=1:size(ncons,2)
        ncon=ncons(k);
        clear inputnames;
        functionstr = 'i1 - i2';
        inputnames{1} = fullfile(data_path,subject,contrasts_path,sprintf('%s%04d.img',contrasts_name,ncon));
        inputnames{2} = fullfile(data_path,subject,contrasts_path,sprintf('flip_%s%04d.img',contrasts_name,ncon));
        outname = fullfile(data_path,subject,contrasts_path,sprintf('diff_%s%04d.img',contrasts_name,ncon));
        spm_imcalc_ui(inputnames,outname,functionstr);
    end

if ~isempty(spm_select('List',fullfile(data_path,subject,anat_path),'^wGM.*\.img')) & ~isempty(spm_select('List',fullfile(data_path,subject,anat_path),'^wWM.*\.img')) % in case of grey / white matter already processed !!

    % Grey matter
    clear inputnames;
    inputnames{1} = fullfile(data_path,subject,anat_path,sprintf('%s.img',gm_name));
    inputnames{2} = fullfile(data_path,subject,anat_path,sprintf('flip_%s.img',gm_name));
    outname = fullfile(data_path,subject,anat_path,sprintf('diff_%s.img',gm_name));
    spm_imcalc_ui(inputnames,outname,functionstr);

    % White matter
    clear inputnames;
    inputnames{1} = fullfile(data_path,subject,anat_path,sprintf('%s.img',wm_name));
    inputnames{2} = fullfile(data_path,subject,anat_path,sprintf('flip_%s.img',wm_name));
    outname = fullfile(data_path,subject,anat_path,sprintf('diff_%s.img',wm_name));
    spm_imcalc_ui(inputnames,outname,functionstr);

    % Smoothed grey matter
    clear inputnames;
    inputnames{1} = fullfile(data_path,subject,anat_path,sprintf('%s.img',sgm_name));
    inputnames{2} = fullfile(data_path,subject,anat_path,sprintf('flip_%s.img',sgm_name));
    outname = fullfile(data_path,subject,anat_path,sprintf('diff_%s.img',sgm_name));
    spm_imcalc_ui(inputnames,outname,functionstr);

    % Smoothed white matter
    clear inputnames;
    inputnames{1} = fullfile(data_path,subject,anat_path,sprintf('%s.img',swm_name));
    inputnames{2} = fullfile(data_path,subject,anat_path,sprintf('flip_%s.img',swm_name));
    outname = fullfile(data_path,subject,anat_path,sprintf('diff_%s.img',swm_name));
    spm_imcalc_ui(inputnames,outname,functionstr);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Smoothing of the difference images of the contrasts

jobs_dsmooth = {};
nbjobs_dsmooth = 0;
logfile_dsmooth = fullfile(data_path,subject,'dsmooth.log');

diff_smoothing = 5;
%diffwc  = '^diff.*\.img$';
diffwc  = '^diff_con.*\.img$';

logmsg(logfile_dsmooth,'Scanning for difference files...');
diff_images = cellstr(spm_select('List', fullfile(data_path,subject,contrasts_path), diffwc));

for n=1:length(diff_images)
    diff_images{n} = fullfile(data_path,subject,contrasts_path,diff_images{n});
end

logmsg(logfile_dsmooth,sprintf('Smoothing %d files ("%s"...) with fwhm = %d mm',sum(cellfun('size',diff_images,1)),diff_images{1}(1,:),diff_smoothing));
nbjobs_dsmooth = nbjobs_dsmooth + 1;
jobs_dsmooth{nbjobs_dsmooth}.spatial{1}.smooth.data = cellstr(strvcat(diff_images));
jobs_dsmooth{nbjobs_dsmooth}.spatial{1}.smooth.fwhm = diff_smoothing;

logmsg(logfile_dsmooth,sprintf('Job batch file saved in %s.',fullfile(data_path,subject,'jobs_diff_smooth.mat')));
save(fullfile(data_path,subject,'jobs_diff_smooth.mat'),'jobs_dsmooth');
spm_jobman('run',jobs_dsmooth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Smoothing of the difference images of the grey and white matters

if ~isempty(spm_select('List',fullfile(data_path,subject,anat_path),'^wGM.*\.img')) & ~isempty(spm_select('List',fullfile(data_path,subject,anat_path),'^wWM.*\.img')) % in case of grey / white matter already processed !!

    jobs_dsmooth = {};
    nbjobs_dsmooth = 0;
    logfile_dsmooth = fullfile(data_path,subject,'dgwmsmooth.log');

    diff_smoothing = 3;
    diffwc  = '^diff.*\.img$';

    logmsg(logfile_dsmooth,'Scanning for difference files...');
    diff_images = cellstr(spm_select('List', fullfile(data_path,subject,anat_path), diffwc));

    for n=1:length(diff_images)
        diff_images{n} = fullfile(data_path,subject,anat_path,diff_images{n});
    end

    logmsg(logfile_dsmooth,sprintf('Smoothing %d files ("%s"...) with fwhm = %d mm',sum(cellfun('size',diff_images,1)),diff_images{1}(1,:),diff_smoothing));
    nbjobs_dsmooth = nbjobs_dsmooth + 1;
    jobs_dsmooth{nbjobs_dsmooth}.spatial{1}.smooth.data = cellstr(strvcat(diff_images));
    jobs_dsmooth{nbjobs_dsmooth}.spatial{1}.smooth.fwhm = diff_smoothing;

    logmsg(logfile_dsmooth,sprintf('Job batch file saved in %s.',fullfile(data_path,subject,'jobs_diffgwm_smooth.mat')));
    save(fullfile(data_path,subject,'jobs_diffgwm_smooth.mat'),'jobs_dsmooth');
    spm_jobman('run',jobs_dsmooth);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creation of directory for flipped images and difference images

complete_path = fullfile(data_path,subject,contrasts_path);
cd(complete_path);

mkdir('symmetry_images');

% Move the flipped and the difference images to the "symmetry_images" directory

ligne=sprintf('!mv flip_* symmetry_images');
eval(ligne)
ligne=sprintf('!mv diff_* symmetry_images');
eval(ligne)
ligne=sprintf('!mv sdiff_* symmetry_images');
eval(ligne)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Symmetry batch finished!');
