function [list_M, list_S, list_L] = global_roi_analysis(todo_AllSubjects, todo_ComputeMeanValues, todo_PrintSummaryFiles, todo_location, data_path, roi_path, output_dir, list_subjects, modality, list_acquisitions, list_sessions, list_roi, list_contrastnumbers, list_contrasts, contrasts_names, threshold, size_extent, EPIname, anat_dir, ANATname, modeldir, modelname, sym_name)

%----------------------------------------------------------------------
% INITIALIZATIONS
%----------------------------------------------------------------------

list_M = [];
list_S = [];
%number_of_minfs = 0; % To avoid taking into account the .minf files created by BrainVISA

% "list_subjects" must be completed if "todo_AllSubjects == 0"
if todo_AllSubjects
    clear list_subjects;
    list_images = dir(data_path);
    for i = 3:length(list_images) % To avoid "." and ".."
        if list_images(i).isdir == 1 % if it is a directory
            list_subjects{i-2} = list_images(i).name;
        end
    end
end

%----------------------------------------------------------------------
% COMPUTATION OF THE MEAN VALUES IN THE SELECTED ROIS 
% (computed from all the subjects)
%----------------------------------------------------------------------
if todo_ComputeMeanValues

    message = sprintf('\nComputing mean values...');
    disp(message);

    for t = 1:length(list_contrastnumbers)

        contrast = list_contrastnumbers(t);
        roi = fullfile(roi_path,list_roi{t});

        [M, S, L] = roi_analysis(data_path, list_subjects, modality, list_acquisitions, modeldir, contrast, roi, threshold);
        list_M(t)=M;
        list_S(t)=S;
        list_L(t).L=L;

    end

end

%----------------------------------------------------------------------
% PRINTING THE SUMMARY FILES WITH THE MEAN VALUES 
% AND THE VALUE FOR EACH SUBJECT
%----------------------------------------------------------------------

if todo_PrintSummaryFiles

    for k=1:length(list_subjects)

        cd(fullfile(data_path,list_subjects{k}));

        % We verify that the directory "fMRI" exists and go inside it if it does.
        list_insubject = dir(pwd);
        for ii = 3:length(list_insubject)
            dirinsubjectname = list_insubject(ii).name;

            if (strcmp(dirinsubjectname,modality))

                cd(dirinsubjectname);
                list_dir = dir(pwd);
                
                nb_acquisition = length(list_acquisitions);

                % FOR EACH ACQUISITION
                for mm = 1:nb_acquisition

                    cd(fullfile(data_path,list_subjects{k},modality,sprintf('acquisition%d',mm)));

                    % THIS IS DONE ONLY FOR ONE MODEL
                    
                    message = sprintf('\nCreating main contrasts page for subject %s, acquisition %s, model %s\n',list_subjects{k},list_acquisitions{mm},modelname);
                    disp(message);
                    output_name = sprintf('summary_%s_%s',list_subjects{k},modelname);

                    work_directory = fullfile(data_path,list_subjects{k},modality,list_acquisitions{mm},modeldir);

                    % The anatomical directory is fixed from the beginning in section "parameters"

                    con_plot(work_directory, threshold, size_extent, list_contrasts, contrasts_names, data_path, list_subjects{k}, modality, list_acquisitions{mm}, list_sessions, todo_location, output_dir, roi_path, list_roi, list_M, list_S, EPIname, anat_dir, ANATname, sym_name, modeldir, modelname, todo_ComputeMeanValues);

                    colormap(gray) % black and white

                    % Creation of a .tif file with the figures
                    fr = getframe(gcf);
                    [img, map] = frame2im(fr);
                    fi = fullfile(output_dir,strcat(output_name,'.tif'));
                    s = sprintf('\nSaving image in file %s',fi);
                    disp(s)
                    imwrite(img, fi);

                    % Creation of .ps and .pdf files with the figures
                    psfilename = fullfile(data_path,list_subjects{k},modality,list_acquisitions{mm},strcat(output_name,'.ps'));
                    pdffilename = fullfile(data_path,list_subjects{k},modality,list_acquisitions{mm},strcat(output_name,'.pdf'));
                    if exist(psfilename,'file')
                        system(sprintf('rm %s',psfilename));
                    end
                    if exist(pdffilename,'file')
                        system(sprintf('rm %s',pdffilename));
                    end
                    cd(fullfile(data_path,list_subjects{k},modality,list_acquisitions{mm}));
                    spm_print(output_name);
                    system(sprintf('ps2pdf %s',psfilename));

                    close

                end
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M, S, L] = roi_analysis(data_path, list_subjects, modality, list_acquisitions, modeldir, contrast, roi, threshold)

% Pour chaque contraste on aura des ROI différentes

% Output : M = mean
%          S = standard error
%          L = list of the number of voxels

roiheader = spm_vol(roi);
roimat = spm_read_vols(roiheader);
selection = find(roimat>0);

%----------------------------------------------------------------------
% MAIN PROGRAM (FOR EACH SUBJECT...)
%----------------------------------------------------------------------
j=1;
ListNbVoxel=[];

for k=1:length(list_subjects)

    cd(fullfile(data_path,list_subjects{k}));

    % We verify that the directory "fMRI" exists and go inside it if it does.
    list_insubject = dir(pwd);
    for ii = 3:length(list_insubject)
        dirinsubjectname = list_insubject(ii).name;

        if (strcmp(dirinsubjectname,modality))
        
            cd(dirinsubjectname);
            list_dir = dir(pwd);

            nb_acquisition = length(list_acquisitions);
            
            % FOR EACH ACQUISITION
            for mm = 1:nb_acquisition

                cd(fullfile(data_path,list_subjects{k},modality,sprintf('acquisition%d',mm),modeldir));

                spmTheader = spm_vol(sprintf('spmT_%04d.img',contrast));
                spmTmat = spm_read_vols(spmTheader);
                spmTselection = spmTmat(selection);

                nbvoxel = length(find(spmTselection>threshold));
                ListNbVoxel(j) = nbvoxel;
                j = j+1;

            end
        end
    end
end

M=mean(ListNbVoxel);
S=std(ListNbVoxel);
L=ListNbVoxel;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function con_plot(work_directory, threshold, size_extent, list_contrasts, contrasts_names, data_path, subject_name, modality, acquisition_name, list_sessions, todo_location, output_dir, roi_path, list_roi, list_M, list_S, EPIname, anat_dir, ANATname, sym_name, modeldir, modelname, todo_ComputeMeanValues)

%%%%% TITLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure_width_px = 840; % Square sheet
%figure_height_px = 840; % Square sheet
figure_width_px = 672; % A4 width = 21cm
figure_height_px = 950; % A4 height = 29.7cm
titlefontsize = 20;
titlefontweight = 'Bold';
subtitlefontsize = 10;

title_position_X = 0.1;
title_position_Y = 0.88; %0.87;
title_width =0.8;
title_height = 0.01;

title_position = [title_position_X title_position_Y title_width title_height];

number_of_contrasts = length(contrasts_names);
number_of_contrast_lines = 2;
%number_of_columns = ceil(number_of_contrasts/number_of_contrast_lines);
number_of_columns = 4; % FIXED FOR GENERAL BATCH

name_of_figure = sprintf('Summary for subject %s', subject_name);

h_figure = figure('Name',name_of_figure,'NumberTitle','off','Position',[0 0 figure_width_px figure_height_px]); % Dimensions of the figure window.
paper_width = 21; % cm - A4
paper_height = 29.7; % cm - A4
min_print_margin = 0.5; % cm
max_print_margin = 1; % cm
set(h_figure,'PaperPositionMode','manual','PaperUnits','centimeters',...
    'PaperPosition',[min_print_margin min_print_margin paper_width-max_print_margin paper_height-max_print_margin])
%'PaperPositionMode','auto'
%'PaperUnits','centimeters'
%'PaperPosition',[0 0 0.5 0.5]

ax=axes('Position',title_position,'Parent',h_figure);
set(ax,'Visible','off')
title_of_figure = sprintf('Summary for subject %s', subject_name);
string_title = {title_of_figure};
set(get(ax,'Title'),'String',string_title,'FontSize',titlefontsize,'FontWeight',titlefontweight,'Visible','on');
x = 0.5;
y = 0.5;
text(x,y,sprintf('in %s', fullfile(data_path,subject_name)),...
    'FontSize',subtitlefontsize,'HorizontalAlignment','center','Interpreter','none','Parent',ax);

% Added to have name of the anatomical image (containing scanner field and
% place of scanning)
complete_anatdir = fullfile(work_directory,'../../../..',anat_dir);
cd(complete_anatdir);
anatname = spm_select('List', complete_anatdir, ANATname);

y = -1;
text(x,y,sprintf('Anatomy acquisition : %s / Anatomy file : %s', anat_dir, anatname(1,:)),...
    'FontSize',subtitlefontsize,'FontWeight',titlefontweight,'HorizontalAlignment','center','Interpreter','none','Parent',ax);

%%%%% CONTRASTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contraststitle_fontsize = 10;
contraststitle_fontweight = 'Bold';

contrast_width = 0.9/number_of_columns;
contrast_height = (contrast_width*figure_width_px)/figure_height_px; % For having square images
margin_contrasts_X = 0.025;
margin_contrasts_Y = 1 - title_position_Y + 0.21; %0.22;
espace_between_contrasts_X = 0.05/(number_of_columns-1);
espace_between_contrasts_Y = 1/(number_of_contrast_lines+3);

score = 0;
sym_image_done = 0;

for n = 1:(length(list_contrasts)+1)

    line_index = floor((n-1)/number_of_columns); % First contrast line has index = 0
    n_in_line = n - line_index*number_of_columns;
    contrast_posX = margin_contrasts_X + (n_in_line-1)*(contrast_width+espace_between_contrasts_X);
    contrast_posY = 1 - margin_contrasts_Y - espace_between_contrasts_Y*line_index;
    contrast_position = [contrast_posX contrast_posY contrast_width contrast_height];

    ax=axes('Position',contrast_position,'Parent',h_figure);

    if (n == 4) % Place for the symmetry image

        sym_image_done = 1;

        %%%%% SYMMETRY IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        symtitle_fontsize = 10;
        symtitle_fontweight = 'Bold';
        sym_fontsize = 10;

        set(ax,'Visible','off')
        x = 0.2;
        y = 0.5;

        conname = fullfile(work_directory, sym_name);
        conheader = spm_vol(conname);
        condata = spm_read_vols(conheader);

        % NEW: for general image size
        dimx=conheader.dim(1);
        dimy=conheader.dim(2);
        dimz=conheader.dim(3);

        % NEW: for axial slice !!
        %slice = reshape(condata(27,:,:), 63, 46);
        %slice = reshape(condata(round(dimx/2),:,:), dimy, dimz);  % for a sagittal slice
        slice = reshape(condata(:,:,round(dimz*0.4)),dimx, dimy); % display slide around 0.4 * nb of slice

        slice_inv = zeros(size(slice));

        for i = 1:size(slice,1)
            slice_inv(size(slice,1)-(i-1),:) = slice(i,:);
        end

        slice = slice_inv';

        slice_inv = zeros(size(slice));
        for i = 1:size(slice,1)
            slice_inv(size(slice,1)-(i-1),:) = slice(i,:);
        end

        imagesc(slice_inv); axis image; axis tight; axis off;
        colormap Gray
        set(get(ax,'Title'),'String','Symmetry image','FontSize',symtitle_fontsize,'FontWeight',symtitle_fontweight,'Visible','on');

        %%%%% END OF SYMMETRY IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    else

        con_name = list_contrasts{n-sym_image_done};

        maskheader = spm_vol(fullfile(work_directory,'mask.img'));
        maskvol = spm_read_vols(maskheader);
        maskcol = maskvol(:)';
        m = find(maskcol>0);

        A = spm_vol(fullfile(work_directory,con_name));
        B = spm_read_vols(A);

        load(fullfile(work_directory,'SPM.mat'));
        vox = SPM.xVol.M;
        XYZ = SPM.xVol.XYZ;

        Tvalue=B(:)';

        X = XYZ(1,:)* A.mat(1,1) + A.mat(1,4);
        Y = XYZ(2,:)* A.mat(2,2) + A.mat(2,4);
        Z = XYZ(3,:)* A.mat(3,3) + A.mat(3,4);
        XYZ2 = [X; Y; Z];

        % Selection based on T-value
        Tvalue_masked = Tvalue(m);
        k = find(Tvalue_masked<threshold);
        Tvalue_masked(k) = 0;

        %spm_mip(Tvalue_masked,XYZ2,vox);
        MatProject = spm_mip_database(Tvalue_masked,XYZ2,vox);
        visu_improvement(MatProject);

        % Selection based on cluster size
        Tvalue_masked_thresholded = zeros(size(Tvalue));

        if todo_location % Generates the XYZ in the whole matrix
            clear XYZ
            k=1;
            for z=1:maskheader.dim(3)
                for y=1:maskheader.dim(2)
                    for x=1:maskheader.dim(1)
                        XYZ(:,k)=[x y z]';
                        k=k+1;
                    end
                end
            end
            save XYZ.mat XYZ
            system(sprintf('mv XYZ.mat %s',output_dir));
        else
            %M = load('XYZ.mat');
            [xyz_name,errmsg] = sprintf('%s/XYZ.mat',output_dir);
            M = load(xyz_name);
            XYZ = M.XYZ;
        end

        i=find(Tvalue>threshold);

        % NEW !!!
        % If there are some voxels with a value (T-value), we select them.
        % Otherwise, we do not select anything.
        if ~isempty(i)
            XYZ_select = XYZ(:,i);
            A = spm_clusters(XYZ_select);
        else
            A = [0];
        end

        Q     = [];
        for k = 1:max(A)
            j = find(A == k);
            if length(j) >= size_extent; Q = [Q j]; end
        end

        i_thr = i(Q);
        Tvalue_masked_thresholded(i_thr) = Tvalue(i_thr);

        Tvalue_masked_thresholded = Tvalue_masked_thresholded(m);

        %spm_mip(Tvalue_masked_thresholded,XYZ2,vox)
        MatProject = spm_mip_database(Tvalue_masked_thresholded,XYZ2,vox);
        visu_improvement(MatProject);

        title(contrasts_names{n-sym_image_done},'FontSize',contraststitle_fontsize,'FontWeight',contraststitle_fontweight);

        % For each contrast we show the mean and the standard deviation -------

        if ( ~isempty(list_M) & ~isempty(list_S) & (n-sym_image_done < length(list_M)+1) )

            % The number of voxels with a value higher than the threshold are
            % computed for each subject

            % Copied from function "roi_analysis" -------
            roi = fullfile(roi_path,list_roi{n-sym_image_done});

            roiheader = spm_vol(roi);
            roimat = spm_read_vols(roiheader);
            selection = find(roimat>0);

            cd(work_directory);

            spmTheader = spm_vol(con_name);
            spmTmat = spm_read_vols(spmTheader);
            spmTselection = spmTmat(selection);

            nb_of_voxels = length(find(spmTselection > threshold));
            % End copied from function "roi_analysis" ---

            %        % If the number of voxels above the threshold is further from the
            %        % mean value than 2 times the standard deviation, this means that
            %        % the subject is an outlier !!
            %        if ( (nb_of_voxels < list_M(n) - 2*list_S(n)) | (nb_of_voxels > list_M(n) + 2*list_S(n)) )
            % If the number of voxels above the threshold less than 10% of the mean value,
            % then we consider this subject as an outlier
            if ( nb_of_voxels < list_M(n-sym_image_done)/10 )
                conclusion_message = 'NOK';
                color_error = [1 0 0];
                conclusion_fontweight = 'Bold';
            else
                conclusion_message = 'OK';
                color_error = [0 0 0];
                conclusion_fontweight = 'Normal';
                score = score + 1;
            end

            % Title
            meantitle_fontsize = 8;
            meantitle_fontweight = 'Bold';
            x = 0.75;
            y = 0.375;
            text(x,y,sprintf('Voxels in ROI\n%d',nb_of_voxels),...
                'Units','normalized','FontSize',meantitle_fontsize,'FontWeight', meantitle_fontweight,...
                'HorizontalAlignment','center',...
                'Interpreter','none','Parent',ax);

            % Information
            mean_fontsize = 8;
            mean_fontweight = 'Normal';
            y = 0.21;
            text(x,y,sprintf('mean = %2.1f\nstd = %2.2f',...
                list_M(n-sym_image_done), list_S(n-sym_image_done)),...
                'Units','normalized','FontSize',mean_fontsize,'FontWeight', mean_fontweight,...
                'HorizontalAlignment','center',...
                'Interpreter','none','Parent',ax);

            % Conclusion
            y = 0.05;
            text(x,y,sprintf('Conclusion %s',conclusion_message),...
                'Units','normalized','FontSize',meantitle_fontsize,'FontWeight', conclusion_fontweight,...
                'HorizontalAlignment','center',...
                'Color', color_error,...
                'Interpreter','none','Parent',ax);

        end

    end

end

%%%%% CONTRAST IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (avant, on l'appelait "EPI image")

epititle_fontsize = 10;
epititle_fontweight = 'Bold';
epi_fontsize = 10;

n_in_line = number_of_contrasts + 1 - line_index*number_of_columns + 1;
epi_posX = margin_contrasts_X + (n_in_line-1)*(contrast_width+espace_between_contrasts_X);
correction_espace_between_contrasts_and_epi = 0.015;
epi_posY = 1 - margin_contrasts_Y - espace_between_contrasts_Y*line_index + correction_espace_between_contrasts_and_epi;
epi_position = [epi_posX epi_posY contrast_width contrast_height];

ax=axes('Position',epi_position,'Parent',h_figure);
set(ax,'Visible','off')
x = 0.2;
y = 0.5;

conname = fullfile(work_directory, EPIname);
conheader = spm_vol(conname);
condata = spm_read_vols(conheader);

% NEW: for general image size
dimx=conheader.dim(1);
dimy=conheader.dim(2);
dimz=conheader.dim(3);

%slice = reshape(condata(27,:,:), 63, 46);
slice = reshape(condata(round(dimx/2),:,:), dimy, dimz);  % for a sagittal slice

slice_inv = zeros(size(slice));

for i = 1:size(slice,1)
    slice_inv(size(slice,1)-(i-1),:) = slice(i,:);
end

slice = slice_inv';

slice_inv = zeros(size(slice));
for i = 1:size(slice,1)
    slice_inv(size(slice,1)-(i-1),:) = slice(i,:);
end

imagesc(slice_inv); axis image; axis tight; axis off;
colormap Gray
set(get(ax,'Title'),'String','Contrast image','FontSize',epititle_fontsize,'FontWeight',epititle_fontweight,'Visible','on');

%%%%% PROCESSING PARAMETERS AND MORE INFORMATION %%%%%%%%%%%%%%%%%%%%%%%%%%
summary_marge_X = 0.075; %0.1;
summary_width = 0.65; %0.8;
summary_height = 0.01;
summarytitle_fontsize = 10;
summarytitle_fontweight='Bold';
summary_fontsize = 10;
summary_color = [0 0 0];
score_color = [0 0 1];
summary_fontweight = 'Normal';
score_fontweight = 'Bold';
            
espace_between_contrasts_and_summary = 0.075; %0.090; %0.095;
summary_posY = epi_posY - espace_between_contrasts_and_summary;
summary_position = [summary_marge_X summary_posY summary_width summary_height];

ax=axes('Position',summary_position,'Parent',h_figure);
set(ax,'Visible','off')
set(get(ax,'Title'),'String',sprintf('Parameters and other information for subject %s:', subject_name),...
    'FontSize',summarytitle_fontsize,'FontWeight',summarytitle_fontweight,'Visible','on');
%y = 0.5;
y = -2;

if todo_ComputeMeanValues
    x = 0.2;
    text(x,y,sprintf('Subject = %s\nAcquisition = %s\nModel = %s',...
        subject_name, acquisition_name, modelname),...
        'HorizontalAlignment','center',...
        'FontSize',summary_fontsize,'FontWeight', summary_fontweight,...
        'HorizontalAlignment','center',...
        'Color', summary_color,...
        'Interpreter','none','Parent',ax);

    x = 0.5;
    text(x,y,sprintf('Threshold = %2.2f\nCluster size = %d',...
        threshold, size_extent),...
        'HorizontalAlignment','center',...
        'FontSize',summary_fontsize,'FontWeight', summary_fontweight,...
        'HorizontalAlignment','center',...
        'Color', summary_color,...
        'Interpreter','none','Parent',ax);

    x = 0.8;
    text(x,y,sprintf('Final score = %d/%d', score, length(list_contrasts)),...
        'HorizontalAlignment','center',...
        'FontSize',summary_fontsize,'FontWeight', score_fontweight,...
        'HorizontalAlignment','center',...
        'Color', score_color,...
        'Interpreter','none','Parent',ax);
else
    x = 0.3;
    text(x,y,sprintf('Subject = %s\nAcquisition = %s\nModel = %s',...
        subject_name, acquisition_name, modelname),...
        'HorizontalAlignment','center',...
        'FontSize',summary_fontsize,'FontWeight', summary_fontweight,...
        'HorizontalAlignment','center',...
        'Color', summary_color,...
        'Interpreter','none','Parent',ax);

    x = 0.7;
    text(x,y,sprintf('Threshold = %2.2f\nCluster size = %d',...
        threshold, size_extent),...
        'HorizontalAlignment','center',...
        'FontSize',summary_fontsize,'FontWeight', summary_fontweight,...
        'HorizontalAlignment','center',...
        'Color', summary_color,...
        'Interpreter','none','Parent',ax);
end

%%%%% ANATOMICAL IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
anattitle_fontsize = 10;
anattitle_fontweight = 'Bold';
anat_fontsize = 10;

anat_posX = epi_posX;
espace_between_contrasts_and_anat = 0.17;
anat_posY = epi_posY - espace_between_contrasts_and_anat;
anat_position = [anat_posX anat_posY contrast_width contrast_height];

ax=axes('Position',anat_position,'Parent',h_figure);
set(ax,'Visible','off')
x = 0.2;
y = 0.5;

complete_anatdir = fullfile(work_directory,'../../../..',anat_dir);
cd(complete_anatdir);
anatname = spm_select('List', complete_anatdir, ANATname);
anatheader = spm_vol(anatname(1,:));
anatdata = spm_read_vols(anatheader);

imagesize = size(anatdata);
slice = reshape(anatdata(floor(imagesize(1)/2 + 1),:,:), imagesize(2), imagesize(3));
slice_inv = zeros(size(slice));

for i = 1:size(slice,1)
    slice_inv(size(slice,1)-(i-1),:) = slice(i,:);
end

slice = slice_inv';

slice_inv = zeros(size(slice));
for i = 1:size(slice,1)
    slice_inv(size(slice,1)-(i-1),:) = slice(i,:);
end

imagesc(slice_inv); axis image; axis tight; axis off;
colormap Gray
set(get(ax,'Title'),'String','Anatomical image','FontSize',anattitle_fontsize,'FontWeight',anattitle_fontweight,'Visible','on');

%%%%% MOVEMENT CURVES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mwv_mat = [];
mvt_txt = '';
for nsess = 1:length(list_sessions)
    list_files = dir(fullfile(data_path,subject_name,modality,acquisition_name,list_sessions{nsess}));
    maxfiles = size(list_files,1);
    for j=3:maxfiles
        nom = list_files(j).name;
        [radical extent]=strtok(nom,'.');
        if strcmp(extent,'.txt')
            mvt_txt = list_files(j).name;
        end
    end

    if isempty(mvt_txt)
        return
    else
        swd = fullfile(data_path,subject_name,modality,acquisition_name,list_sessions{nsess},mvt_txt);
        mwv_temp = load(swd);
        mwv_mat = [mwv_mat; mwv_temp];
    end
end

x = mwv_mat(:,1);
y = mwv_mat(:,2);
z = mwv_mat(:,3);
pt = mwv_mat(:,4)*(180/pi); % pi radians = 180 degrees --> valeur * (180/pi)
rl = mwv_mat(:,5)*(180/pi);
ya = mwv_mat(:,6)*(180/pi);
datamat(i,1) = abs(max(x)-min(x));
datamat(i,2) = abs(max(y)-min(y));
datamat(i,3) = abs(max(z)-min(z));
datamat(i,4) = abs(max(pt)-min(pt));
datamat(i,5) = abs(max(rl)-min(rl));
datamat(i,6) = abs(max(ya)-min(ya));

bilan = sprintf('ampli x=%.2f, ampli y=%.2f, ampli z=%.2f, ampli pt=%.2f, ampli rl=%.2f, ampli ya=%.2f',datamat(i,1),datamat(i,2),datamat(i,3),datamat(i,4),datamat(i,5),datamat(i,6));

bilan2x = sprintf('max x=%.2f, min x=%.2f',max(x), min(x));
bilan2y = sprintf('max y=%.2f, min y=%.2f',max(y), min(y));
bilan2z = sprintf('max z=%.2f, min z=%.2f',max(z), min(z));
bilan2pt = sprintf('max pt=%.2f, min pt=%.2f',max(pt), min(pt));
bilan2rl = sprintf('max rl=%.2f, min rl=%.2f',max(rl), min(rl));
bilan2ya = sprintf('max ya=%.2f, min ya=%.2f',max(ya), min(ya));

bilan3x = sprintf('mean x=%.2f, std x=%.2f',mean(x), std(x));
bilan3y = sprintf('mean y=%.2f, std y=%.2f',mean(y), std(y));
bilan3z = sprintf('mean z=%.2f, std z=%.2f',mean(z), std(z));
bilan3pt = sprintf('mean pt=%.2f, std pt=%.2f',mean(pt), std(pt));
bilan3rl = sprintf('mean rl=%.2f, std rl=%.2f',mean(rl), std(rl));
bilan3ya = sprintf('mean ya=%.2f, std ya=%.2f',mean(ya), std(ya));

amplitude_x = max(x)-min(x);
amplitude_y = max(y)-min(y);
amplitude_z = max(z)-min(z);
amplitude_pt = max(pt)-min(pt);
amplitude_rl = max(rl)-min(rl);
amplitude_ya = max(ya)-min(ya);

% disp(bilan);
% disp(bilan3x);disp(bilan2x);disp(bilan3y);disp(bilan2y);disp(bilan3z);disp(bilan2z);
% disp(bilan3pt);disp(bilan2pt);disp(bilan3rl);disp(bilan2rl);disp(bilan3ya);disp(bilan2ya);

%---- Display curves ------------------------------------------------------

curves_height = 0.1;
curves_width = 0.65;
margin_curves_X = 0.075;
espace_between_summary_and_curves = 0.19; %0.145; % 0.15;
espace_between_curves = 0.17;

%---- Translation curve ---------------------------------------------------

curve1_posY = summary_posY - espace_between_summary_and_curves;
curve1_position = [margin_curves_X curve1_posY curves_width curves_height];

ax=axes('Position',curve1_position,'Parent',h_figure);

plot(ax,x)
hold on
plot(ax,y,'g')
plot(ax,z,'r')
s = ['x translation';'y translation';'z translation'];
legend(ax, s, 0) % gca = current axes handle
set(get(ax,'Title'),'String','translation','FontSize',10,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','images');
set(get(ax,'Ylabel'),'String','mm');
set(ax,'YLim',[-4 4])

%---- Summary values

summariestitle_fontweight = 'Bold';
summaries_fontsize = 10;

summaries_posX = margin_curves_X + curves_width + 0.05;
summaries_width = (1 - summaries_posX) - 0.05;
summaries_height = curves_height;

summaries_bias = 0.04;

summary1_posY = curve1_posY - summaries_bias;
summary1_position = [summaries_posX summary1_posY summaries_width summaries_height];

ax=axes('Position',summary1_position,'Parent',h_figure);
set(ax,'Visible','off')
set(get(ax,'Title'),'String','Translation amplitudes','FontSize',summaries_fontsize,'FontWeight',summariestitle_fontweight,'Visible','on');
x = 0.5;
y = 0.7;
text(x,y,sprintf('x amplitude = %.2f\ny amplitude = %.2f\nz amplitude = %.2f\n(in mm)',...
    amplitude_x, amplitude_y, amplitude_z),...
    'FontSize',summaries_fontsize,'HorizontalAlignment','center','Interpreter','none','Parent',ax);

%---- Rotation curve ------------------------------------------------------

curve2_posY = curve1_posY - espace_between_curves;
curve2_position = [margin_curves_X curve2_posY curves_width curves_height];

ax=axes('Position',curve2_position,'Parent',h_figure);

plot(ax,pt)
hold on
plot(ax,rl,'g')
plot(ax,ya,'r')
s = ['pitch';'roll ';'yaw  '];
legend(ax, s, 0) % gca = current axes handle
set(get(ax,'Title'),'String','rotation','FontSize',10,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','images');
set(get(ax,'Ylabel'),'String','degrees');
set(ax,'YLim',[-4 4])

%---- Summary values

summary2_posY = curve2_posY - summaries_bias;
summary2_position = [summaries_posX summary2_posY summaries_width summaries_height];

ax=axes('Position',summary2_position,'Parent',h_figure);
set(ax,'Visible','off')
set(get(ax,'Title'),'String','Rotation amplitudes','FontSize',summaries_fontsize,'FontWeight',summariestitle_fontweight,'Visible','on');
x = 0.5;
y = 0.7;
text(x,y,sprintf('pitch amplitude = %.2f\n roll amplitude = %.2f\nyaw amplitude = %.2f\n(in degrees)',...
    amplitude_pt, amplitude_rl, amplitude_ya),...
    'FontSize',summaries_fontsize,'HorizontalAlignment','center','Interpreter','none','Parent',ax);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MatrixProjection = spm_mip_database(Z,XYZ,M,units)

% By Ph. Pinel and A. Moreno - 21/08/2008
% spm_mip adapted function: in order to have as an output the projection
% matrix.

% SPM maximum intensity projection
% FORMAT spm_mip(Z,XYZ,M);
% Z       - vector point list of SPM values for MIP
% XYZ     - matrix of coordinates of points (Talairach coordinates)
% M       - voxels - > mm matrix or size of voxels (mm)
% units   - defining space     [default {'mm' 'mm' 'mm'}]
%         - Scalar specifies intensity of grid
%_______________________________________________________________________
%
% If the data are 2 dimensional [DIM(3) = 1] the projection is simply an
% image, otherwise:
%
% spm_mip creates and displays a maximum intensity projection of a point
% list of voxel values (Z) and their location (XYZ) in three orthogonal
% views of the brain.  It is assumed voxel locations conform to the space
% defined in the atlas of Talairach and Tournoux (1988); unless the third
% dimnesion is time.
%
% This routine loads a mip putline from MIP.mat. This is an image with
% contours and grids defining the space of Talairach & Tournoux (1988).
% mip05 corresponds to the Talairach atlas, mip96 to the MNI templates.
% The outline and grid are superimposed at intensity defaults.grid,
% defaulting to 0.4.
%
% A default colormap of 64 levels is assumed. The pointlist image is
% scaled to fit in the interval [0.25,1]*64 for display. Flat images
% are scaled to 1*64.
%
% If M or DIM are not specified, it is assumed the XYZ locations are
% in Talairach mm.
%
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston et al.
% $Id: spm_mip.m 1063 2008-01-04 12:12:54Z guillaume $

%-Get units and grid scaling
%--------------------------------------------------------------------------
global defaults
try, Grid = defaults.grid; catch, Grid = 0.4;               end
try, units;                catch, units = {'mm' 'mm' 'mm'}; end

% transpose locations if necessary
%--------------------------------------------------------------------------
if size(XYZ,1) ~= 3, XYZ = XYZ';         end
if size(Z,1)   ~= 1, Z   = Z';           end
if size(M,1)   == 1, M   = speye(4,4)*M; end

%-Scale & offset point list values to fit in [0.25,1]
%==========================================================================
Z    = Z - min(Z);
mx   = max(Z);
Scal = 8;
if isempty(mx),
    Z = [];
elseif isfinite(mx) && (numel(Z) ~= 1),
    Z = (1 + Scal*Z/mx)/(Scal + 1);
else
    Z = ones(1,length(Z));
end

%-Display format
%==========================================================================
load('MIP.mat');

%-Single slice case
%--------------------------------------------------------------------------
if isempty(units{3})

    %-2d case
    %----------------------------------------------------------------------
    mip = 4*grid_trans + mask_trans;
        
elseif units{3} == '%'
    
    %-3d case: Space-time
    %----------------------------------------------------------------------
    mip = 4*grid_time + mask_trans;

else
    %-3d case: Space
    %----------------------------------------------------------------------
    mip = 4*grid_all + mask_all;
end

% Load mip and create maximum intensity projection
%--------------------------------------------------------------------------
mip  = mip/max(mip(:));
c    = [0 0 0 ;
        0 0 1 ;
        0 1 0 ;
        0 1 1 ;
        1 0 0 ;
        1 0 1 ; 
        1 1 0 ; 
        1 1 1 ] - 0.5;
c    = c*M(1:3,1:3);
dim  = [(max(c) - min(c)) size(mip)];
d    = spm_project(Z,round(XYZ),dim);
mip  = max(d,Grid*mip);
%image(rot90((1 - mip)*64)); axis tight; axis off;
MatrixProjection = rot90((1 - mip)*64);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function visu_improvement(MatProject)

% Visualization improvement :
% In order to improve the display of the results, we put a transparent
% brain (white), a grey background (grey level = 50), and a white
% surrounding background.

for k = 1:size(MatProject,2)
    kk = MatProject(:,k);
    if ~isempty(find(kk<255))
        maxcol(k) = max(kk(find(kk<255)));
    else
        maxcol(k)=0;
    end
end
M = max(maxcol);
% We color the surrounding background in white.
MatProject(find(MatProject==M)) = 255;

for k = 1:size(MatProject,2)
    kk = MatProject(:,k);
    if ~isempty(find(kk<255))
        maxcol(k) = max(kk(find(kk<255)));
    else
        maxcol(k)=0;
    end
end
M = max(maxcol);
% We color the background in grey (level=50).
MatProject(find(MatProject==M)) = 50;

for k = 1:size(MatProject,2)
    kk = MatProject(:,k);
    if ~isempty(find(kk<255))
        maxcol(k) = max(kk(find(kk<255)));
    else
        maxcol(k)=0;
    end
end
M = max(maxcol);
% We color the not-active brain in white.
MatProject(find(MatProject==M)) = 255;

image(MatProject); axis tight; axis off;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
