function rename_files(data_path, subject, anat_path, anat_name, session_path, study_name, origin_name)


switch origin_name

    case '15T_shfj' % s***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % RENAME THE ANATOMICAL IMAGES

        cd(fullfile(data_path,subject,anat_path));
        list_images=dir(fullfile(data_path,subject,anat_path));

        % Rename all old-named files (of kind "anat_*****" or "anat-*****")
        for i=1:length(list_images)
            [nom, extension]=strtok(list_images(i).name,'.');
            if  (strcmp(extension,'.img') | strcmp(extension,'.hdr') | ...
                 strcmp(extension,'.ima') | strcmp(extension,'.dim') | strcmp(extension,'.mat'))
                image=list_images(i).name;
                typeimage=nom(1);
                if typeimage=='w'
                    prefix='w' ;
                elseif typeimage=='n'
                    prefix='nobias_';
                else prefix='' ;
                end
                image1=sprintf('%s%s%s',prefix,anat_name,extension);
                if ~strcmp(image,image1)
                    ligne=sprintf('!mv %s %s',image,image1);
                    eval(ligne);
                end
            end
            if  strcmp(extension,'.ima.minf')
                image=list_images(i).name;
                if typeimage=='w'
                    prefix='w' ;
                elseif typeimage=='n'
                    prefix='nobias_';
                else prefix='' ;
                end
                image1=sprintf('%s%s%s',prefix,anat_name,extension);
                if ~strcmp(image,image1)
                    ligne=sprintf('!mv %s %s',image,image1);
                    eval(ligne);
                end
            end
        end
        
        % RENAME THE FUNCTIONAL IMAGES

        cd(fullfile(data_path,subject,session_path));
        list_images=dir(fullfile(data_path,subject,session_path));

        % Erase files of non interest
        ligne='!rm *aim*.*';
        eval(ligne);
        ligne='!rm *.mat';
        eval(ligne);

        % Rename all old-named files (containing "localizer_corr****")
        for i=1:length(list_images)
            [nom, extension]=strtok(list_images(i).name,'.');
            if  strcmp(extension,'.img') | strcmp(extension,'.hdr')
                image=list_images(i).name;
                len2=length(image);
                len1=len2-8;
                typeimage=nom(1);
                if typeimage=='s'
                    prefix='swa';
                elseif typeimage=='w'
                    prefix='wa' ;
                elseif typeimage=='a'
                    prefix='a';
                else prefix='' ;
                end
                image1=sprintf('%s%s_%s',prefix,study_name,image(len1+1:len2));
                ligne=sprintf('!mv %s %s',image,image1);
                eval(ligne);
            end
        end


    case '3T_shfj' % bru**** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % RENAME THE ANATOMICAL IMAGES

        cd(fullfile(data_path,subject,anat_path));
        list_images=dir(fullfile(data_path,subject,anat_path));

        % Rename all old-named files (containing only the old NIP)
        for i=1:length(list_images)
            [nom, extension]=strtok(list_images(i).name,'.');
            if  (strcmp(extension,'.img') | strcmp(extension,'.hdr') | ...
                 strcmp(extension,'.ima') | strcmp(extension,'.dim'))
                image=list_images(i).name;
                len2=length(image);
                len1=len2 - (length(subject)+4);
                image1=sprintf('%s%s%s',image(1:len1),anat_name,extension);
                if ~strcmp(image,image1)
                    ligne=sprintf('!mv %s %s',image,image1);
                    eval(ligne);
                end
            end
            if  strcmp(extension,'.ima.minf')
                image=list_images(i).name;
                len2=length(image);
                len1=len2 - (length(subject)+9);
                image1=sprintf('%s%s%s',image(1:len1),anat_name,extension);
                if ~strcmp(image,image1)
                    ligne=sprintf('!mv %s %s',image,image1);
                    eval(ligne);
                end
            end
            if  strcmp(extension,'.mat')
                image=list_images(i).name;
                len2=length(image);
                len1=length(subject)+1;
                image1=sprintf('%s%s',anat_name,image(len1:len2));
                if ~strcmp(image,image1)
                    ligne=sprintf('!mv %s %s',image,image1);
                    eval(ligne);
                end
            end
        end

        % RENAME THE FUNCTIONAL IMAGES

        cd(fullfile(data_path,subject,session_path));
        list_images=dir(fullfile(data_path,subject,session_path));

        % Erase files of non interest
        ligne='!rm *aim*.*';
        eval(ligne);
        ligne='!rm *.mat';
        eval(ligne);

        % Rename all old-named files (containing "ima1")
        for i=1:length(list_images)
            [nom, extension]=strtok(list_images(i).name,'.');
            if  strcmp(extension,'.img') | strcmp(extension,'.hdr')
                image=list_images(i).name;
                len2=length(image);
                len1=len2-8;
                typeimage=nom(1);
                if typeimage=='s'
                    prefix='swa';
                elseif typeimage=='w'
                    prefix='wa' ;
                elseif typeimage=='a'
                    prefix='a';
                else prefix='' ;
                end
                image1=sprintf('%s%s_%s',prefix,study_name,image(len1+1:len2));
                ligne=sprintf('!mv %s %s',image,image1);
                eval(ligne);
            end
        end

        
    case '3T_neurospin' % XX****** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % RENAME THE ANATOMICAL IMAGES

        cd(fullfile(data_path,subject,anat_path));
        list_images=dir(fullfile(data_path,subject,anat_path));

        % Rename files recovered with NmrServer
        for i=1:length(list_images)
            [nom, extension]=strtok(list_images(i).name,'.');
            if  (strcmp(extension,'.img') | strcmp(extension,'.hdr') | ...
                 strcmp(extension,'.ima') | strcmp(extension,'.dim'))
                image=list_images(i).name;
                typeimage=nom(1);
                if typeimage=='c'
                    if nom(2)=='1'
                        prefix='c1';
                    elseif nom(2)=='2'
                        prefix='c2';
                    end
                elseif typeimage=='w'
                    prefix='wm' ;
                elseif typeimage=='m'
                    prefix='m';
                else prefix='' ; % This also takes into account images called anat.*
                end
                image1=sprintf('%s%s%s',prefix,anat_name,extension);
                if ~strcmp(image,image1)
                    ligne=sprintf('!mv %s %s',image,image1);
                    eval(ligne);
                end
            end
            if  strcmp(extension,'.ima.minf')
                image=list_images(i).name;
                if typeimage=='c'
                    if nom(2)=='1'
                        prefix='c1';
                    elseif nom(2)=='2'
                        prefix='c2';
                    end
                elseif typeimage=='w'
                    prefix='wm' ;
                elseif typeimage=='m'
                    prefix='m';
                else prefix='' ;
                end
                image1=sprintf('%s%s%s',prefix,anat_name,extension);
                if ~strcmp(image,image1)
                    ligne=sprintf('!mv %s %s',image,image1);
                    eval(ligne);
                end
            end
            if  strcmp(extension,'.mat')
                image=list_images(i).name;
                len2=length(image);
                if nom(length(image)-7)=='v' % 8 = 4 (extension) + 3
                    len1=len2 - (4+10); % 4 - extension; 10 - "_seg_inv_sn"
                    image1=sprintf('%s%s',anat_name,image(len1:len2));
                    if ~strcmp(image,image1)
                        ligne=sprintf('!mv %s %s',image,image1);
                        eval(ligne);
                    end
                elseif nom(length(image)-7)=='g' % 8 = 4 (extension) + 3
                    len1=len2 - (4+6); % 4 - extension; 6 - "_seg_sn"
                    image1=sprintf('%s%s',anat_name,image(len1:len2));
                    if ~strcmp(image,image1)
                        ligne=sprintf('!mv %s %s',image,image1);
                        eval(ligne);
                    end
                else
                end
            end
        end

        % RENAME THE FUNCTIONAL IMAGES

        cd(fullfile(data_path,subject,session_path));
        list_images=dir(fullfile(data_path,subject,session_path));

        % Erase files of non interest
        %ligne='!rm *aim*.*';
        %eval(ligne);
        %ligne='!rm *.mat';
        %eval(ligne);

        % Rename files recovered with NmrServer
        for i=1:length(list_images)
            [nom, extension]=strtok(list_images(i).name,'.');
            if  (strcmp(extension,'.img') | strcmp(extension,'.hdr') | ...
                 strcmp(extension,'.txt'))
                typeimage=nom(1);
                if typeimage=='s'
                    prefix='swa';
                elseif typeimage=='w'
                    prefix='wa' ;
                elseif typeimage=='a'
                    prefix='a';
                elseif typeimage=='m'
                    prefix='meana';
                elseif typeimage=='r'
                    prefix='rp_a';
                else prefix='' ;
                end
                image=list_images(i).name;
                len2=length(image);
%                if nom(length(image)-7)=='-' % 7 = 4 (extension) + 3 ("-01")
                if nom(length(image)-6)=='-' % 7 = 4 (extension) + 3 ("-01") - 1
                    len1=len2 - (4+3); % 4 - extension; 3 - "-01"
                    image1=sprintf('%s%s_%s%s',prefix,study_name,image(len1-5:len1),extension);
                else % if the name of the images is localizer_****
                    len1=len2 - 4; % 4 - extension
                    image1=sprintf('%s%s_00%s%s',prefix,study_name,image(len1-3:len1),extension);
                end
                ligne=sprintf('!mv %s %s',image,image1);
                eval(ligne);
            end
        end
end

%otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nothing

cd(data_path);

end
