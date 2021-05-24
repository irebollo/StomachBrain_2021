
%{
This script will 

1- resample the gastric network to the same voxel resolution as Margulies gradients
2- Load the two gradients from Margulies et al 2017 PNAS, and the gastric network in the same space 
3- Divide each gradient in 100 bins and compute overlap with the gastric
network
4- Shift the gastric network randomly through the cortex and recompute
chance level overlap
5- Plot


Ignacio Rebollo
05-05-2021

%}


%% resample the gastric network to the same voxel resolution as Margulies gradients

% List of open inputs
nrun = 1; % enter the number of runs here
jobfile = {'/media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/batchResampleGastricnetwork2Margulies_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});


%% First gradient


% Load gastric network (in the same resolution as margulies gradient)
gasnet_mask = ft_read_mri('/media/irebollo/storage/projects/Navigastric/scripts/newMAIN_paper/otherParcellLocalizers/gradientsMargulies//media/irebollo/storage/projects/Navigastric/scripts_4_github/newMAIN_paper/otherParcellLocalizers/gradientsMargulies/gasnet_marguliesSpace.nii.nii')
gasnet_mask=gasnet_mask.anatomy;
gasindex = gasnet_mask>0.1;

% Load gradient 1
gradient1=ft_read_mri('/media/irebollo/storage/projects/Navigastric/scripts/newMAIN_paper/otherParcellLocalizers/gradientsMargulies/volume.cort.0.nii')
gradient1=gradient1.anatomy;
gradIdx=find(gradient1(:)); % 
gradient1_diff=gradient1(:)~=0;


% Conjunction between gastric network and cortical gradient
conj = (gasindex(:) + gradient1_diff(:))==2;

nBins=100; % numbers of bins into which divide the gradients
idx = linspace(min(gradient1(:)),max(gradient1(:)),nBins)

% Loop through each bin
idx_currentgradient=[]
table1_g1=zeros(length(idx)-1);
for i=1:length(idx)-1
idx_currentgradient{i} = gradient1(:)>idx(i)&gradient1(:)<idx(i+1);% find voxels of gradient corresponding to that bin
idx_currentgradient_size{i} = sum(gradient1(:)>idx(i)&gradient1(:)<idx(i+1)); % count how many voxels in that bin

table1_g1(i) = sum(conj(idx_currentgradient{i}))/idx_currentgradient_size{i}; % 

end





%% Estimate chance level overlap for first gradient



% create size gradient nonzero

gasnet_gradients =gasindex(gradient1_diff);


% permute the location of the gastric network across the cortex
nPermutations=1000
gasnet_gradients_perm=zeros(nPermutations,length(gasnet_gradients));
for i=1:nPermutations
gasnet_gradients_perm(i,:) = gasnet_gradients(randperm(length(gasnet_gradients)));
end


% preallocate. 1000 randomizations, 99 bins
table1=zeros(nPermutations,nBins-1);
table2=zeros(nPermutations,nBins-1);

for iPerm=1:nPermutations

iPerm

gasnetperm = zeros(91,109,91);
gasnetperm(gradient1_diff)=gasnet_gradients_perm(iPerm,:); % pick one permutation
gasindex_p = gasnetperm>0.1; % threshold
gasindex_masked = gasindex_p;
gasindex_masked(~gradient1_diff) = 0;
idx_gradient1_diff = find(gradient1_diff);


% conjunction surrogate gastric network and first gradient 
conj = (gasindex_p(:) + gradient1_diff(:))==2;


% Separate bins into griadient
idx = linspace(min(gradient1(:)),max(gradient1(:)),nBins);
idx_currentgradient=[];
for i=1:length(idx)-1
idx_currentgradient{i} = gradient1(:)>idx(i)&gradient1(:)<idx(i+1);
idx_currentgradient_size{i} = sum(gradient1(:)>idx(i)&gradient1(:)<idx(i+1));


% find gastric network voxels in each gradient
table1(iPerm,i) = sum(conj(idx_currentgradient{i}))/idx_currentgradient_size{i};
table2(iPerm,i) = sum(conj(idx_currentgradient{i}))/sum(conj);

end


end

table1_surr = table1;
table2_surr = table2;



figure
plot(table1_surr','ok')
hold on
plot(table1_g1,'r')
xticks(linspace(1,99,11))
xticklabels ((round(linspace(min(gradient1(:)),max(gradient1(:)),11),3)))
xlabel('gradient 1')


%% Second gradient
gradient2= ft_read_mri('/media/irebollo/storage/projects/Navigastric/scripts/newMAIN_paper/otherParcellLocalizers/gradientsMargulies/volume.cort.1.nii')
% gradient2= ft_read_mri('/media/irebollo/30C9-0B0E/current/Situating the default-mode network along a principal gradient of macroscale cortical organization/volume.cort.1.nii')

gradient2=gradient2.anatomy;

conj = (gasindex(:) + gradient1_diff(:))==2; % we can use gradient1 as this is just to remvoe non cortical voxels

idx = linspace(min(gradient2(:)),max(gradient2(:)),nBins)
idx_currentgradient=[]
idx_currentgradient_size=[]
table1_g2=zeros(length(idx)-1);
for i=1:length(idx)-1
idx_currentgradient{i} = gradient2(:)>idx(i)&gradient2(:)<idx(i+1);
idx_currentgradient_size{i} = sum(gradient2(:)>idx(i)&gradient2(:)<idx(i+1));

table1_g2(i) = sum(conj(idx_currentgradient{i}))/idx_currentgradient_size{i};

end



%% gradient 2 permute chance


nPermutations=1000

table1=zeros(nPermutations,nBins-1);
table2=zeros(nPermutations,nBins-1);

for iPerm=1:nPermutations

iPerm
% gasindex_masked(~gradient1_diff) = 0;

gasnetperm = zeros(91,109,91);
gasnetperm(gradient1_diff)=gasnet_gradients_perm(iPerm,:);
gasindex_p = gasnetperm>0.1;


gasindex_masked = gasindex_p;
gasindex_masked(~gradient1_diff) = 0;
idx_gradient1_diff = find(gradient1_diff);



conj = (gasindex_p(:) + gradient1_diff(:))==2;
%gasnet_mask(gradIdx)=0;

idx = linspace(min(gradient2(:)),max(gradient2(:)),nBins);

idx_currentgradient=[];
for i=1:length(idx)-1
idx_currentgradient{i} = gradient2(:)>idx(i)&gradient2(:)<idx(i+1);
idx_currentgradient_size{i} = sum(gradient2(:)>idx(i)&gradient2(:)<idx(i+1));

table1(iPerm,i) = sum(conj(idx_currentgradient{i}))/idx_currentgradient_size{i};
table2(iPerm,i) = sum(conj(idx_currentgradient{i}))/sum(conj);

end


end

table1_surr_g2 = table1;
table2_surr_g2 = table2;

figure
plot(table1_surr_g2','ok')
hold on
plot(table1_g2,'r')
xticks(linspace(1,99,11))
xticklabels ((round(linspace(min(gradient2(:)),max(gradient2(:)),11),3)))
xlabel('gradient 2')


%% plots 


figure

scatterhist(gradient2(conj),gradient1(conj),'MarkerSize',1,'Kernel','on','color','r')
hold on

plot(gradient2(gradIdx),gradient1(gradIdx),'ok','MarkerSize',1)
plot(gradient2(conj),gradient1(conj),'or','MarkerSize',1)

h= gcf
set(h,'PaperOrientation','landscape');
set(h,'PaperPositionMode','auto');  
%     print ('-dpdf', '-painters', '/media/Navigastric/scripts/plots&figures/GradientsGasnet')
    
figure

scatterhist(gradient2(gradIdx),gradient1(gradIdx),'MarkerSize',1,'Kernel','on','color','k')
hold on

h= gcf
set(h,'PaperOrientation','landscape');
set(h,'PaperPositionMode','auto');  


% plot mesh

  
editedParula=colormap(parula)

figure
 [H,wx, wy, N]= densityplot(gradient2(conj),gradient1(conj),'nbins',[250,250]);
I = imagesc(log(N),[-5 5]);
colormap(editedParula)

 figure
 [H,wx, wy, N]= densityplot(gradient2(:),gradient1(:),'nbins',[250,250]);
I = imagesc(log(N),[-5 5]);
colormap(editedParula)


 