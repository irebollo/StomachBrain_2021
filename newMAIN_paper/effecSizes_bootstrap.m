%{
This script will 

1-compute the effect size (Cohen's d) of coupling strenght
(empirical vs chance PLV). 
2- obtain a metric of the variability
of effect sizes by bootstrapping the effect sizes across participants
3- write effecti sizes and variability into disk as a nifti


Ignacio Rebollo
03-05-2021

%}

%% Set cfg and subjects (this reapeats from main script)


subjectsPhysiens = [8 9 10 13 15 16 17 18 19 20 21 22 23 24 26 29 31 32 33 34 35 36 37 38 39 40 41 43 44];
subjectsNavigastric = [1:10 12 13 15:19 21 22 23 25 26 27  29 30 31 32 33 34 35 36 38 91 96];
cfgMain = global_getcfgmain


inside = tools_getIndexBrain('inside'); % index of voxels inside the brain
cs_all = zeros(length(subjectsPhysiens)+length(subjectsNavigastric),length(inside)); % preallocate
cfgMain.sample = 1 % 1st sample navigastric

% load the coupling strenght images of all subjects
counter = 0
for iS=1:length(subjectsNavigastric)
    counter = counter +1
    subj_idx = subjectsNavigastric(iS);

couplingStrenghtFilename = global_filename(subj_idx,cfgMain,'couplingStrenghtFilename')

tempcs = readnifti([couplingStrenghtFilename '.nii']);
cs_all(counter,:) = tempcs(inside);
    
end

    cfgMain.sample = 2 % second sample, physiens
for iS=1:length(subjectsPhysiens)
        counter = counter +1
    subj_idx = subjectsPhysiens(iS);
couplingStrenghtFilename = global_filename(subj_idx,cfgMain,'couplingStrenghtFilename')
tempcs = readnifti([couplingStrenghtFilename '.nii']);
cs_all(counter,:) = tempcs(inside);
end

%% calculate empirical effect size

[h p ci stats_emp] = ttest(cs_all);
figure;nhist(stats_emp.tstat)
cd=  stats_emp.tstat ./sqrt(length(subjectsNavigastric)+length(subjectsPhysiens)); % cd = Cohen's d
figure;nhist(cd)



%% compute bootstraps

bootstrapped_cd = zeros(1000,length(inside));

bootstrap_samples = randi(63,63,1000);


for iBs=1:1000
    iBs
    bootstrap_cs = cs_all(bootstrap_samples(:,iBs),:);
    [h p ci stats_surr] = ttest(bootstrap_cs);
cd_surr=  stats_surr.tstat ./sqrt(length(subjectsNavigastric)+length(subjectsPhysiens));
bootstrapped_cd(iBs,:)  = cd_surr;  
end

std_bootstrapped_cd = std(bootstrapped_cd,[],1);

figure;nhist(std_bootstrapped_cd)

low_prct = prctile(bootstrapped_cd,2.5);
high_prct = prctile(bootstrapped_cd,97.5);
figure;nhist(high_prct-low_prct)

figure
plot(high_prct-low_prct,std_bootstrapped_cd,'ok');lsline


figure
plot(cd(cd>0),std_bootstrapped_cd(cd>0),'ok');lsline
vline(0.3)

hold on
plot(cd(cd>0.3),std_bootstrapped_cd(cd>0.3),'or');lsline
vline(0.3)
xlabel('Cohen D')
ylabel('Cohen D STD')

[rcdstd,p]=corrcoef(cd(cd>0.3),std_bootstrapped_cd(cd>0.3))
r_lm = rcdstd(3)
bf_lm= corrbf(r_lm,sum(cd>0.3))


[r,p]=corrcoef(cd,std_bootstrapped_cd)
r = r(3)
bf_lm= corrbf(r,length(cd))



%% write into disk 


effectsize_fn = '/media/Navigastric/ClusterResults/effecsizes/cd' % Cohen's d
empty_brain_cd = zeros(53,63,46);
empty_brain_cd(inside) = cd;
tools_writeMri(empty_brain_cd,effectsize_fn)

empty_brain_cd_td = empty_brain_cd>0.3;

variability_effectsize_fn = '/media/Navigastric/ClusterResults/effecsizes/bs_ci' % bootstraps
empty_brain = zeros(53,63,46);
empty_brain(inside) = std_bootstrapped_cd;
% empty_brain(~empty_brain_cd_td) = 0;
tools_writeMri(empty_brain,variability_effectsize_fn)
