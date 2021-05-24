subjectsPhysiens = [8 9 10 13 15 16 17 18 19 20 21 22 23 24 26 29 31 32 33 34 35 36 37 38 39 40 41 43 44];
subjectsNavigastric = [1:10 12 13 15:19 21 22 23 25 26 27  29 30 31 32 33 34 35 36 38 91 96];



dataAll_phys_R = zeros(length(subjectsPhysiens),180)
dataAll_phys_L = zeros(length(subjectsPhysiens),180)


% Z:\ClusterResults\BOLDvariablity\models\model1_genderStaiBMIGroupandTime


filename = 'D:\SADREST\files\lh.HCPMMP1.annot'
[verticesHCP, labelHCP, colortableHC_LH] = read_annotation(filename);
filename = 'D:\SADREST\files\rh.HCPMMP1.annot'
[verticesHCP, labelHCP, colortableHC_RH] = read_annotation(filename);


for iSubj = 1:length(subjectsPhysiens)
    subj_idx=subjectsPhysiens(iSubj)
    subj_path= ['Z:\dataPhysiens\subj',sprintf('%.2d',subj_idx),'\Timeseries\fsparcels\']
%     data_lh = load ([subj_path '\cs_hcp_par_RFANTS.lh.dat']);
%     data_rh = load ([subj_path '\cs_hcp_par_RFANTS.rh.dat']);
    data_lh = load ([subj_path '\cs_hcp_par_RFANTS_nt.lh.dat']);
    data_rh = load ([subj_path '\cs_hcp_par_RFANTS_nt.rh.dat']);
dataAll_phys_R(iSubj,:) = [data_rh];
dataAll_phys_L(iSubj,:) = [data_lh];

end


dataAll_navi_R = zeros(length(subjectsNavigastric),180)
dataAll_navi_L = zeros(length(subjectsNavigastric),180)

for iSubj = 1:length(subjectsNavigastric)
    subj_idx=subjectsNavigastric(iSubj)
    subj_path= ['Z:\subjects\subj',sprintf('%.2d',subj_idx),'\Timeseries\fsparcels\']
%     data_lh = load ([subj_path '\cs_hcp_par_RFANTS.lh.dat']);
%     data_rh = load ([subj_path '\cs_hcp_par_RFANTS.rh.dat']);
        data_lh = load ([subj_path '\cs_hcp_par_RFANTS_nt.lh.dat']);
    data_rh = load ([subj_path '\cs_hcp_par_RFANTS_nt.rh.dat']);
dataAll_navi_R(iSubj,:) = [data_rh];
dataAll_navi_L(iSubj,:) = [data_lh];

end

data_All_RH = [dataAll_phys_R;dataAll_navi_R]
data_All_LH = [dataAll_phys_L;dataAll_navi_L]

%% read
table = csvread('Z:\Rfun\BFglasserNONAMES.csv',1,0)

p_values = table(:,3:3:33)


p_values_bonf = p_values*11;
sum(p_values_bonf<0.05)

[Y,X]=find(p_values_bonf<0.05)
AA=[X,Y]

p_values_fdr=zeros(size(p_values))
for i=1:11
    
    
[A B C] = fdr(p_values(:,i));
p_values_fdr(:,i) = B;
end

sum(p_values_fdr<0.05)
find(p_values_fdr<0.05)
p_values(2373)