for i=1:length(namesTablePaper)
    
    if isnumeric(namesTablePaper{i})
        
    namesTablePaper{i}=num2str(namesTablePaper{i})
    
    end
end

for i=1:length(allTablePaper)
    
    if isnumeric(allTablePaper{i})
        
    allTablePaper{i}=num2str(allTablePaper{i})
    
    end
end


% Find the corresponding name of each StrucIndex

index_label=ones(1,length(allTablePaper))*200
for iROI=1:length(namesTablePaper)
index_label(iROI)=find(ismember(allTablePaper,namesTablePaper{iROI}))
end

index_label = index_label'


index_label=ones(1,length(allTablePaper))*200
for iROI=1:length(allTablePaper)
    
    try
index_label(iROI)=find(ismember(namesTablePaper,allTablePaper{iROI}))

    catch
        index_label(iROI)=900

    end
end

index_label = index_label'