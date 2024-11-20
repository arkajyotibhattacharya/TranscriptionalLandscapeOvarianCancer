
library(Seurat)
library(data.table)

cell_annotations = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Single_cell_data/GSE158722.cell_annotations.txt'))
rownames(cell_annotations) = cell_annotations$Cell

full_data_rownames = NULL
full_data_colnames = NULL

for(i in c(1:21,23:24))
{
  if(i ==1)
  {
    current_data = data.frame(fread(paste0("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Single_cell_data/GSE158722_P0",i,".counts.txt.gz")), row.names = 1)
    print(length(which(current_data$Gene.ID=="")))
    current_data = current_data[which(current_data$Gene.ID!=""),]
    current_data = current_data[which(current_data$Gene.ID%in%names(which(table(current_data$Gene.ID)==1))),]
    rownames(current_data) = current_data$Gene.ID
    current_data$Gene.ID = NULL
    current_data$Gene.Symbol = NULL
    
    common_profiles = colnames(current_data)[which(!(colnames(current_data) %in%rownames(cell_annotations)))]
    current_data = current_data[,common_profiles]
    
    full_data_rownames = rownames(current_data)
    full_data_colnames = colnames(current_data)
    print(head(rownames(current_data)))
  }else if(i <10)
  {
    current_data = data.frame(fread(paste0("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Single_cell_data/GSE158722_P0",i,".counts.txt.gz")), row.names = 1)
    print(length(which(current_data$Gene.ID=="")))
    current_data = current_data[which(current_data$Gene.ID!=""),]
    current_data = current_data[which(current_data$Gene.ID%in%names(which(table(current_data$Gene.ID)==1))),]
    rownames(current_data) = current_data$Gene.ID
    current_data$Gene.ID = NULL
    current_data$Gene.Symbol = NULL
    
    common_profiles = colnames(current_data)[which(!(colnames(current_data) %in%rownames(cell_annotations)))]
    current_data = current_data[,common_profiles]
    
    
    full_data_rownames = unique(c(rownames(current_data), full_data_rownames))
    full_data_colnames = unique(c(colnames(current_data), full_data_colnames))
    print(head(rownames(current_data)))
  }else{
    current_data = data.frame(fread(paste0("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Single_cell_data/GSE158722_P",i,".counts.txt.gz")), row.names = 1)
    print(length(which(current_data$Gene.ID=="")))
    current_data = current_data[which(current_data$Gene.ID!=""),]
    current_data = current_data[which(current_data$Gene.ID%in%names(which(table(current_data$Gene.ID)==1))),]
    rownames(current_data) = current_data$Gene.ID
    current_data$Gene.ID = NULL
    current_data$Gene.Symbol = NULL
    
    common_profiles = colnames(current_data)[which(!(colnames(current_data) %in%rownames(cell_annotations)))]
    current_data = current_data[,common_profiles]
    
    full_data_rownames = unique(c(rownames(current_data), full_data_rownames))
    full_data_colnames = unique(c(colnames(current_data), full_data_colnames))
    print(head(rownames(current_data)))
  }
  
  print(i)
}


n_cols = round(length(full_data_colnames)*0.1)

full_data = as.data.frame(matrix(NA, length(full_data_rownames), n_cols))
rownames(full_data) = full_data_rownames

set.seed(1234)
colnames(full_data) = sample(full_data_colnames, n_cols, replace = FALSE)

for(i in c(1:21,23:24))
{
  if(i ==1)
  {
    current_data = data.frame(fread(paste0("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Single_cell_data/GSE158722_P0",i,".counts.txt.gz")), row.names = 1)
    length(which(current_data$Gene.ID==""))
    current_data = current_data[which(current_data$Gene.ID!=""),]
    current_data = current_data[which(current_data$Gene.ID%in%names(which(table(current_data$Gene.ID)==1))),]
    rownames(current_data) = current_data$Gene.ID
    current_data$Gene.ID = NULL
    current_data$Gene.Symbol = NULL
    
    common_samples = intersect(colnames(full_data), colnames(current_data))
    current_data = current_data[,common_samples]
    full_data[rownames(current_data),  colnames(current_data)] = current_data
    
  }else if(i <10)
  {
    current_data = data.frame(fread(paste0("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Single_cell_data/GSE158722_P0",i,".counts.txt.gz")), row.names = 1)
    length(which(current_data$Gene.ID==""))
    current_data = current_data[which(current_data$Gene.ID!=""),]
    current_data = current_data[which(current_data$Gene.ID%in%names(which(table(current_data$Gene.ID)==1))),]
    rownames(current_data) = current_data$Gene.ID
    current_data$Gene.ID = NULL
    current_data$Gene.Symbol = NULL
    
    common_samples = intersect(colnames(full_data), colnames(current_data))
    current_data = current_data[,common_samples]
    full_data[rownames(current_data),  colnames(current_data)] = current_data
    
  }else{
    current_data = data.frame(fread(paste0("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Single_cell_data/GSE158722_P",i,".counts.txt.gz")), row.names = 1)
    length(which(current_data$Gene.ID==""))
    current_data = current_data[which(current_data$Gene.ID!=""),]
    current_data = current_data[which(current_data$Gene.ID%in%names(which(table(current_data$Gene.ID)==1))),]
    rownames(current_data) = current_data$Gene.ID
    current_data$Gene.ID = NULL
    current_data$Gene.Symbol = NULL
    
    common_samples = intersect(colnames(full_data), colnames(current_data))
    current_data = current_data[,common_samples]
    full_data[rownames(current_data),  colnames(current_data)] = current_data
    
  }
  
  print(i)
}

# write.table(full_data, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Single_cell_data/Ten_percent_subset_of_single_cells_ovarian_cancer.txt", sep = "\t", quote = FALSE)

na_colwise = apply(full_data,2, function(x){sum(is.na(x))})
table(na_colwise)
na_rowwise = apply(full_data,1, function(x){sum(is.na(x))})
table(na_rowwise)

full_data = full_data[which(na_rowwise==0),]

sum(is.na(full_data))

# write.table(full_data, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Single_cell_data/Ten_percent_subset_of_single_cells_ovarian_cancer.txt", sep = "\t", quote = FALSE)

sd_colwise = apply(full_data,2, function(x){sd(x)})
table(sd_colwise)
sd_rowwise = apply(full_data,1, function(x){sd(x)})
table(sd_rowwise)
length(which(sd_rowwise==0))
length(which(sd_colwise==0))

full_data = full_data[which(sd_rowwise!=0),]


write.table(full_data, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Single_cell_data/Ten_percent_subset_of_single_cells_ovarian_cancer_unannotated.txt", sep = "\t", quote = FALSE)

full_data_annnotated = data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Single_cell_data/Ten_percent_subset_of_single_cells_ovarian_cancer.txt"), row.names = 1)

common_genes = intersect(rownames(full_data), rownames(full_data_annnotated))
full_data_annnotated = full_data_annnotated[common_genes,]
full_data = full_data[common_genes,]

combined_data = as.data.frame(cbind(full_data, full_data_annnotated))
write.table(combined_data, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Single_cell_data/Ten_percent_subset_of_single_cells_ovarian_cancer_annotated_unannotated_combined.txt", sep = "\t", quote = FALSE)

