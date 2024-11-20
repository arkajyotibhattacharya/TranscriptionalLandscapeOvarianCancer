

library(DescTools)
library(parallel)


mix_mat = read.table("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/ICA/Consensus_Mix_Matrix_14112018_Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_.txt", sep = "\t", header = TRUE)

sample_annotation = read.table("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Data/Annotation\ ovarian\ cancer\ samples\ for\ plot.txt", sep = "\t", header = TRUE)

sample_annotation = sample_annotation[which(sample_annotation$Identifier%in%colnames(mix_mat)),]

t_mix_mat = t(as.matrix(mix_mat))

t_mix_mat = t_mix_mat[as.character(sample_annotation$Identifier),]

aov1 = aov(t_mix_mat[,5] ~ sample_annotation$Series)
summary(aov1)

aov2 = aov(t_mix_mat[,5] ~ sample_annotation$Type)
summary(aov2)


unique_series = as.character(unique(sample_annotation$Series))
for(i in 1:length(unique_series))
{
  new_series = ifelse(sample_annotation$Series==unique_series[i], 1,0)
  sample_annotation = cbind(sample_annotation,new_series)
  colnames(sample_annotation)[dim(sample_annotation)[2]] = unique_series[i]
}


unique_tissue_type = as.character(unique(sample_annotation$Type))
for(i in 1:length(unique_tissue_type))
{
  new_tissue_type = ifelse(sample_annotation$Type==unique_tissue_type[i], 1,0)
  sample_annotation = cbind(sample_annotation,new_tissue_type)
  colnames(sample_annotation)[dim(sample_annotation)[2]] = unique_tissue_type[i]
}


aov_pvalues = apply(t_mix_mat, 2, function(x)
  {
  aov1 = aov(as.numeric(x) ~ sample_annotation$`Serous carcinoma - low grade`)
  return(summary(aov1)[[1]][["Pr(>F)"]][1])
  
})

association_with_tissue_types = as.data.frame(matrix(NA, dim(t_mix_mat)[2], length(unique_tissue_type)))
colnames(association_with_tissue_types) = unique_tissue_type
rownames(association_with_tissue_types) = colnames(t_mix_mat)

for(i in 1:length(unique_tissue_type))
{
  association_with_tissue_types[,i] = apply(t_mix_mat, 2, function(x)
  {
    aov1 = aov(as.numeric(x) ~ sample_annotation[,unique_tissue_type[i]])
    return(summary(aov1)[[1]][["Pr(>F)"]][1])
    
  })
  
}

association_with_series = as.data.frame(matrix(NA, dim(t_mix_mat)[2], length(unique_series)))
colnames(association_with_series) = unique_series
rownames(association_with_series) = colnames(t_mix_mat)

for(i in 1:length(unique_series))
{
  association_with_series[,i] = apply(t_mix_mat, 2, function(x)
  {
    aov1 = aov(as.numeric(x) ~ sample_annotation[,unique_series[i]])
    return(summary(aov1)[[1]][["Pr(>F)"]][1])
    
  })
  
}

association_with_tissue_types = -log10(association_with_tissue_types)
association_with_series = -log10(association_with_series)

CES = paste("GEO consensus estimated source", c(1:dim(association_with_tissue_types)[1]))
association_with_tissue_types = as.data.frame(cbind(CES,association_with_tissue_types))
association_with_series = as.data.frame(cbind(CES,association_with_series))

write.table(association_with_tissue_types, file = "/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/ICA/Minus_log10pvalue_of_association_of_CES_with_tissue_types.txt", sep = "\t", row.names = FALSE)

write.table(association_with_series, file = "/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/ICA/Minus_log10pvalue_of_association_of_CES_with_series.txt", sep = "\t", row.names = FALSE)


