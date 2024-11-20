library(data.table)
tab5rows  = read.table("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/ICA/Consensus_Independent_Components_14112018_Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_.txt", sep = "\t" , header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
geo_ovarian_tf = read.table("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/ICA/Consensus_Independent_Components_14112018_Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_.txt", sep = "\t", header = TRUE, colClasses = classes)
flip_source = read.table("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/Survival\ analysis/OS_univariate_flip_Y_N_with_series.txt", sep = "\t", header = TRUE)

flip_source_vec = ifelse(as.character(flip_source$Flip_1_yes)==1, -1, 1)

geo_ovarian_tf = t(t(geo_ovarian_tf)*flip_source_vec)



tab5rows  = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch.txt",header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
ccle_mrna_expression_data = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch.txt",header = TRUE, colClasses = classes)
tab5rows  = read.table("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Data/Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_12112018.txt",header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
geo_ovarian_mrna_expression_data = read.table("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Data/Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_12112018.txt",header = TRUE, colClasses = classes)

common_probesets = unique(intersect(rownames(geo_ovarian_tf), rownames(ccle_mrna_expression_data)))
common_probesets = unique(intersect(common_probesets, rownames(geo_ovarian_mrna_expression_data)))

geo_ovarian_tf = geo_ovarian_tf[common_probesets,]
geo_ovarian_mrna_expression_data = geo_ovarian_mrna_expression_data[common_probesets,]
ccle_mrna_expression_data = ccle_mrna_expression_data[common_probesets,]


ccle_mrna_expression_data_standardized = t(scale(t(ccle_mrna_expression_data)))

geo_ovarian_mrna_expression_data_standardized = t(scale(t(geo_ovarian_mrna_expression_data)))


combined_expression_dataset = cbind(geo_ovarian_mrna_expression_data_standardized, ccle_mrna_expression_data_standardized)
geo_ovarian_tf = as.matrix(geo_ovarian_tf)
dot_prod = t(geo_ovarian_tf)%*%geo_ovarian_tf
V = solve(dot_prod)
mix_matrix = V%*%t(geo_ovarian_tf)%*%as.matrix(combined_expression_dataset)
rownames(mix_matrix) = colnames(geo_ovarian_tf)

write.table(mix_matrix, file = "/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/Projection_on_CCLE_DEPMAP/Mixing_matrix_GEO_ovcar_CCLE_combined_using_GEO_ovcar_TF.txt", sep = "\t", row.names = TRUE)
