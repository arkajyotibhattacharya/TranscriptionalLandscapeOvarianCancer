library(data.table)
geo_ovarian_tf = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/ICA/Genelevel_using_jetset_Consensus_Independent_Components_14112018_Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_after_flip.txt", sep = "\t"))
rownames(geo_ovarian_tf) = geo_ovarian_tf$ENTREZID
geo_ovarian_tf = geo_ovarian_tf[,-c(1:12)]
flip_source = read.table("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/Survival\ analysis/OS_univariate_flip_Y_N_with_series.txt", sep = "\t", header = TRUE)

flip_source_vec = ifelse(as.character(flip_source$Flip_1_yes)==1, -1, 1)

geo_ovarian_tf = t(t(geo_ovarian_tf)*flip_source_vec)



tab5rows  = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ GDSC/Data/GDSC__Affy_hgu219_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetset.txt",sep = "\t",header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
gdsc_mrna_expression_data = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ GDSC/Data/GDSC__Affy_hgu219_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetset.txt",sep = "\t",header = TRUE, colClasses = classes)

tab5rows  = read.table("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Data/Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_12112018.txt",header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
geo_ovarian_mrna_expression_data = read.table("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Data/Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_12112018.txt",header = TRUE, colClasses = classes)

probe_gene_info = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/ICA/Genelevel_using_jetset_Consensus_Independent_Components_14112018_Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_after_flip.txt", sep = "\t"))
rownames(probe_gene_info) = probe_gene_info$PROBESET

common_genes = unique(intersect(rownames(probe_gene_info), rownames(geo_ovarian_mrna_expression_data)))

geo_ovarian_mrna_expression_data = geo_ovarian_mrna_expression_data[rownames(probe_gene_info),]
rownames(geo_ovarian_mrna_expression_data) = probe_gene_info$ENTREZID
rm(probe_gene_info)



common_genes = unique(intersect(rownames(geo_ovarian_tf), rownames(gdsc_mrna_expression_data)))
common_genes = unique(intersect(common_genes, rownames(geo_ovarian_mrna_expression_data)))

geo_ovarian_tf = geo_ovarian_tf[common_genes,]
geo_ovarian_mrna_expression_data = geo_ovarian_mrna_expression_data[common_genes,]
gdsc_mrna_expression_data = gdsc_mrna_expression_data[common_genes,]


gdsc_mrna_expression_data_standardized = t(scale(t(gdsc_mrna_expression_data)))

geo_ovarian_mrna_expression_data_standardized = t(scale(t(geo_ovarian_mrna_expression_data)))


combined_expression_dataset = cbind(geo_ovarian_mrna_expression_data_standardized, gdsc_mrna_expression_data_standardized)
geo_ovarian_tf = as.matrix(geo_ovarian_tf)
dot_prod = t(geo_ovarian_tf)%*%geo_ovarian_tf
V = solve(dot_prod)
mix_matrix = V%*%t(geo_ovarian_tf)%*%as.matrix(combined_expression_dataset)
rownames(mix_matrix) = colnames(geo_ovarian_tf)

write.table(mix_matrix, file = "/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/Projection_on_GDSC/Mixing_matrix_GEO_ovcar_GDSC_combined_using_GEO_ovcar_TF.txt", sep = "\t", row.names = TRUE)
