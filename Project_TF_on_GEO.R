library(data.table)
tab5rows  = read.table("Downloads/Consensus_Independent_Components_14112018_Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_.txt", sep = "\t" , header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
geo_ovarian_tf = read.table("Downloads/Consensus_Independent_Components_14112018_Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_.txt", sep = "\t", header = TRUE, colClasses = classes)
flip_source = read.table("Downloads/OS_univariate_flip_Y_N_with_series.txt", sep = "\t", header = TRUE)

flip_source_vec = ifelse(as.character(flip_source$Flip_1_yes)==1, -1, 1)

geo_ovarian_tf = t(t(geo_ovarian_tf)*flip_source_vec)


geo_mrna_expression_data = as.data.frame(fread("Downloads/GPL570__Affy_hgu133plus2_NoDuplicatesamples_CleanedIdentifiers_RMA-sketch_NormalAndCancer_QCed.txt"))
rownames(geo_mrna_expression_data) = geo_mrna_expression_data$V1
geo_mrna_expression_data$V1 = NULL
tab5rows  = read.table("Downloads/Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_12112018.txt",header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
geo_ovarian_mrna_expression_data = read.table("Downloads/Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_12112018.txt",header = TRUE, colClasses = classes)

common_probesets = unique(intersect(rownames(geo_ovarian_tf), rownames(geo_mrna_expression_data)))
common_probesets = unique(intersect(common_probesets, rownames(geo_ovarian_mrna_expression_data)))

geo_ovarian_tf = geo_ovarian_tf[common_probesets,]
geo_ovarian_mrna_expression_data = geo_ovarian_mrna_expression_data[common_probesets,]
geo_mrna_expression_data = geo_mrna_expression_data[common_probesets,]

# geo_mrna_expression_data_standardized = geo_mrna_expression_data
# 
# for(i in 1:dim(geo_mrna_expression_data)[1])
# {
#   a = geo_mrna_expression_data[i,]
#   geo_mrna_expression_data_standardized[i,] = scale(t(a))
#   print(i)
# }

geo_mrna_expression_data_standardized = t(scale(t(geo_mrna_expression_data)))

geo_ovarian_mrna_expression_data_standardized = t(scale(t(geo_ovarian_mrna_expression_data)))


combined_expression_dataset = cbind(geo_ovarian_mrna_expression_data_standardized, geo_mrna_expression_data_standardized)
geo_ovarian_tf = as.matrix(geo_ovarian_tf)
dot_prod = t(geo_ovarian_tf)%*%geo_ovarian_tf
V = solve(dot_prod)
mix_matrix = V%*%t(geo_ovarian_tf)%*%as.matrix(combined_expression_dataset)
rownames(mix_matrix) = colnames(geo_ovarian_tf)

write.table(mix_matrix, file = "Downloads/Mixing_matrix_GEO_ovcar_geo_combined_using_GEO_ovcar_TF.txt", sep = "\t", row.names = TRUE)
