library(data.table)
tab5rows  = read.table("Downloads/Consensus_Independent_Components_14112018_Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_.txt", sep = "\t" , header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
geo_ovarian_tf = read.table("Downloads/Consensus_Independent_Components_14112018_Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_.txt", sep = "\t", header = TRUE, colClasses = classes)
flip_source = read.table("Downloads/OS_univariate_flip_Y_N_with_series.txt", sep = "\t", header = TRUE)

flip_source_vec = ifelse(as.character(flip_source$Flip_1_yes)==1, -1, 1)

geo_ovarian_tf = t(t(geo_ovarian_tf)*flip_source_vec)

genomic_mapping_file = as.data.frame(fread("Downloads/Genomic_Mapping_hgu133plus2_using_jetscore_3003201_v1.txt"))
genomic_mapping_file = genomic_mapping_file[which(genomic_mapping_file$top_probe_indicator==1),]

geo_ovarian_tf_entrezid = geo_ovarian_tf[genomic_mapping_file$PROBESET,]
rownames(geo_ovarian_tf_entrezid) = genomic_mapping_file$ENTREZID

tcga_mrna_expression_data = as.data.frame(fread("Downloads/TCGA__RSEM_genes_RNAseq__duplicate_samples_removed__genes_with_all_zeroes_removed.txt"))
rownames(tcga_mrna_expression_data) = tcga_mrna_expression_data$V1
tcga_mrna_expression_data$V1 = NULL
tab5rows  = read.table("Downloads/Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_12112018.txt",header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
geo_ovarian_mrna_expression_data = read.table("Downloads/Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_12112018.txt",header = TRUE, colClasses = classes)
geo_ovarian_mrna_expression_data = geo_ovarian_mrna_expression_data[genomic_mapping_file$PROBESET,]
rownames(geo_ovarian_mrna_expression_data) = genomic_mapping_file$ENTREZID

common_probesets = unique(intersect(rownames(geo_ovarian_tf_entrezid), rownames(tcga_mrna_expression_data)))
common_probesets = unique(intersect(common_probesets, rownames(geo_ovarian_mrna_expression_data)))

geo_ovarian_tf_entrezid = geo_ovarian_tf_entrezid[common_probesets,]
geo_ovarian_mrna_expression_data = geo_ovarian_mrna_expression_data[common_probesets,]
tcga_mrna_expression_data = tcga_mrna_expression_data[common_probesets,]

# tcga_mrna_expression_data_standardized = tcga_mrna_expression_data
# 
# for(i in 1:dim(tcga_mrna_expression_data)[1])
# {
#   a = tcga_mrna_expression_data[i,]
#   tcga_mrna_expression_data_standardized[i,] = scale(t(a))
#   print(i)
# }

tcga_mrna_expression_data_standardized = t(scale(t(tcga_mrna_expression_data)))

geo_ovarian_mrna_expression_data_standardized = t(scale(t(geo_ovarian_mrna_expression_data)))


combined_expression_dataset = cbind(geo_ovarian_mrna_expression_data_standardized, tcga_mrna_expression_data_standardized)
geo_ovarian_tf = as.matrix(geo_ovarian_tf)
dot_prod = t(geo_ovarian_tf)%*%geo_ovarian_tf
V = solve(dot_prod)
mix_matrix = V%*%t(geo_ovarian_tf_entrezid)%*%as.matrix(combined_expression_dataset)
rownames(mix_matrix) = colnames(geo_ovarian_tf_entrezid)

write.table(mix_matrix, file = "Downloads/Mixing_matrix_GEO_ovcar_tcga_combined_using_GEO_ovcar_TF.txt", sep = "\t", row.names = TRUE)
