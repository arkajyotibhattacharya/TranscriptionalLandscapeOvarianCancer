library(data.table)
library(survival)
expression_data = readRDS('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Bolton_OCCC-master/processing/CCOC.RNAseq.normalized.rds')
gene_info = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Laptop backup/Work outside projects 09102022/Work outside projects/Tertiary Lymphoid structures/Data/hgnc_gene_anno_merged.tsv'))
gene_info <- gene_info[!duplicated(gene_info$ensembl_id) & !duplicated(gene_info$ensembl_id, fromLast = TRUE), ]
gene_info <- gene_info[!duplicated(gene_info$NCBI_id) & !duplicated(gene_info$NCBI_id, fromLast = TRUE), ]
rownames(gene_info) = gene_info$ensembl_id
length(unique(gene_info$NCBI_id))

rownames_expression_data = sapply(rownames(expression_data), function(x){strsplit(x, "\\.")[[1]][1]})

expression_data_subset = as.data.frame(expression_data)

expression_data_subset$rownames_expression_data = rownames_expression_data

expression_data_subset <- expression_data_subset[!duplicated(expression_data_subset$rownames_expression_data) & !duplicated(expression_data_subset$rownames_expression_data, fromLast = TRUE), ]

rownames(expression_data_subset) = expression_data_subset$rownames_expression_data
expression_data_subset$rownames_expression_data = NULL
common_genes = intersect(rownames(expression_data_subset) ,gene_info$ensembl_id )
expression_data_subset = expression_data_subset[common_genes,]
gene_info = gene_info[common_genes,]
rownames(expression_data_subset) = gene_info$NCBI_id
expression_data = expression_data_subset

sd_genes = apply(expression_data,1,sd)
length(which(sd_genes==0))
expression_data = expression_data[which(sd_genes!=0),]
sd_samples = apply(expression_data,2,sd)
length(which(sd_samples==0))
expression_data = expression_data[,which(sd_samples!=0)]

expression_data_sample_id = sapply(colnames(expression_data), function(x){strsplit(x, "_")[[1]][1]})

length(unique(expression_data_sample_id))

clinical_info = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Bolton_OCCC-master/processing/clinical_fin_minimum.txt'))

clinical_info_subset = clinical_info[which(clinical_info$leukgen_individualid%in%expression_data_sample_id),]

independent_components = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/ICA/Genelevel_using_jetset_Consensus_Independent_Components_14112018_Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_reformatted_for_analyzer_tool_gsea_after_flip_03062020.txt'))
colnames(independent_components) = gsub("GEO.ovarian.cancer...consensus.estimated.source.", "TC", colnames(independent_components))
rownames(independent_components) = independent_components$ENTREZID
TCs = paste0("TC", c(76, 14, 121, 250, 78, 320, 253, 197, 146, 239, 247, 138, 220, 166))
independent_components = independent_components[, TCs]


tab5rows  = read.table("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian\ cancer\ project/Data/Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_12112018.txt",header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
geo_ovarian_mrna_expression_data = read.table("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian\ cancer\ project/Data/Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_12112018.txt",header = TRUE, colClasses = classes)

probe_gene_info = as.data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian\ cancer\ project/Results/ICA/Genelevel_using_jetset_Consensus_Independent_Components_14112018_Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_after_flip.txt", sep = "\t"))
rownames(probe_gene_info) = probe_gene_info$PROBESET

common_genes = unique(intersect(rownames(probe_gene_info), rownames(geo_ovarian_mrna_expression_data)))

geo_ovarian_mrna_expression_data = geo_ovarian_mrna_expression_data[rownames(probe_gene_info),]
rownames(geo_ovarian_mrna_expression_data) = probe_gene_info$ENTREZID
rm(probe_gene_info)



common_genes = unique(intersect(rownames(independent_components), rownames(expression_data)))
common_genes = unique(intersect(common_genes, rownames(geo_ovarian_mrna_expression_data)))

independent_components = independent_components[common_genes,]
geo_ovarian_mrna_expression_data = geo_ovarian_mrna_expression_data[common_genes,]
expression_data = expression_data[common_genes,]


expression_data_standardized = t(scale(t(expression_data)))

geo_ovarian_mrna_expression_data_standardized = t(scale(t(geo_ovarian_mrna_expression_data)))


combined_expression_dataset = cbind(geo_ovarian_mrna_expression_data_standardized, expression_data_standardized)
independent_components = as.matrix(independent_components)
dot_prod = t(independent_components)%*%independent_components
V = solve(dot_prod)
mix_matrix = V%*%t(independent_components)%*%as.matrix(combined_expression_dataset)
rownames(mix_matrix) = colnames(independent_components)

mix_matrix = mix_matrix[, which(colnames(mix_matrix)%in%colnames(expression_data_standardized))]

write.table(mix_matrix, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian\ cancer\ project/Results/Projection_on_OCCC_Mixing_matrix_GEO_ovcar_OCCC_samples.txt", sep = "\t", row.names = TRUE)


colnames(mix_matrix) = sapply(colnames(mix_matrix), function(x){strsplit(x, "_")[[1]][1]})
common_samples = intersect(colnames(mix_matrix), clinical_info_subset$leukgen_individualid)
rownames(clinical_info_subset) = clinical_info_subset$leukgen_individualid
clinical_info_subset = clinical_info_subset[common_samples,]
mix_matrix = mix_matrix[,common_samples]

table(clinical_info_subset$primary_therapy_outcome)


data = t(mix_matrix)
follow_up= clinical_info_subset$timelastfu
event = ifelse(clinical_info_subset$vitalstatus=="alive",0,1)
list_of_predictors=colnames(data)
newdata = as.data.frame(data)			
FDR = 0.01
nPerm = 10000
conf_level = 0.8


# Coxph loop
result <- lapply(newdata, function(x){coxph(as.formula(paste("Surv(follow_up,event)~",paste("x"))))})

# Create matrix for results  
full_coef_matrix = NULL
for(i in names(result)){
  each_result_summary = summary(result[[i]])
  full_coef_matrix = rbind(full_coef_matrix, cbind(each_result_summary$coefficients,each_result_summary$conf.int)[1,])
  rownames(full_coef_matrix)[dim(full_coef_matrix)[1]] = i}

write.table(full_coef_matrix, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian\ cancer\ project/Results/Association_with_OS_based_on_Projection_on_OCCC_Mixing_matrix_GEO_ovcar_OCCC_samples.txt", sep = "\t", row.names = TRUE)


