library(data.table)
library(survival)
tcga_survival_data = data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/denseDataOnlyDownload.tsv"), row.names = 1)
class(tcga_survival_data$OS)
sum(is.na(tcga_survival_data$OS))
class(tcga_survival_data$OS.time)
sum(is.na(tcga_survival_data$OS.time))
tcga_survival_data = tcga_survival_data[which(!is.na(tcga_survival_data$OS)),]
tcga_survival_data = tcga_survival_data[which(!is.na(tcga_survival_data$OS.time)),]

tcga_activity_scores = data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/Projection_on_GEO_TCGA/Mixing_matrix_GEO_ovcar_tcga_combined_using_GEO_ovcar_TF.txt"), row.names = 1)
tcga_activity_scores = tcga_activity_scores[,c(1126:dim(tcga_activity_scores)[2])]

colnames(tcga_activity_scores) = gsub("\\.", "-", colnames(tcga_activity_scores))
colnames(tcga_activity_scores) = substr(colnames(tcga_activity_scores),1,15)
tcga_activity_scores = as.data.frame(t(tcga_activity_scores))
rownames(tcga_activity_scores) = gsub("\\.", "-", rownames(tcga_activity_scores))

common_samples = intersect(rownames(tcga_activity_scores), rownames(tcga_survival_data))

tcga_activity_scores = tcga_activity_scores[common_samples,]
tcga_survival_data = tcga_survival_data[common_samples,]

unique_subtypes = sort(unique(tcga_survival_data$cancer.type.abbreviation))

summary_matrix = matrix(NA, dim(tcga_activity_scores)[2], length(unique_subtypes))
colnames(summary_matrix) = unique_subtypes
rownames(summary_matrix) = paste("TC", 1:dim(summary_matrix)[1], sep = "")

for( j in 1:length(unique_subtypes))
{
  current_subtype = unique_subtypes[j]
  data = tcga_activity_scores[ rownames(tcga_survival_data)[which(tcga_survival_data$cancer.type.abbreviation==current_subtype)],]
  follow_up= tcga_survival_data$OS.time[which(tcga_survival_data$cancer.type.abbreviation==current_subtype)]
  event = as.numeric(tcga_survival_data$OS)[which(tcga_survival_data$cancer.type.abbreviation==current_subtype)]
  outcome_variable_range=c(1:374)
  list_of_predictors=names(data[outcome_variable_range])
  newdata = data[,outcome_variable_range]						
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
  
  summary_matrix[,j] = -log10(full_coef_matrix[,"Pr(>|z|)"])*sign(full_coef_matrix[,"coef"])
  
  print(j)
}



write.table(summary_matrix, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/Projection_on_GEO_TCGA/Association_with_OS_and_TC_activity_TCGA_all_cancer_types.txt", sep = "\t", quote = FALSE)
summary_matrix = t(summary_matrix)
write.table(summary_matrix, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/Projection_on_GEO_TCGA/Association_with_OS_and_TC_activity_TCGA_all_cancer_types_transposed.txt", sep = "\t", quote = FALSE)

summary_matrix = data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/Projection_on_GEO_TCGA/Association_with_OS_and_TC_activity_TCGA_all_cancer_types_transposed.txt"), row.names = 1)
TCs = c(76, 14, 121, 250, 78, 320, 253, 197, 146, 239, 247, 138, 220, 166)

summary_matrix = summary_matrix[,TCs]
custom_breaks <- sort(unique(c(seq(-6, -2, length.out = 25), 
                               seq(-2, 2, length.out = 50), 
                               seq(2, 6, length.out = 25))))
custom_colors <- c(colorRampPalette(c("blue", "white"))(100)[1:25], 
                   colorRampPalette(c("blue", "white", "red"))(200)[76:125], 
                   colorRampPalette(c("white", "red"))(100)[76:100])

pdf("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Plots/Univariate_OS_association_TCGA_all_cancer_types.pdf", width = 9, height = 8)
pheatmap::pheatmap(
  mat               = summary_matrix,
  # color             = custom_colors,
  breaks            = custom_breaks,
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  # annotation_row    = tcga_clusters,
  # annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  cluster_rows=TRUE, cluster_cols=FALSE
)
dev.off()


