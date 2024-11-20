library(data.table)
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6130818/#S9
# https://bitbucket.org/britroc/cnsignatures
# tcga_clusters = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/OV-TP.bestclus.txt'), row.names = 1)
tcga_clusters = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/CN signature smith paper on TCGA/matrix_to_compare.txt'))
colnames(tcga_clusters) = gsub("\\.", "-", colnames(tcga_clusters))
tcga_clusters = t(tcga_clusters)
mixing_matrix = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/Projection_on_GEO_TCGA/Mixing_matrix_GEO_ovcar_tcga_combined_using_GEO_ovcar_TF.txt'), row.names = 1)

colnames(mixing_matrix) = gsub("\\.", "-", colnames(mixing_matrix))

common_samples = intersect(substr(colnames(mixing_matrix),1,19), substr(rownames(tcga_clusters),1,19))                           

rownames(tcga_clusters) = substr(rownames(tcga_clusters),1,19)

colnames(mixing_matrix) = substr(colnames(mixing_matrix),1,19)
mixing_matrix = mixing_matrix[, common_samples]
tcga_clusters = tcga_clusters[common_samples,]

# tcga_clusters_all = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/TCGA_clus.membership.txt'), row.names = 1)
# rownames(tcga_clusters_all) = substr(rownames(tcga_clusters_all),1,20)
# tcga_clusters_all = tcga_clusters_all[common_samples,]
# 
# table(tcga_clusters$cluster, tcga_clusters_all$clu.3)

TCs = c(76, 14, 121, 250, 78, 320, 253, 197, 146, 239, 247, 138, 220, 166)

mixing_matrix = mixing_matrix[ TCs,]
mixing_matrix = as.data.frame(t(mixing_matrix))

# tcga_clusters = tcga_clusters[order(tcga_clusters$clu.2, tcga_clusters$clu.3, tcga_clusters$clu.4, tcga_clusters$clu.5
#                                     , tcga_clusters$clu.6, tcga_clusters$clu.7, tcga_clusters$clu.8, tcga_clusters$clu.9
#                                     , tcga_clusters$clu.10),]
# mixing_matrix = mixing_matrix[rownames(tcga_clusters),]

# for(i in 1:9)
# {
#   tcga_clusters[,i] = as.factor(tcga_clusters[,i])
# }


cor_mat = matrix(NA, 7, 14)

for(i in 1:7)
{
  for(j in 1:14)
  {
    cor_mat[i,j] = cor(mixing_matrix[,j], tcga_clusters[,i], method = "spearman")
  }
}

colnames(cor_mat) = colnames(mixing_matrix)
rownames(cor_mat) = paste0("CN_signature_", c(1:7))

pdf("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Plots/Heatmap_of_correlation_between_CN_signatures_and_TCs_in_all_TCGA_OV_samples.pdf", width = 9, height = 8)
pheatmap::pheatmap(
  mat               = cor_mat,
  # color             = colorRampPalette(c("navy", "white", "red"))(100),
  # breaks            = seq(-0.1, 0.1, length.out = 100),
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  # annotation_row    = tcga_clusters,
  # annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  cluster_rows=FALSE, cluster_cols=FALSE
)
dev.off()




