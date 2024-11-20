library(data.table)

# tcga_clusters = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/OV-TP.bestclus.txt'), row.names = 1)
tcga_clusters = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/TCGA_clus.membership.txt'), row.names = 1)

mixing_matrix = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/Projection_on_GEO_TCGA/Mixing_matrix_GEO_ovcar_tcga_combined_using_GEO_ovcar_TF.txt'), row.names = 1)

colnames(mixing_matrix) = gsub("\\.", "-", colnames(mixing_matrix))

common_samples = intersect(substr(colnames(mixing_matrix),1,20), substr(rownames(tcga_clusters),1,20))                           

rownames(tcga_clusters) = substr(rownames(tcga_clusters),1,20)

colnames(mixing_matrix) = substr(colnames(mixing_matrix),1,20)
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

tcga_clusters = tcga_clusters[order(tcga_clusters$clu.2, tcga_clusters$clu.3, tcga_clusters$clu.4, tcga_clusters$clu.5
                                    , tcga_clusters$clu.6, tcga_clusters$clu.7, tcga_clusters$clu.8, tcga_clusters$clu.9
                                    , tcga_clusters$clu.10),]
mixing_matrix = mixing_matrix[rownames(tcga_clusters),]

for(i in 1:9)
{
  tcga_clusters[,i] = as.factor(tcga_clusters[,i])
}


pdf("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Plots/Heatmap_of_TCs_in_all_TCGA_OV_samples.pdf", width = 9, height = 8)
pheatmap::pheatmap(
  mat               = mixing_matrix,
  color             = colorRampPalette(c("navy", "white", "red"))(100),
  breaks            = seq(-0.1, 0.1, length.out = 100),
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  annotation_row    = tcga_clusters,
  # annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  cluster_rows=FALSE, cluster_cols=FALSE
)
dev.off()

data = as.data.frame(cbind(mixing_matrix, tcga_clusters))
colnames(data) = gsub("\\.", "_", colnames(data))
kruskal_result <- kruskal.test(V76 ~ clu_2, data = data)
kruskal_result <- kruskal.test(V76 ~ clu_3, data = data)
kruskal_result$p.value
kruskal_result <- kruskal.test(V76 ~ clu_4, data = data)
kruskal_result$p.value

association_matrix = matrix(NA, 9, 14)
colnames(association_matrix) = colnames(mixing_matrix)
rownames(association_matrix) = colnames(tcga_clusters)

for(i in 1:9)
{
  for(j in 1:14)
  {
    var1 = mixing_matrix[,j]
    var2 = tcga_clusters[,i]
    data = as.data.frame(cbind(var1, var2))
    kruskal_result <- kruskal.test(var1 ~ var2, data = data)
    association_matrix[i, j] = -log10(kruskal_result$p.value)
    
    
    
  }
}

write.table(association_matrix, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/KW_between_TCGA_clusters_and_TCs_in_all_TCGA_OV_samples.txt", sep = "\t", quote = FALSE)



pdf("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Plots/Heatmap_of_KW_between_TCGA_clusters_and_TCs_in_all_TCGA_OV_samples.pdf", width = 9, height = 8)
pheatmap::pheatmap(
  mat               = association_matrix,
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

