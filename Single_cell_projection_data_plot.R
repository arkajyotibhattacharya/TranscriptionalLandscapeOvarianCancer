library(data.table)

mixing_matrix = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/single_cell_data/Project_data_on_independent_components_{e9738ccc-e862-44c4-9245-cf1826464c1d}/mixing_matrix.tsv'), row.names = 1)

cell_annotations = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Single_cell_data/GSE158722.cell_annotations.txt'))
rownames(cell_annotations) = cell_annotations$Cell

common_cells = intersect(rownames(mixing_matrix), rownames(cell_annotations))

mixing_matrix_rownames_not_present = rownames(mixing_matrix)[which(!rownames(mixing_matrix)%in%common_cells)]
cell_annotations_not_present = rownames(cell_annotations)[which(!rownames(cell_annotations)%in%common_cells)]
cell_annotations_not_present = list()

cell_annotations_not_present$Cell = mixing_matrix_rownames_not_present
cell_annotations_not_present$Normal.or.Cancer = "Unknown"
cell_annotations_not_present$SingleR_Encode_main_type = "Unannotated"
cell_annotations_not_present = as.data.frame(cell_annotations_not_present)

rownames(cell_annotations_not_present) = cell_annotations_not_present$Cell
cell_annotations$Platform = NULL
cell_annotations$Sample.name = NULL
cell_annotations$time = NULL
cell_annotations$Cluster = NULL
cell_annotations$UMAP_1 = NULL
cell_annotations$UMAP_2 = NULL

cell_annotations = as.data.frame(rbind(cell_annotations, cell_annotations_not_present))
common_cells = intersect(rownames(mixing_matrix), rownames(cell_annotations))

mixing_matrix_subset = mixing_matrix[common_cells,]
cell_annotations_subset = cell_annotations[common_cells,]
  
colnames(mixing_matrix_subset) = gsub("GEO.ovarian.cancer...consensus.estimated.source.", "TC", colnames(mixing_matrix_subset))

corrected_mixing_matrix = as.data.frame(cbind(mixing_matrix_subset,cell_annotations_subset ))

corrected_mixing_matrix = corrected_mixing_matrix[order(corrected_mixing_matrix$SingleR_Encode_main_type),]


pdf("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Plots/corrected_mixing_matrix_per_celltype.pdf", width = 36,height = 11)
par(mar = c(23, 5, 4, 2) )
par(mfrow = c(1,3))
for(i in 1:14)
{
  print(boxplot(as.formula(paste(colnames(corrected_mixing_matrix)[i], " ~ SingleR_Encode_main_type")), data = corrected_mixing_matrix
                , ylim = c(-5,5)
                , main = colnames(corrected_mixing_matrix)[i], ylab = "Activity scores",xlab = "", las = 2, col = "lightblue", notch = TRUE, pch = 19))
}

dev.off()


pdf("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Plots/corrected_mixing_matrix_per_normal_or_cancer.pdf", width = 36,height = 11)
par(mar = c(23, 5, 4, 2) )
par(mfrow = c(1,3))
for(i in 1:14)
{
  print(boxplot(as.formula(paste(colnames(corrected_mixing_matrix)[i], " ~ Normal.or.Cancer")), data = corrected_mixing_matrix
                , ylim = c(-5,5)
                , main = colnames(corrected_mixing_matrix)[i], ylab = "Activity scores",xlab = "", las = 2, col = "lightblue", notch = TRUE, pch = 19))
}

dev.off()

pdf("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Plots/corrected_mixing_matrix_per_celltype_all_TCs.pdf", width = 36,height = 11)
par(mar = c(23, 5, 4, 2) )
par(mfrow = c(1,3))
unique_cell_types= unique(corrected_mixing_matrix$SingleR_Encode_main_type)
for(i in 1:length(unique_cell_types))
{
  print(boxplot(corrected_mixing_matrix[which(corrected_mixing_matrix$SingleR_Encode_main_type==unique_cell_types[i]),1:14]
                , ylim = c(-5,5)
                , main = unique_cell_types[i], ylab = "Activity scores",xlab = "", las = 2, col = "lightblue", notch = TRUE, pch = 19))
  # print(boxplot(as.formula(paste(colnames(corrected_mixing_matrix)[i], " ~ Normal.or.Cancer")), data = corrected_mixing_matrix
  #               , ylim = c(-5,5)
  #               , main = colnames(corrected_mixing_matrix)[i], ylab = "Activity scores",xlab = "", las = 2, col = "lightblue", notch = TRUE, pch = 19))
}

dev.off()

boxplot(corrected_mixing_matrix[,1:14])

