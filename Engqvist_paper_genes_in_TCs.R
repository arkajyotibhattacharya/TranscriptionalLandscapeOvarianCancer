library(data.table)

independent_components = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/ICA/Genelevel_using_jetset_Consensus_Independent_Components_14112018_Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_reformatted_for_analyzer_tool_gsea_after_flip_03062020.txt'))
colnames(independent_components) = gsub("GEO.ovarian.cancer...consensus.estimated.source.", "TC", colnames(independent_components))
de_genes_found_in_engqvist= c('FAM20A',
                              'LAMB1',
                              
                              'CLMN',
                              
                              'CACNA1A',
                              'CACNB1',
                              'CELF4',
                              'CLMN',
                              'COL14A1',
                              'EBF4',
                              'EHF',
                              'HMGA2',
                              'IDH3B',
                              'KCNMB2',
                              'LINC00578',
                              'MAPK4',
                              'MEIS2',
                              'MTBP',
                              'MYLK2',
                              'NRSN2',
                              'PDYN',
                              'PPP1R1B',
                              'PROKR2',
                              'RASSF2',
                              'RPL22L1',
                              'TP63',
                              'ELP3',
                              
                              'AARD',
                              'CCNE1',
                              'CTCFL',
                              'LINC00578',
                              'LINC01532',
                              'RBM38',
                              'RSPO4',
                              'UQCRFS1',
                              'URI1',
                              'PDE8B',
                              
                              'ANKS1B',
                              'COLEC10',
                              'EGFEM1P',
                              'KCNMB2-AS1',
                              'LINC00578',
                              'MYEF2',
                              'PDYN',
                              'RBFOX1',
                              'SNAP25',
                              'STON2',
                              'SULF1',
                              'TSHR',
                              'TGFBR2')

common_genes = intersect(independent_components$SYMBOL,de_genes_found_in_engqvist )
TCs = paste0("TC", c(76, 14, 121, 250, 78, 320, 253, 197, 146, 239, 247, 138, 220, 166))
independent_components = independent_components[which(independent_components$SYMBOL%in%common_genes), ]

rownames(independent_components) = independent_components$SYMBOL

independent_components = independent_components[,TCs ]
custom_breaks <- sort(unique(c(seq(-5, -3, length.out = 25), 
                   seq(-3, 3, length.out = 50), 
                   seq(3, 5, length.out = 25))))
custom_colors <- c(colorRampPalette(c("blue", "white"))(100)[1:25], 
                   colorRampPalette(c("blue", "white", "red"))(200)[76:125], 
                   colorRampPalette(c("white", "red"))(100)[76:100])

pdf("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Plots/Heatmap_of_weights_of_de_genes_engqvist_in_TCs.pdf", width = 9, height = 8)
pheatmap::pheatmap(
  mat               = independent_components,
  color             = custom_colors,
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


