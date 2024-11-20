library(parallel)
library(data.table)
library(sqldf)
GEO_sample_info = as.data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Annotations_sample_to_series_subtype.txt', sep = "\t", header = TRUE))
GEO_sample_info = GEO_sample_info[which(GEO_sample_info$Type_updated!="Normal"),]

tab5rows  = read.table('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/ICA/Consensus_Mix_Matrix_14112018_Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_reformatted_after_flip.txt', sep = "\t" , header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
geo_ovarian_mm = read.table('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/ICA/Consensus_Mix_Matrix_14112018_Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_reformatted_after_flip.txt', sep = "\t", header = TRUE, colClasses = classes)
rownames(geo_ovarian_mm) = geo_ovarian_mm$CES
geo_ovarian_mm$CES = NULL


common_id = intersect(colnames(geo_ovarian_mm), GEO_sample_info$Identifier)
rownames(GEO_sample_info) = GEO_sample_info$Identifier

GEO_sample_info = GEO_sample_info[common_id,]
geo_ovarian_mm = geo_ovarian_mm[,common_id]


geo_ovarian_mm_t = as.data.frame(t(geo_ovarian_mm))
geo_ovarian_mm_t$cancer_type = GEO_sample_info$Type_updated
table(geo_ovarian_mm_t$cancer_type )
library(ggplot2)

fac <- with(geo_ovarian_mm_t, reorder(cancer_type, V121, median, order = TRUE))
geo_ovarian_mm_t$cancer_type <- factor(geo_ovarian_mm_t$cancer_type, levels = levels(fac))

pdf(paste('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Plots/OVcar_all_subtypes_V121_mm_projection_distribution.pdf', sep = "")
    , height = 4, width = 8)
ggplot(geo_ovarian_mm_t, aes(x=cancer_type, y=V121)) + 
  geom_dotplot(binaxis="y",stackdir="center", dotsize=0.1, color=rgb(0,0,0,0.5), binwidth = 0.0015) + 
  geom_boxplot(fill=rgb(0,0,0,0.5), color=rgb(0,0,0,0.5), alpha=0.4) + scale_y_continuous(name = "Activity score", limits = c(-0.15,0.15), breaks = c(-0.15,-0.1,-0.05,0, 0.05, 0.1, 0.15)) + 
  scale_x_discrete(name = "") + 
  # geom_hline(yintercept=0) + 
  theme_classic()+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  ggtitle(paste("TC121 activity scores"))+ coord_flip()
dev.off()



