library(data.table)
library(survival)
library(MST)
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

# summary_matrix = matrix(NA, dim(tcga_activity_scores)[2], length(unique_subtypes))
# colnames(summary_matrix) = unique_subtypes
# rownames(summary_matrix) = paste("TC", 1:dim(summary_matrix)[1], sep = "")
full_tree = list()
for( j in 1:length(unique_subtypes))
{
  current_subtype = unique_subtypes[j]
  data = tcga_activity_scores[ rownames(tcga_survival_data)[which(tcga_survival_data$cancer.type.abbreviation==current_subtype)],]
  follow_up= tcga_survival_data$OS.time[which(tcga_survival_data$cancer.type.abbreviation==current_subtype)]
  event = as.numeric(tcga_survival_data$OS)[which(tcga_survival_data$cancer.type.abbreviation==current_subtype)]
  outcome_variable_range=c(1:374)
  list_of_predictors=names(data[outcome_variable_range])
  newdata = as.data.frame(cbind(follow_up, cbind(event,data[,outcome_variable_range])))						
  newdata$identifier = rownames(newdata)
  
  minsplits = 50
  
  full_tree[[j]] <- MST(formula = Surv(follow_up, event) ~ V14
                   +V76
                   +V78
                   +V121
                   +V138
                   +V146
                   +V197
                   +V220
                   +V239
                   +V247
                   +V250
                   +V253
                   +V320
                   +V166
                   | identifier,
                   data = newdata,
                   test = newdata,
                   method = "independence",
                   minsplit = minsplits,
                   minevents = ceiling(minsplits/2),
                   minbucket = ceiling(minsplits/3),
                   selection.method = "test.sample",
                   # LeBlanc = TRUE,
                   plot.Ga = TRUE,
                   sortTrees = TRUE,
                   details = FALSE)
  
  (tree_final <- getTree(full_tree[[j]], "0"))
  
  plot(tree_final)
  print(j)
  
}



save(full_tree, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/Projection_on_GEO_TCGA/Association_with_OS_and_TC_activity_TCGA_all_cancer_types_all_trees.RDS")

pdf("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Plots/all_trees.pdf", height = 8, width = 12)
# par(mfrow = c(2,3))
# par(mar = c(23, 5, 4, 2) )
n_rows <- 7
n_cols <- 5
layout(matrix(1:(n_rows * n_cols), nrow = n_rows, byrow = TRUE))
par(mar = c(2, 2, 2, 1))  # Adjust margins for tighter spacing

for(j in 1:33)
{
  tree_final <- getTree(full_tree[[j]], "0")
  print(plot(tree_final, main = unique_subtypes[j]))
  # plot(tree_final, main = unique_subtypes[j])
}

dev.off()

