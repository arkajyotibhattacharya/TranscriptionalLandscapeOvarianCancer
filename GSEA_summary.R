library(data.table)
#hallmark
hallmark = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/GSEA/Hallmark.txt", sep = "\t", header = TRUE))
colnames(hallmark)[3:dim(hallmark)[2]] = paste("TF", c(1:374))
ces_association = read.table("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/Survival\ analysis/CES_association_summary_with_clinical_pathological_variables_and_uni_multivariate_survival_significance_27062020.txt", sep = "\t", header = TRUE)

ces_association$ces_v1 = sapply(as.character(ces_association$CES),function(x){as.numeric(substr(x,2,nchar(x)))})

hallmark = hallmark[,c(1,2,(ces_association$ces_v1+2))]

hallmark = hallmark[,c(1:16)]
rownames(hallmark) = hallmark$`GeneSet identifier`
hallmark = hallmark[,-c(1:2)]

cutoff = -log10(0.01/14)

# hallmark = as.data.frame(t(hallmark))
# 
# hclust_hallmark = hclust(as.dist(1-cor(hallmark)), method = "ward.D2")
# hallmark = hallmark[,hclust_hallmark$order]
# hallmark$TF = rownames(hallmark)
# hallmark = hallmark[,c(dim(hallmark)[2],1:(dim(hallmark)[2]-1))]
# write.table(hallmark, file = "/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/GSEA/Hallmark_before_screening.txt", sep = "\t", row.names = FALSE)


significant_gs_summary = apply(hallmark, 1, function(x){length(which(abs(x)>cutoff))})
significant_gs_summary = significant_gs_summary[which(significant_gs_summary!=0)]
hallmark = hallmark[names(significant_gs_summary),]
hallmark = as.data.frame(t(hallmark))

hclust_hallmark = hclust(as.dist(1-cor(hallmark)), method = "ward.D2")
hallmark = hallmark[,hclust_hallmark$order]

hallmark_v1 = apply(hallmark, 1, function(x){ifelse(x>5, 5, ifelse(x< -5, -5, x))})

hallmark = as.data.frame(t(hallmark_v1))

hallmark$TF = rownames(hallmark)
hallmark = hallmark[,c(dim(hallmark)[2],1:(dim(hallmark)[2]-1))]
write.table(hallmark, file = "/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/GSEA/Hallmark_after_screening.txt", sep = "\t", row.names = FALSE)


#reactome
reactome = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/GSEA/Reactome.txt", sep = "\t", header = TRUE))
colnames(reactome)[3:dim(reactome)[2]] = paste("TF", c(1:374))
ces_association = read.table("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/Survival\ analysis/CES_association_summary_with_clinical_pathological_variables_and_uni_multivariate_survival_significance_27062020.txt", sep = "\t", header = TRUE)

ces_association$ces_v1 = sapply(as.character(ces_association$CES),function(x){as.numeric(substr(x,2,nchar(x)))})

reactome = reactome[,c(1,2,(ces_association$ces_v1+2))]

reactome = reactome[,c(1:16)]
rownames(reactome) = reactome$`GeneSet identifier`
reactome = reactome[,-c(1:2)]

cutoff = -log10(0.01/14)

significant_gs_summary = apply(reactome, 1, function(x){length(which(abs(x)>cutoff))})
significant_gs_summary = significant_gs_summary[which(significant_gs_summary!=0)]
reactome = reactome[names(significant_gs_summary),]
reactome = as.data.frame(t(reactome))

hclust_reactome = hclust(as.dist(1-cor(reactome)), method = "ward.D2")
reactome = reactome[,hclust_reactome$order]

max_reactome = apply(reactome, 1, function(x){max(x)})
reactome_v1 = apply(reactome, 1, function(x){ifelse(x>5, 5, ifelse(x< -5, -5, x))})

reactome = as.data.frame(t(reactome_v1))


 reactome$TF = rownames(reactome)
reactome = reactome[,c(dim(reactome)[2],1:(dim(reactome)[2]-1))]

colnames(reactome) = gsub("REACTOME_", "", colnames(reactome))
colnames(reactome) = gsub("_", " ", colnames(reactome))
colnames(reactome) = tolower(colnames(reactome))
colnames(reactome) = sapply(colnames(reactome), function(x){paste(toupper(substr(x,1,1)), substr(x,2,nchar(x)), sep = "")})
write.table(reactome, file = "/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/GSEA/reactome_after_screening.txt", sep = "\t", row.names = FALSE)




#GO_BP
GO_BP = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/GSEA/GO_BP.txt", sep = "\t", header = TRUE))
colnames(GO_BP)[3:dim(GO_BP)[2]] = paste("TF", c(1:374))
ces_association = read.table("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/Survival\ analysis/CES_association_summary_with_clinical_pathological_variables_and_uni_multivariate_survival_significance_27062020.txt", sep = "\t", header = TRUE)

ces_association$ces_v1 = sapply(as.character(ces_association$CES),function(x){as.numeric(substr(x,2,nchar(x)))})

GO_BP = GO_BP[,c(1,2,(ces_association$ces_v1+2))]

GO_BP = GO_BP[,c(1:16)]
rownames(GO_BP) = GO_BP$`GeneSet identifier`
GO_BP = GO_BP[,-c(1:2)]

cutoff = -log10(0.01/14)

significant_gs_summary = apply(GO_BP, 1, function(x){length(which(abs(x)>cutoff))})
significant_gs_summary = significant_gs_summary[which(significant_gs_summary!=0)]
GO_BP = GO_BP[names(significant_gs_summary),]
GO_BP = as.data.frame(t(GO_BP))

hclust_GO_BP = hclust(as.dist(1-cor(GO_BP)), method = "ward.D2")
GO_BP = GO_BP[,hclust_GO_BP$order]
GO_BP$TF = rownames(GO_BP)
GO_BP = GO_BP[,c(dim(GO_BP)[2],1:(dim(GO_BP)[2]-1))]
write.table(GO_BP, file = "/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/GSEA/GO_BP_after_screening.txt", sep = "\t", row.names = FALSE)

#KEGG
KEGG = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/GSEA/KEGG.txt", sep = "\t", header = TRUE))
colnames(KEGG)[3:dim(KEGG)[2]] = paste("TF", c(1:374))
ces_association = read.table("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/Survival\ analysis/CES_association_summary_with_clinical_pathological_variables_and_uni_multivariate_survival_significance_27062020.txt", sep = "\t", header = TRUE)

ces_association$ces_v1 = sapply(as.character(ces_association$CES),function(x){as.numeric(substr(x,2,nchar(x)))})

KEGG = KEGG[,c(1,2,(ces_association$ces_v1+2))]

KEGG = KEGG[,c(1:16)]
rownames(KEGG) = KEGG$`GeneSet identifier`
KEGG = KEGG[,-c(1:2)]

cutoff = -log10(0.01/14)

significant_gs_summary = apply(KEGG, 1, function(x){length(which(abs(x)>cutoff))})
significant_gs_summary = significant_gs_summary[which(significant_gs_summary!=0)]
KEGG = KEGG[names(significant_gs_summary),]
KEGG = as.data.frame(t(KEGG))

hclust_KEGG = hclust(as.dist(1-cor(KEGG)), method = "ward.D2")
KEGG = KEGG[,hclust_KEGG$order]

KEGG_v1 = apply(KEGG, 1, function(x){ifelse(x>5, 5, ifelse(x< -5, -5, x))})

KEGG = as.data.frame(t(KEGG_v1))

KEGG$TF = rownames(KEGG)
KEGG = KEGG[,c(dim(KEGG)[2],1:(dim(KEGG)[2]-1))]
colnames(KEGG) = gsub("KEGG_", "", colnames(KEGG))
colnames(KEGG) = gsub("_", " ", colnames(KEGG))
colnames(KEGG) = tolower(colnames(KEGG))
colnames(KEGG) = sapply(colnames(KEGG), function(x){paste(toupper(substr(x,1,1)), substr(x,2,nchar(x)), sep = "")})

write.table(KEGG, file = "/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/GSEA/KEGG_after_screening.txt", sep = "\t", row.names = FALSE)
