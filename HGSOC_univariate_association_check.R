
library(MST)
library(parallel)
library(aods3)
#Read in any of the subset-multiv_OS_inps ("OS_inp_[series number].txt")
multiv_OS_inp = read.table( "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Clinical_pathological_information_of_samples_with_CES_activity_scores_for_multivariate_analysis_flipped.txt", sep = "\t", header = TRUE)
colnames(multiv_OS_inp) = gsub("\\.","-", colnames(multiv_OS_inp))
# multiv_OS_inp = multiv_OS_inp[which(multiv_OS_inp$Platinum==1),]
multiv_OS_inp = multiv_OS_inp[which(multiv_OS_inp$Type=="Serous"),]
multiv_OS_inp = multiv_OS_inp[which(multiv_OS_inp$Grade!=1),]
multiv_OS_inp$OS_binary = multiv_OS_inp$`OS-binary`
multiv_OS_inp$Grade <- as.factor(multiv_OS_inp$Grade)
multiv_OS_inp$Stage <- as.factor(multiv_OS_inp$Stage)
multiv_OS_inp$Debulking <- as.factor(multiv_OS_inp$Debulking)


data = multiv_OS_inp
follow_up= data$OS
event = as.numeric(data$`OS-binary`)
outcome_variable_range=c(21:394)
list_of_predictors=names(data[outcome_variable_range])
newdata = data[,outcome_variable_range]						
FDR = 0.01
nPerm = 10000
conf_level = 0.8
multiv_OS_inp$Series_original <- as.factor(multiv_OS_inp$Series_original)

# Coxph loop
result <- lapply(newdata, function(x){coxph(as.formula(paste("Surv(follow_up,event)~",paste("x"))))})

# Create matrix for results  
full_coef_matrix = NULL
for(i in names(result)){
  each_result_summary = summary(result[[i]])
  full_coef_matrix = rbind(full_coef_matrix, cbind(each_result_summary$coefficients,each_result_summary$conf.int)[1,])
  rownames(full_coef_matrix)[dim(full_coef_matrix)[1]] = i}



# Repeat loop for nPerm permutations								                                          #like this the loops works OK for univariate analysis
time1=Sys.time()
perm_data=NULL										                                                          #Removed "if" statement, will always run mp_testing
options(warn=1)                                                                             #essential for withCallingHandlers to catch warnings/permutation
warn_total = NULL
# for (i in 1:nPerm){
multiple_testing_survival = function(i){
  
  adjusted_Pval = NULL
  time2=Sys.time()
  set.seed(12345+i)										                                                    #makes sure the RNG in "sample" does exactly thesame when rerunning the code (so, purely for reproducibility)
  perm_data = newdata[sample(nrow(newdata)),]			                                        #creates permutation of the data, works
  perm_result = NULL
  warn_perm = NULL
  
  # withCallingHandlers({
    for(ces in list_of_predictors){
      perm_result[[ces]]$survival = coxph(Surv(follow_up,event)~ perm_data[[ces]])
      # assign("last.warning", NULL, envir = baseenv())
    }
  # }, warning = function(w){
  #   warn_perm <<- c(ces,i, w$message)
  #   warn_total <<- rbind(warn_total, warn_perm)
  # })
  
  for (j in names(perm_result)){								
    each_perm_summary = summary(perm_result[[j]]$survival)
    adjusted_Pval = rbind(adjusted_Pval, each_perm_summary$coefficients[1,5])
    rownames(adjusted_Pval)[dim(adjusted_Pval)[1]] = j}
  # adjusted_Pval_combined = cbind(adjusted_Pval_combined,adjusted_Pval)
  output = list()
  output$adjusted_Pval = adjusted_Pval
  # output$warn_total = warn_total
  return(output)
  print(c("Time cox permutation round",i,":",Sys.time()-time2))
}

time1 = proc.time()[3]
no_cores = 10
cl <- makeCluster(no_cores, type = "FORK")

outputs_combined = parLapply(cl, 1:nPerm,multiple_testing_survival)
# outputs_combined = parLapply(cl, 1:20,multiple_testing_survival)
stopCluster(cl)
print((proc.time()[3] - time1)/60)

adjusted_Pval_combined = matrix(NA,dim(outputs_combined[[1]]$adjusted_Pval)[1], nPerm)

for(i in 1:nPerm)
{
  adjusted_Pval_combined[,i] = outputs_combined[[i]]$adjusted_Pval
  
}


print(c("Time cox permutations total:", Sys.time()-time1))

# Calculate MVP threshold
time1=Sys.time()
full_coef_matrix = as.data.frame(full_coef_matrix)
full_coef_matrix$mlogpva = -log10(full_coef_matrix$`Pr(>|z|)`)				          #works correctly

original_mlogpval_ordered = sort(full_coef_matrix$mlogpva,decreasing=TRUE )			#becomes a sorted list, no longer in the DF
adjusted_Pval_combined_mlogpval = -log10(adjusted_Pval_combined)
adjusted_Pval_combined_mlogpval_v1 = apply(abs(adjusted_Pval_combined_mlogpval),2,function(x) sort(x,decreasing=TRUE))      #Should remain a DF in this step. Watch out, any problems in the earlier steps (eg. warnings) may introduce NA's, forcing this into layered list that cannot be processed in the subsequent script.

cutoffs = array(0,nPerm)									#creates a list including all cutoff results
for(j in 1:nPerm)
{
  for(row_number in 1:dim(adjusted_Pval_combined_mlogpval_v1)[1])
  {
    if(length(which(adjusted_Pval_combined_mlogpval_v1[,j]>original_mlogpval_ordered[row_number]))/row_number>=FDR)
    {
      cutoffs[j] = original_mlogpval_ordered[row_number]				#Corrected: now takes p-value from original p-values instead of permutation-p-values
      break
    }
  }
}												

#Calculate cutoff based on confidence level and use this to create a results matix with clear FDR TRUE/FALSE based on cutoff
list_of_cutoffs = quantile(cutoffs, conf_level)						#uses this list (in combination with the confidence level to mark a single cutoff value), works as expected
output = ifelse(full_coef_matrix$mlogpva<list_of_cutoffs,"False","True")			#checks every -log10pvalue in the original data for the cutoff, only values more extreme should become "TRUE"
results_final = cbind(full_coef_matrix, output)
colnames(results_final)[11]=c(paste("MVP_reject_H0", "nPerm=", nPerm, "FDR=", FDR, "conf=", conf_level))

#Save data after analysis to prevent loss during startup PC/R next day
fileName_p = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian\ cancer\ project/Results/Survival analysis/Univariate_summary_14052024.txt"
fileName_warn = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian\ cancer\ project/Results/Survival analysis/Univariate_warning_summary_14052024.txt"
fileName_p_perms = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian\ cancer\ project/Results/Survival analysis/Univariate_permutation_pvalues_14052024.txt"

write.table(results_final, file = fileName_p, quote = FALSE, sep = "\t")
write.table(warn_total, file= fileName_warn, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(adjusted_Pval_combined, file= fileName_p_perms, quote = FALSE, sep = "\t", row.names = FALSE)  

# Information on the coxph loop
print(c("Total:", each_result_summary$n))
print(c("Events:", each_result_summary$nevent))
print(c("Censored:", (each_result_summary$n)-(each_result_summary$nevent)))
print(warnings())


