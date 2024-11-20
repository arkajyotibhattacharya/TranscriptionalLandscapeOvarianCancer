#tree
library(MST)
library(parallel)
library(sqldf)
#Read in any of the subset-multiv_OS_inps ("OS_inp_[series number].txt")
multiv_OS_inp = read.table( "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Clinical_pathological_information_of_samples_with_CES_activity_scores_for_multivariate_analysis_flipped.txt", sep = "\t", header = TRUE)
colnames(multiv_OS_inp) = gsub("\\.","-", colnames(multiv_OS_inp))
multiv_OS_inp = multiv_OS_inp[which(multiv_OS_inp$Platinum==1),]
multiv_OS_inp = multiv_OS_inp[which(multiv_OS_inp$Type=="Serous"),]
multiv_OS_inp = multiv_OS_inp[which(multiv_OS_inp$Grade!=1),]

multiv_OS_inp$OS_binary = multiv_OS_inp$`OS-binary`
# multiv_OS_inp$Debulking = as.character(multiv_OS_inp$Debulking)
# multiv_OS_inp$Debulking[which(is.na(multiv_OS_inp$Debulking))] = "unknown"
multiv_OS_inp$Debulking = as.factor(multiv_OS_inp$Debulking)

# multiv_OS_inp$Grade = as.character(multiv_OS_inp$Grade)
# multiv_OS_inp$Grade[which(is.na(multiv_OS_inp$Grade))] = "unknown"
multiv_OS_inp$Grade = as.factor(multiv_OS_inp$Grade)

# multiv_OS_inp$Stage = as.character(multiv_OS_inp$Stage)
# multiv_OS_inp$Stage[which(is.na(multiv_OS_inp$Stage))] = "unknown"
multiv_OS_inp$Stage = as.factor(multiv_OS_inp$Stage)



minsplits = 50

full_tree <- MST(formula = Surv(OS, OS_binary) ~ Debulking + Grade + Stage + Age
                 +V14
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
                 | Identifier,
                 data = multiv_OS_inp,
                 test = multiv_OS_inp,
                 method = "independence",
                 minsplit = minsplits,
                 minevents = ceiling(minsplits/2),
                 minbucket = ceiling(minsplits/3),
                 selection.method = "test.sample",
                 # LeBlanc = TRUE,
                 plot.Ga = TRUE,
                 sortTrees = TRUE,
                 details = FALSE)


make_tree_summary = function(tree)
{
  library(partykit)
  library(stringr)
  
  varlist = c("Debulking",
              "Grade",
              "Stage",
              "Age",
              "V14",
              "V76",
              "V78",
              "V121",
              "V138",
              "V146",
              "V197",
              "V220",
              "V239",
              "V247",
              "V250",
              "V253",
              "V320",
              "V166")
  varid_locations = str_locate_all(pattern ='varid',  as.character(tree)[1])[[1]]
  mean_location = (varid_locations[,1]+varid_locations[,2])/2
  varid_locations = cbind(varid_locations,mean_location)
  id_locations = str_locate_all(pattern ='\\(id = ',  as.character(tree)[1])[[1]]
  mean_location = (id_locations[,1]+id_locations[,2])/2
  id_locations = cbind(id_locations,mean_location)
  
  text_to_search_in = as.character(tree)[1]
  tree_summary = as.data.frame(matrix(NA, dim(varid_locations)[1],3))
  colnames(tree_summary) = c("ID_number", "Var_number", "Var_name")
  for(i in 1:dim(varid_locations)[1])
  {
    tree_summary$ID_number[i] = which(id_locations[,3]>varid_locations[i,3])[1]-1
    tree_summary$Var_number[i] = substr(text_to_search_in, str_locate_all(pattern ='varid = ',  text_to_search_in)[[1]][i,2]+1
                                        ,str_locate_all(pattern ='varid = ',  text_to_search_in)[[1]][i,2]+2)
    tree_summary$Var_number[i] = as.numeric(as.character(gsub(",","",tree_summary$Var_number[i])))
    
  }
  tree_summary$Var_name = varlist[as.numeric(tree_summary$Var_number)-1]
  tab <- sapply(tree_summary$ID_number, function(id) {
    y <- dim(data_party(tree[id]))[1]
    
  })
  
  tree_summary$data_size = tab
  
  tree_summary = sqldf("select Var_name, max(data_size) as data_size from tree_summary group by Var_name")
  rownames(tree_summary) = tree_summary$Var_name
  tree_summary_for_non_signficant_varlist = as.data.frame(matrix(NA, length(varlist) - dim(tree_summary)[1] ,2))
  colnames(tree_summary_for_non_signficant_varlist) = c("Var_name", "data_size")
  rownames(tree_summary_for_non_signficant_varlist) = varlist[which(!varlist%in%rownames(tree_summary))]
  tree_summary_for_non_signficant_varlist$Var_name = rownames(tree_summary_for_non_signficant_varlist) 
  tree_summary_for_non_signficant_varlist$data_size = 0
  tree_summary = rbind(tree_summary, tree_summary_for_non_signficant_varlist )
  return(tree_summary)
}

tree_final <- getTree(full_tree, "0")
summary_file_main_tree = make_tree_summary(tree_final)

all_tree_summary  = list()

rank_correlation_data_size_per_classifier_between_main_and_current_tree = array(NA, 20000)
try2 <- function(code, silent = FALSE) {
  tryCatch(code, error = function(c) {
    if (!silent) {-1}
    else{code}})}


for(i in 1:10)
{
  time1 = proc.time()[3]
  
  load(paste("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/Survival analysis/trees_multivariate_with_significant_CESs_HGSOC_iterations_",i,".RData", sep = ""))
  
  extract_cor = function(j){
    current_data = make_tree_summary(fit_results[[j]]$tree0)
    
    common_var_name = intersect(current_data$Var_name, summary_file_main_tree$Var_name)
    
    summary_file_main_tree_temp = summary_file_main_tree[common_var_name,]
    current_data = current_data[common_var_name,]
    
    return(try2(cor(current_data$data_size, summary_file_main_tree_temp$data_size, method = c("spearman"), use = "pairwise.complete.obs")))
    
  }
  
  cl = makeCluster(10, type = "FORK")
  rank_correlation_data_size_per_classifier_between_main_and_current_tree[((i-1)*2000+1):((i-1)*2000+2000)] <- unlist(parLapply(cl,1:2000,extract_cor))
  stopCluster(cl)
  
  
  rm(fit_results)
  
  # for(j in 1:2000)
  # {
  #   current_data = make_tree_summary(fit_results[[(i-1)*2000+j]]$tree0)
  #   
  #   common_var_name = intersect(current_data, summary_file_main_tree$Var_name)
  #   
  #   summary_file_main_tree_temp = summary_file_main_tree[common_var_name,]
  #   current_data = current_data[common_var_name,]
  # 
  #   rank_correlation_data_size_per_classifier_between_main_and_current_tree[(i-1)*2000+j] = try2(cor(current_data$data_size, summary_file_main_tree_temp$data_size, method = c("spearman"), use = "pairwise.complete.obs"))
  #   
  # }
  print((proc.time()-time1)/60)
  print(i)
}

summary(rank_correlation_data_size_per_classifier_between_main_and_current_tree)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.4480  0.3527  0.5219  0.5095  0.6849  0.9975 
# > quantile(rank_correlation_data_size_per_classifier_between_main_and_current_tree, na.rm = T)
# 0%        25%        50%        75%       100% 
# -1.0000000  0.1428571  0.6571429  0.9000000  1.0000000 

make_tree_summary = function(tree)
{
  library(stringr)
  library(partykit)
  
  varlist = c("Debulking",
              "Grade",
              "Stage",
              "Age",
              "V14",
              "V76",
              "V78",
              "V121",
              "V138",
              "V146",
              "V197",
              "V220",
              "V239",
              "V247",
              "V250",
              "V253",
              "V320",
              "V166")
  varid_locations = str_locate_all(pattern ='varid',  as.character(tree)[1])[[1]]
  mean_location = (varid_locations[,1]+varid_locations[,2])/2
  varid_locations = cbind(varid_locations,mean_location)
  id_locations = str_locate_all(pattern ='\\(id = ',  as.character(tree)[1])[[1]]
  mean_location = (id_locations[,1]+id_locations[,2])/2
  id_locations = cbind(id_locations,mean_location)
  
  text_to_search_in = as.character(tree)[1]
  tree_summary = as.data.frame(matrix(NA, dim(varid_locations)[1],3))
  colnames(tree_summary) = c("ID_number", "Var_number", "Var_name")
  for(i in 1:dim(varid_locations)[1])
  {
    tree_summary$ID_number[i] = which(id_locations[,3]>varid_locations[i,3])[1]-1
    tree_summary$Var_number[i] = substr(text_to_search_in, str_locate_all(pattern ='varid = ',  text_to_search_in)[[1]][i,2]+1
                                        ,str_locate_all(pattern ='varid = ',  text_to_search_in)[[1]][i,2]+2)
    tree_summary$Var_number[i] = as.numeric(as.character(gsub(",","",tree_summary$Var_number[i])))
    
  }
  tree_summary$Var_name = varlist[as.numeric(tree_summary$Var_number)-1]
  tab <- sapply(tree_summary$ID_number, function(id) {
    y <- dim(data_party(tree[id]))[1]
    
  })
  
  tree_summary$data_size = tab
  
  tree_summary = sqldf("select Var_name, max(data_size) as data_size from tree_summary group by Var_name")
  
  rownames(tree_summary) = tree_summary$Var_name
  
  tree_summary_for_non_signficant_varlist = as.data.frame(matrix(NA, length(varlist) - dim(tree_summary)[1] ,2))
  colnames(tree_summary_for_non_signficant_varlist) = c("Var_name", "data_size")
  rownames(tree_summary_for_non_signficant_varlist) = varlist[which(!varlist%in%rownames(tree_summary))]
  tree_summary_for_non_signficant_varlist$Var_name = rownames(tree_summary_for_non_signficant_varlist) 
  tree_summary_for_non_signficant_varlist$data_size = 0
  tree_summary = rbind(tree_summary, tree_summary_for_non_signficant_varlist )
  tree_summary$rank = rank(-tree_summary$data_size, ties.method = "min")
  
  return(tree_summary)
}


all_tree_summary_in_matrix = NULL

for(i in 1:10)
{
  time1 = proc.time()[3]
  
  load(paste("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/Survival analysis/trees_multivariate_with_significant_CESs_HGSOC_iterations_",i,".RData", sep = ""))
  
  for(j in 1:2000)
  {
    temp_mat = make_tree_summary(fit_results[[j]]$tree0)
    temp_mat$iteration = (i-1)*2000+j
    all_tree_summary_in_matrix = rbind(all_tree_summary_in_matrix, temp_mat)
  }
  print((proc.time()-time1)/60)
  print(i)
  rm(fit_results)
}

for(i in 1:max(all_tree_summary_in_matrix$rank))
{
  print(i)
  print(sort(table(all_tree_summary_in_matrix$Var_name[which(all_tree_summary_in_matrix$rank==i)]), decreasing = TRUE))
}

write.table(all_tree_summary_in_matrix, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Results/Survival\ analysis/Robustness_details_all_iterations_HGSOC.txt", sep = "\t", row.names = FALSE)
###############################################################

all_tree_summary  = list()
original_tree_summary = make_tree_summary(full_tree$tree0)


for(i in 1:10)
  {
  time1 = proc.time()[3]
  
  rm(fit_results)
  load(paste("/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/Survival\ analysis/trees_multivariate_with_significant_CESs_iterations_",i,".RData", sep = ""))
  
  for(j in 1:2000)
  {
    all_tree_summary[[(i-1)*2000+j]] = make_tree_summary(fit_results[[(i-1)*2000+j]]$tree0)
  }
  print((proc.time()-time1)/60)
  print(i)
  }

all_tree_summary_in_matrix = NULL
for(i in 1:20000)
{
  temp_mat = all_tree_summary[[i]]
  temp_mat$iteration = i
  
  all_tree_summary_in_matrix = rbind(all_tree_summary_in_matrix, temp_mat)
  
}

length(which(all_tree_summary_in_matrix$ID_number[which(all_tree_summary_in_matrix$Var_name=="V121")]==1))

table(all_tree_summary_in_matrix$Var_name[which(all_tree_summary_in_matrix$ID_number==1)])

all_tree_summary_in_matrix_v121 = all_tree_summary_in_matrix[intersect(which(all_tree_summary_in_matrix$ID_number==1)
                                                                       ,which(all_tree_summary_in_matrix$Var_name=="V121")),]
all_tree_summary_in_matrix_v121 = all_tree_summary_in_matrix[which(all_tree_summary_in_matrix$iteration%in%all_tree_summary_in_matrix_v121$iteration),]


length(which(all_tree_summary_in_matrix_v121$ID_number[which(all_tree_summary_in_matrix_v121$Var_name=="Stage")]==2))

table(all_tree_summary_in_matrix$Var_name[which(all_tree_summary_in_matrix$ID_number==1)])
all_tree_summary_in_matrix_to_store = all_tree_summary_in_matrix


all_tree_summary_in_matrix = all_tree_summary_in_matrix_to_store
for(i in 1:dim(original_tree_summary)[1])
{
  print(original_tree_summary$Var_name[i])
  print(length(which(all_tree_summary_in_matrix$ID_number[which(all_tree_summary_in_matrix$Var_name==original_tree_summary$Var_name[i])]==original_tree_summary$ID_number[i]))/length(unique(all_tree_summary_in_matrix$iteration)))
  print(table(all_tree_summary_in_matrix$Var_name[which(all_tree_summary_in_matrix$ID_number==original_tree_summary$ID_number[i])]))
  temp_data = all_tree_summary_in_matrix[intersect(which(all_tree_summary_in_matrix$ID_number==original_tree_summary$ID_number[i])
                                                                         ,which(all_tree_summary_in_matrix$Var_name==original_tree_summary$Var_name[i])),]
  
  all_tree_summary_in_matrix = all_tree_summary_in_matrix[which(all_tree_summary_in_matrix$iteration%in%temp_data$iteration),]
  
}

save(all_tree_summary_in_matrix_to_store, file = "/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/Survival\ analysis/trees_multivariate_with_significant_CESs_iterations_all.RData")


#####################################################################pipeline###################################################################
make_tree_summary = function(tree, varlist)
{
  library(stringr)
  
  
  varid_locations = str_locate_all(pattern ='varid',  as.character(tree)[1])[[1]]
  mean_location = (varid_locations[,1]+varid_locations[,2])/2
  varid_locations = cbind(varid_locations,mean_location)
  id_locations = str_locate_all(pattern ='\\(id = ',  as.character(tree)[1])[[1]]
  mean_location = (id_locations[,1]+id_locations[,2])/2
  id_locations = cbind(id_locations,mean_location)
  
  text_to_search_in = as.character(tree)[1]
  tree_summary = as.data.frame(matrix(NA, dim(varid_locations)[1],3))
  colnames(tree_summary) = c("ID_number", "Var_number", "Var_name")
  for(i in 1:dim(varid_locations)[1])
  {
    tree_summary$ID_number[i] = which(id_locations[,3]>varid_locations[i,3])[1]-1
    tree_summary$Var_number[i] = substr(text_to_search_in, str_locate_all(pattern ='varid = ',  text_to_search_in)[[1]][i,2]+1
                                        ,str_locate_all(pattern ='varid = ',  text_to_search_in)[[1]][i,2]+2)
    tree_summary$Var_number[i] = as.numeric(as.character(gsub(",","",tree_summary$Var_number[i])))
    
  }
  tree_summary$Var_name = varlist[as.numeric(tree_summary$Var_number)-1]
  
  return(tree_summary)
}

tree_robustness_scoring = function(original_tree = full_tree$tree0
                                   , RDataname_with_path_without_iteration_number = "/Users/arkajyotibhattacharya/Projects/Ovarian\ cancer\ project/Results/Survival\ analysis/trees_multivariate_with_significant_CESs_iterations_"
                                   , number_of_RData_files = 10
                                   , number_of_iterations_in_each_RData = 2000
                                   , varlist = c("Debulking",
                                                 "Grade",
                                                 "Stage",
                                                 "Age",
                                                 "V14",
                                                 "V76",
                                                 "V78",
                                                 "V121",
                                                 "V138",
                                                 "V146",
                                                 "V197",
                                                 "V220",
                                                 "V239",
                                                 "V247",
                                                 "V250",
                                                 "V253",
                                                 "V320",
                                                 "V166"))
{
  all_tree_summary  = list()
  original_tree_summary = make_tree_summary(original_tree, varlist)
  
  for(i in 1:number_of_RData_files)
  {
    time1 = proc.time()[3]
    
    load(paste(RDataname_with_path_without_iteration_number,i,".RData", sep = ""))
    
    for(j in 1:number_of_iterations_in_each_RData)
    {
      all_tree_summary[[(i-1)*number_of_iterations_in_each_RData+j]] = make_tree_summary(fit_results[[(i-1)*number_of_iterations_in_each_RData+j]]$tree0, varlist)
    }
    rm(fit_results)
    
    print((proc.time()-time1)/60)
    print(i)
  }
  
  all_tree_summary_in_matrix = NULL
  for(i in 1:(number_of_iterations_in_each_RData*number_of_RData_files))
  {
    temp_mat = all_tree_summary[[i]]
    temp_mat$iteration = i
    
    all_tree_summary_in_matrix = rbind(all_tree_summary_in_matrix, temp_mat)
    
  }

  
  
  all_tree_summary_in_matrix_to_store = all_tree_summary_in_matrix
  
  robustness_scores = list()
  for(i in 1:dim(original_tree_summary)[1])
  {
    robustness_scores[[i]] = list()
    robustness_scores[[i]]$Var_name = original_tree_summary$Var_name[i]
    robustness_scores[[i]]$occurrence_percentage = length(which(all_tree_summary_in_matrix$ID_number[which(all_tree_summary_in_matrix$Var_name==original_tree_summary$Var_name[i])]==original_tree_summary$ID_number[i]))/length(unique(all_tree_summary_in_matrix$iteration))
    robustness_scores[[i]]$occurrence_distributions = table(all_tree_summary_in_matrix$Var_name[which(all_tree_summary_in_matrix$ID_number==original_tree_summary$ID_number[i])])
    
    # print(original_tree_summary$Var_name[i])
    # print(length(which(all_tree_summary_in_matrix$ID_number[which(all_tree_summary_in_matrix$Var_name==original_tree_summary$Var_name[i])]==original_tree_summary$ID_number[i]))/length(unique(all_tree_summary_in_matrix$iteration)))
    # print(table(all_tree_summary_in_matrix$Var_name[which(all_tree_summary_in_matrix$ID_number==original_tree_summary$ID_number[i])]))
    temp_data = all_tree_summary_in_matrix[intersect(which(all_tree_summary_in_matrix$ID_number==original_tree_summary$ID_number[i])
                                                     ,which(all_tree_summary_in_matrix$Var_name==original_tree_summary$Var_name[i])),]
    
    all_tree_summary_in_matrix = all_tree_summary_in_matrix[which(all_tree_summary_in_matrix$iteration%in%temp_data$iteration),]
    
  }
  
  temp_mat = original_tree_summary
  temp_mat$iteration = 0
  
  all_tree_summary_in_matrix_to_store = rbind(all_tree_summary_in_matrix_to_store, temp_mat)
  
  save(all_tree_summary_in_matrix_to_store, file = paste(RDataname_with_path_without_iteration_number,"_all_tree_summary.RData", sep = ""))
  return(robustness_scores)
}


tree_robustness_scoring()

