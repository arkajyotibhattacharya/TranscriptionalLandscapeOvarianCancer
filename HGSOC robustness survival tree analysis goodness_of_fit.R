    library(MST)
      library(parallel)
      library(aods3)
      #Read in any of the subset-multiv_OS_inps ("OS_inp_[series number].txt")
      multiv_OS_inp = read.table( "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Ovarian cancer project/Data/Clinical_pathological_information_of_samples_with_CES_activity_scores_for_multivariate_analysis_flipped.txt", sep = "\t", header = TRUE)
      colnames(multiv_OS_inp) = gsub("\\.","-", colnames(multiv_OS_inp))
      multiv_OS_inp = multiv_OS_inp[which(multiv_OS_inp$Platinum==1),]
      multiv_OS_inp = multiv_OS_inp[which(multiv_OS_inp$Type=="Serous"),]
      multiv_OS_inp = multiv_OS_inp[which(multiv_OS_inp$Grade!=1),]
      multiv_OS_inp$OS_binary = multiv_OS_inp$`OS-binary`
      multiv_OS_inp$Grade <- as.factor(multiv_OS_inp$Grade)
      multiv_OS_inp$Stage <- as.factor(multiv_OS_inp$Stage)
      multiv_OS_inp$Debulking <- as.factor(multiv_OS_inp$Debulking)
      
      
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
      tree_final <- getTree(full_tree, "0")
      
      multiv_OS_inp$term_nodes <- as.factor(predict(tree_final, newdata = multiv_OS_inp, type = 'node'))
      cox_reg = coxph(Surv(OS, OS_binary) ~ term_nodes + cluster(Identifier), data = multiv_OS_inp)
      summary(cox_reg)
      
      