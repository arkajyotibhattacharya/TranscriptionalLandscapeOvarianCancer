    library(MST)
      library(parallel)
      
      #Read in any of the subset-multiv_OS_inps ("OS_inp_[series number].txt")
      multiv_OS_inp = read.table( "F:/My Drive/Projects_27092022/Ovarian cancer project/Data/Clinical_pathological_information_of_samples_with_CES_activity_scores_for_multivariate_analysis_flipped.txt", sep = "\t", header = TRUE)
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
      
      (tree_final <- getTree(full_tree, "0"))
      plot(tree_final)
      
      samples_used_to_build_tree = list()
      fit_results = list()
      
      time2 = proc.time()[3]
      
      for(iterations in 1:2000)
      {
        
        time1 = proc.time()[3]
        
        minsplits = ceiling(50)
        set.seed(42345+iterations)
        samples_used_to_build_tree[[iterations]] = sample(1:nrow(multiv_OS_inp),floor(nrow(multiv_OS_inp)*0.8), replace = FALSE)
        data_to_build_tree = multiv_OS_inp[samples_used_to_build_tree[[iterations]],]
        
        fit_results[[iterations]] = MST(formula = Surv(OS, OS_binary) ~ Debulking + Grade + Stage + Age
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
                                        data = data_to_build_tree,
                                        test = data_to_build_tree,
                                        method = "independence",
                                        minsplit = minsplits,
                                        minevents = ceiling(minsplits/2),
                                        minbucket = ceiling(minsplits/3),
                                        selection.method = "test.sample",
                                        # LeBlanc = TRUE,
                                        plot.Ga = TRUE,
                                        sortTrees = TRUE,
                                        details = FALSE)
        print((proc.time()[3] - time1)/60)
        print(paste("iteration", iterations))
        
      }
      print((proc.time()[3] - time2)/60)
      print(paste("Full"))
      save(fit_results, file = "F:/My Drive/Projects_27092022/Ovarian cancer project/Results/Survival\ analysis/trees_multivariate_with_significant_CESs_HGSOC_iterations_4.RData")
      
      
      
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
      
      (tree_final <- getTree(full_tree, "0"))
      pdf(paste(" F:/My Drive/Projects_27092022/Ovarian cancer project/Plots/Survival\ related\ plots/Tree_multivariate_HGSOC_full_tree_iteration_4.pdf", sep = ""), height = 10, width = 20)
      print(plot(tree_final))
      dev.off()
      
      save.image(" F:/My Drive/Projects_27092022/Ovarian cancer project/Results/Survival\ analysis/trees_multivariate_with_HGSOC_significant_CESs_iteration_4.RData")
      