library(sqldf)

total_time = proc.time()[3]
time1 = proc.time()[3]
setwd("/data/bioinfo-fehrmann/Ovarian Cancer project with Thijs/Results/ICA/")
expression_data = read.table("/data/bioinfo-fehrmann/Ovarian Cancer project with Thijs/Data/Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_12112018.txt", header = TRUE)
data_with_chr_bp_mapping = read.csv("/data/bioinfo-fehrmann/Ovarian Cancer project with Thijs/Data/Genomic_Mapping_hgu133plus2_using_jetscore_30032018.csv")

print(paste("time taken to load the data is", round((proc.time()[3]-time1)/60,2),"mins"))
source("/data/bioinfo-fehrmann/Ovarian Cancer project with Thijs/Codes/combat_consensus_ica_clustering_parallel_20170208.R")


combat_consensus_ica_clustering(expression_data = expression_data,probelevel_standardization = TRUE,
                                combat_identifier = FALSE,
                                Title="14112018_Ovarian_cancer_GEO_TS_QCed_Duplicate_removed",send_email_indicator = FALSE,
                                email_id="arkajyoti.bhattacharya@gmail.com",
                                IC_clustering_algo1 = TRUE,
                                batch_correction_check =FALSE,n_princomp_check = TRUE,
                                prin_comp_cor = FALSE,choose_fastICA_ncomp = FALSE,
                                var_cutoff_to_choose_ncomp = 0.90,ncomp_without_pca=846,
                                nsample_fastica = 25,
                                seed = 12345678,fastica_alg_typ = "parallel",
                                fastica_fun = "logcosh", fastica_alpha = 1,fastica_method = "R",
                                fastica_row_norm = FALSE, fastica_maxit = 1000,
                                fastica_tol = 0.0001, fastica_verbose = FALSE,ICA_Clustering = TRUE,
                                consensus_clustering_of_data = FALSE,
                                distance_function = 'pearson',no_mad_ge = 5000,
                                consensus_maxK=20,consensus_reps=100,consensus_pItem=0.8,
                                consensus_pFeature=1,consensus_clusterAlg="hc",
                                permutation_testing_to_find_state_amp_del_indicator = FALSE,
                                data_with_chr_bp_mapping=data_with_chr_bp_mapping,
                                intervals_for_gaus_kernel = seq(10000,2000000,by=10000),
                                probe_no_for_gaus_kernel = 10,
                                FDR = c(0.01,0.05,0.1,0.2),
                                CL = 0.5,
                                state_deciding_cutoff = 0.95,
                                min_probesets = 5,
                                set_seed = 12345678,
                                gene_level_mapping_indicator = TRUE, ICA_profiling = FALSE, no_cores_v1 = 1)



print(paste("time taken to do the whole analysis is", round((proc.time()[3]-total_time)/60,2),"mins"))

