library(sqldf)
library(truncnorm)

ces_data = read.table("/data/bioinfo-fehrmann/Ovarian Cancer project with Thijs/Results/ICA/ICA_0.9_Explained_variance_25_iterations_Probelevel_standardized__parallel/ICA_Consensus_results/Genelevel_using_jetset_Consensus_Independent_Components_14112018_Ovarian_cancer_GEO_TS_QCed_Duplicate_removed_after_flip.txt", sep = "\t", header = TRUE)

ces_data = ces_data[which(!(ces_data$CHR_Mapping%in%c("X","Y"))),]

ces_data$CHR_Mapping = as.numeric(as.character(ces_data$CHR_Mapping))

ces_data = sqldf("select * from ces_data order by CHR_Mapping, BP_Mapping")
rownames(ces_data) = ces_data$ENTREZID
mapping_file = ces_data[,c(2:12)]
ces_data = ces_data[,c(13:dim(ces_data)[2])]




dataset=as.matrix(ces_data)
data_with_chr_bp_mapping = mapping_file
intervals_for_gaus_kernel = seq(10000,1000000,by=10000)
probe_no_for_gaus_kernel = 10
Title = "DEGR"
FDR = 0.05
CL = 0.5
state_deciding_cutoff = 0.85
min_probesets=5
set_seed = 123456
original_mapping_file = data_with_chr_bp_mapping


full_time = proc.time()[3]
row_num = dim(dataset)[1]
# if(row_num!=dim(data_with_chr_bp_mapping)[1])
# {
#   print("two datasets are not of equal rows")
# }else{
#   if(is.null(colnames(dataset)))
#   {
#     colnames(dataset) = paste("V",1:dim(dataset)[2],sep ="")
#   }
  ##########################
  #Adjustment for Plot
  ###########################
  
  file_SCNA = file.path("/data/bioinfo-fehrmann/Ovarian Cancer project with Thijs/Results/ICA/ICA_0.9_Explained_variance_25_iterations_Probelevel_standardized__parallel/", paste("CNA",
                                         "FDR",FDR,
                                         "CL", CL, "state_deciding_cutoff",state_deciding_cutoff,
                                         "probe_no_for_gaus_kernel",probe_no_for_gaus_kernel,sep = "_"))
  dir.create(file_SCNA, showWarnings = FALSE)
  setwd(file_SCNA)
  
  data_with_chr_bp_mapping$CHR_Mapping = as.numeric(as.character(data_with_chr_bp_mapping$CHR_Mapping))
  chromosome_seq = c(0,cumsum(table(as.numeric(as.character(data_with_chr_bp_mapping$CHR_Mapping)))))
  
  data_with_chr_bp_mapping$BP_Mapping_v1 = data_with_chr_bp_mapping$BP_Mapping/10
  
  for(i in 2:max(data_with_chr_bp_mapping$CHR_Mapping))
  {
    data_with_chr_bp_mapping$BP_Mapping_v1[(chromosome_seq[i]+1):chromosome_seq[i+1]]=data_with_chr_bp_mapping$BP_Mapping_v1[(chromosome_seq[i]+1):chromosome_seq[i+1]]+data_with_chr_bp_mapping$BP_Mapping_v1[chromosome_seq[i]]
  }
  
  label_pos = sqldf("select distinct CHR_Mapping, avg(BP_Mapping_v1) as BP_Mapping from data_with_chr_bp_mapping group by 1 order by 1")
  
  
  
  ########################
  #Finding proper interval for sliding Gaussian Kernel
  #Choosing that interval where 10 or more no. of probesets 
  #corresponding to that chromosome are there in +/- 3*interval for 95% of the cases
  ########################
  quantile_5 = list()
  for( k in 1:max(data_with_chr_bp_mapping$CHR_Mapping))
  {
    quantile_5[[k]] = array(0,length(intervals_for_gaus_kernel))
    for(i in 1:length(intervals_for_gaus_kernel))
    {
      frequency_of_neigh_probes = NULL
      rows = list()
      for(j in (chromosome_seq[k]+1):(chromosome_seq[k+1]))
      {
        frequency_of_neigh_probes[j-chromosome_seq[k]] = length(which(abs(data_with_chr_bp_mapping$BP_Mapping[j] - data_with_chr_bp_mapping$BP_Mapping[c((chromosome_seq[k]+1):(chromosome_seq[k+1]))]) < 3*intervals_for_gaus_kernel[i]))
        
      }
      quantile_5[[k]][i] = quantile(frequency_of_neigh_probes,0.05)%/%1
      
      if(quantile_5[[k]][i]>=probe_no_for_gaus_kernel){
        print(paste("found",probe_no_for_gaus_kernel, "or more number of probe sets in 95% of the chromosome at interval", 
                    intervals_for_gaus_kernel[i],"for Chromosome",k))
        quantile_5[[k]][which(quantile_5[[k]]==0)] = NA
        break
      } 
    }
    
  }
  
  
  ###########################
  #Creating the Density Matrix to get sliding gaussian Kernel
  ###########################
  
  density_matrix = matrix(0, row_num, row_num)
  
  for(j in 1:row_num)
  {
    rows[[j]] = which(data_with_chr_bp_mapping$CHR_Mapping==data_with_chr_bp_mapping$CHR_Mapping[j])
    k = data_with_chr_bp_mapping$CHR_Mapping[j]
    s_x = as.numeric(data_with_chr_bp_mapping$BP_Mapping[rows[[j]]])
    s_mean = as.numeric(data_with_chr_bp_mapping$BP_Mapping[j]) 
    
    den_tnorm = dtruncnorm(s_x,
                           a = s_x[1],
                           b = s_x[length(s_x)],
                           mean = s_mean, 
                           sd = intervals_for_gaus_kernel[max(which(!is.na(quantile_5[[k]])))])
    
    density_matrix[j,rows[[j]]] = den_tnorm/sum(den_tnorm)
  }
  
  
  ###########################
  #Assigning intervals to a vector for each chromosome
  ###########################
  
  final_intervals = NULL
  
  for(chr_no in 1:max(data_with_chr_bp_mapping$CHR_Mapping)){
    final_intervals[chr_no] = intervals_for_gaus_kernel[max(which(!is.na(quantile_5[[chr_no]])))]
  }
  
  
  ###########################
  #Permutation Test for obtaining Amplified or Deleted genomic positions
  ###########################
  
  ################
  #Creating different directories for different outputs
  ################
  
  file_genelevel_bp = file.path(getwd(), "SCNA_sorted_by_base_pair_number_genelevel")
  file_genelevel_value = file.path(getwd(), "SCNA_sorted_by_component_value_genelevel")
  dir.create(file_genelevel_bp, showWarnings = FALSE)
  dir.create(file_genelevel_value, showWarnings = FALSE)  
  
  
  
  ampl_del_probesets_summary = NULL
  amplified_or_deleted_section_summary = NULL
  for(IC_no in 1:dim(dataset)[2])
  {
    set.seed(set_seed+IC_no)
    time1 = proc.time()[3]
    N <- cbind(dataset[,IC_no],replicate(1000, sample(dataset[,IC_no])))
    
    permute_1 = density_matrix%*%N
    
    permute_1_v1 = apply(abs(permute_1),2,function(x) sort(x,decreasing=TRUE))
    
    
    cutoffs = array(0,1000)
    for(j in 2:1001)
    {
      for(i in 1:dim(permute_1_v1)[1])
      {
        if(length(which(permute_1_v1[,j]>permute_1_v1[i,1]))/i>=FDR)
        {
          cutoffs[j-1] = permute_1_v1[i,1]
          break
        }
      }
    }
    
    
    
    indicator_ica_2 = ifelse(permute_1[,1]>quantile(cutoffs,CL),1,ifelse(permute_1[,1]< -quantile(cutoffs,CL),-1,0))
    indicator_ica_3 = density_matrix%*%indicator_ica_2
    indicator_ica_4 = ifelse(indicator_ica_3>state_deciding_cutoff,1, ifelse(indicator_ica_3< -state_deciding_cutoff,-1,0))
    
    
    if(min_probesets>0)
    {
      amp_or_del_pos = which(abs(indicator_ica_4)>0)
      
      for(no in 1:length(amp_or_del_pos)){
        if(length(intersect(which(data_with_chr_bp_mapping$BP_Mapping[data_with_chr_bp_mapping$CHR_Mapping==data_with_chr_bp_mapping$CHR_Mapping[amp_or_del_pos[no]]] >= data_with_chr_bp_mapping$BP_Mapping[amp_or_del_pos[no]]-3*final_intervals[data_with_chr_bp_mapping$CHR_Mapping[amp_or_del_pos[no]]]),
                            which(data_with_chr_bp_mapping$BP_Mapping[data_with_chr_bp_mapping$CHR_Mapping==data_with_chr_bp_mapping$CHR_Mapping[amp_or_del_pos[no]]] <= data_with_chr_bp_mapping$BP_Mapping[amp_or_del_pos[no]]+3*final_intervals[data_with_chr_bp_mapping$CHR_Mapping[amp_or_del_pos[no]]])))<min_probesets){
          indicator_ica_4[amp_or_del_pos[no]] = 0
        }
      }
      
    }
    indicator_run = rle(as.vector(indicator_ica_4)) 
    
    limit = max(abs(min(dataset[,IC_no])),abs(max(dataset[,IC_no]))) 
    
    
    
    
    
    ######################
    #getting the amplified/deleted genomic positions
    ######################
    which_probe_amp_del = which(abs(indicator_ica_4)>0)
    if(length(which_probe_amp_del)>0){
      ampl_del_probesets = as.data.frame(cbind(
        dataset[which_probe_amp_del,IC_no],
        row.names(dataset)[which_probe_amp_del],
        unlist(lapply(data_with_chr_bp_mapping$GENETITLE[which_probe_amp_del],as.character)),
        unlist(lapply(data_with_chr_bp_mapping$SYMBOL[which_probe_amp_del],as.character)),
        data_with_chr_bp_mapping$CHR_Mapping[which_probe_amp_del],
        data_with_chr_bp_mapping$BP_Mapping[which_probe_amp_del],
        data_with_chr_bp_mapping$BP_Mapping_v1[which_probe_amp_del],
        
        indicator_ica_4[which_probe_amp_del],
        
        rep(IC_no,length(which_probe_amp_del))))
      
      colnames(ampl_del_probesets) = c( "Value", "ENTREZID", 'GENETITLE', 'SYMBOL',"CHR_Mapping", "BP_Mapping", "BP_Mapping_for_plot", "Amp_or_Del","Component")
      
      ampl_del_probesets_summary = rbind(ampl_del_probesets_summary,ampl_del_probesets)
      
      run_vector = c(0, cumsum(indicator_run$lengths))
      
      amp_runs = which(abs(indicator_run$values)>0)
      
      for(num in 1:length(amp_runs)){
        
        rows_vec = which(original_mapping_file$ENTREZID%in%data_with_chr_bp_mapping$ENTREZID[(run_vector[amp_runs[num]]+1):run_vector[amp_runs[num]+1]])
        amplified_or_deleted_section_gene_level = as.data.frame(cbind(dataset[(run_vector[amp_runs[num]]+1):run_vector[amp_runs[num]+1],colnames(dataset)[IC_no]],data_with_chr_bp_mapping[(run_vector[amp_runs[num]]+1):run_vector[amp_runs[num]+1],])) 
        
        colnames(amplified_or_deleted_section_gene_level)[1] = "Value"
        
        if(indicator_run$values[amp_runs[num]]>0)
        {
          amplified_or_deleted_section_gene_level_sort_by_value = sqldf("select * from amplified_or_deleted_section_gene_level order by value DESC")
          
        }else{
          amplified_or_deleted_section_gene_level_sort_by_value = sqldf("select * from amplified_or_deleted_section_gene_level order by value")
          
        }
        
        write.table(amplified_or_deleted_section_gene_level,
                    file = paste(file_genelevel_bp,paste(Title,"_CHR", amplified_or_deleted_section_gene_level$CHR_Mapping[1],"_BPSTART",round(amplified_or_deleted_section_gene_level$BP_Mapping[1]), "_BPEND",
                                                         round(amplified_or_deleted_section_gene_level$BP_Mapping[dim(amplified_or_deleted_section_gene_level)[1]]),"_Component",colnames(dataset)[IC_no], ".txt", sep = ""),sep = "/"),
                    sep = "\t", row.names = FALSE)
        
        write.table(amplified_or_deleted_section_gene_level_sort_by_value,
                    file = paste(file_genelevel_value,paste(Title,"_CHR", amplified_or_deleted_section_gene_level_sort_by_value$CHR_Mapping[1],"_BPSTART",round(amplified_or_deleted_section_gene_level_sort_by_value$BP_Mapping[1]), "_BPEND",
                                                            round(amplified_or_deleted_section_gene_level_sort_by_value$BP_Mapping[dim(amplified_or_deleted_section_gene_level_sort_by_value)[1]]),"_Component",colnames(dataset)[IC_no], ".txt", sep = ""),sep = "/"),
                    sep = "\t", row.names = FALSE)
        
        
      }
      
      
      print(paste("Amplification or Deletion identified in some genomic postiions of Component",colnames(dataset)[IC_no]))
    }else{
      print(paste("No Genomic position has got amplification or deletion for Component",colnames(dataset)[IC_no]))
    }
    print(paste("Processing of Component", colnames(dataset)[IC_no], "FDR", FDR,  "CL", CL, "state_deciding_cutoff",state_deciding_cutoff, "is done"))
    
    print(paste("Time taken for this iteration is",round((proc.time()[3]-time1)/60,3),"mins"))
    
  }
  
  write.table(ampl_del_probesets_summary, file = paste("Genelevel",Title,
                                                       paste("All_Components_Amplified_or_Deleted_details",  "FDR", FDR, 
                                                             "CL",CL, "state_deciding_cutoff",
                                                             state_deciding_cutoff, sep = "_"),".txt", sep = "_"), sep = "\t", row.names = FALSE)
  
  
  
  
  
  