## Part2 Simulation of different sample pools 
## using parallel computing to test phylogenetic signal

## function for parallel computing for resampling tests
source("scripts/func_phylosig_parallel.R")

## libraries for parallel computing for resampling tests
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(stringr)

NumberOfCluster = detectCores() - 4
registerDoParallel(NumberOfCluster)


##########################################################################################
#################### 0. Species Pool and Resampling Database #############################
##########################################################################################

## define species pool

species_pool = TreeAllMatrix$tip.label



## Sample from 10 - 40 species, each sampling number have 100 replication sampling

resample_spp_min = 10
resample_spp_max = length(species_pool) - 5
resample_spp_rpt = 100

resample_spp_list = lapply(resample_spp_min:resample_spp_max, function(x) lapply(1:resample_spp_rpt, function(i) 
  sample(species_pool, x, replace = F)))

resample_spp_list_3 = sapply(1: (resample_spp_max - resample_spp_min + 1), function(i) as.data.frame(resample_spp_list[[i]]))

colnames(resample_spp_list_3) = sapply(1: (resample_spp_max - resample_spp_min + 1), function(i) colnames = paste0("resample_sp_", (resample_spp_min + i - 1)))

rownames(resample_spp_list_3) = paste0("repeat_", 1:resample_spp_rpt)

## Calculate Phylogenetic Signal of each sample pool

res = foreach (i = 1:(resample_spp_max - resample_spp_min + 1), .combine = rbind) %:% 
  ## i is the related withdraw number of species, column of resample_spp_list_3
  
  foreach(j=1:resample_spp_rpt, .combine = rbind) %dopar% {
    ## j is the repeat withdraw numbers
    resampling_spp = resample_spp_list_3[[100*(i-1)+j]]
    resample_spp = resample_spp_min + i - 1
    resample_rpt = j
    
    Tree = keep.tip(TreeAllMatrix, resampling_spp)
    
    dat_sp = dat %>% filter(dat$Species %in% Tree$tip.label)
    
    PhyloSig_para =  PhyloSig_para(Tree, dat_sp) %>% mutate(resample_spp = resample_spp) %>% mutate(resample_spp_rpt = resample_rpt) %>% mutate(resampling_spp = list(resampling_spp))
    
    as.data.frame(PhyloSig_para)
    
  } 


## explore the data in dplyr

phyloSig = res 

phyloSig$Blomberg.s.K = as.numeric(phyloSig$Blomberg.s.K)
phyloSig$P.value = as.numeric(phyloSig$P.value)

resam_ttl_var = (phyloSig %>% nrow())/4

## generate figures for Figure 3 and Figure 4

foreach(i = seq_along(unique(phyloSig$Variables)), .combine = rbind) %dopar% {
  
  trait = unique(phyloSig$Variables) [[i]]
  
  ########## character: all data
  
  ## i is the seed trait in analysis:
  phyloSig_i_full = phyloSig %>% dplyr::filter(Variables == trait) 
  
  phyloSig_i = phyloSig %>% dplyr::filter(Variables == trait) %>% 
    dplyr::filter(P.value <0.05) 
  #saveRDS(phyloSig_i, file = paste0('phyloSig_i_', trait, '.rds'))
  ## portion distribution along the increase number of sampling species
  ## portion distribution along the increase number of sampling species
  p_dis = phyloSig_i %>% group_by(resample_spp) %>% dplyr::summarise(cnt = n()) %>% mutate(perc = cnt/100) %>% arrange(resample_spp)
  
  pdf(file = file.path("results", paste0("resample_", trait, ".pdf")))
  print( ggplot(p_dis, aes(resample_spp, perc)) 
         + geom_point()
         + theme_bw()
         + ylim(0,1) +
           labs(x = "number of resample species",
                y = "purportion of phylogenetic signalin resamples (%)"))
  
  graphics.off()
  
  ## distribution of K-value
  ## p <= 0.05, col = "black", shape = 19
  ## p > 0.05, col = "grey70", shape = 1
  
  pdf(file = file.path("results", paste0("resample_k_distr_", trait, ".pdf")))
  
  print(ggplot(phyloSig_i_full, aes(resample_spp, Blomberg.s.K)) +
          geom_point(col = ifelse(phyloSig_i_full$P.value <= 0.05, 'black', 'grey70'), shape = ifelse(phyloSig_i_full$P.value <= 0.05, 19, 1), size = 1, position = "jitter") +
          theme_bw() + ylim(0, 1) +
          labs(x = "number of resample species",
               y = "the value of Blomberg's K"))
  
  
  graphics.off()
  
}  

## unbundle the cores

stopImplicitCluster()
