## set working directory
#setwd("~/Desktop/Research Projects/Traits/GitRepo/Traits4") 

## source the scripts: data and functions
source("Scripts/lib_dat_prep.R")
source("Scripts/func_prun_replac.R")
source("Scripts/func_phylosig.R")
source("Scripts/func_phylosig_parallel.R")
source("Scripts/func_model_aicc.R")


library(dplyr)
library(ggplot2)
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

species_pool = c("Eryngium_leavenworthii","Polytaenia_nuttalli","Asclepias_asperula",
                 "Centaurea_americana","Coreopsis_lanceolata","Coreopsis_tinctoria",
                 "Echinacea_angustifolia","Gutierrezia_sarothrae","Helianthus_annuus",
                 "Liatris_mucronata","Ratibida_columnifera","Tradescantia_occidentalis",
                 "Astragalus_crassicarpus","Desmanthus_illinoensis","Corydalis_curvisiliqua",
                 "Phacelia_congesta","Herbertia_lahue","Monarda_citriodora","Salvia_azurea",
                 "Salvia_coccinea","Salvia_farinacea","Salvia_lyrata","Linum_rigidum",
                 "Callirhoe_involucrata","Pavonia_lasiopetala","Oenothera_rhombipetala",
                 "Argemone_albiflora","Phytolacca_americana","Rivina_humilis","Andropogon_gerardii",
                 "Aristida_purpurea","Bouteloua_curtipendula","Bouteloua_gracilis", 
                 "Chasmanthium_atifolium","Chloris_cucullata",
                 "Digitaria_ciliaris","Eragrostis_trichodes","Schizachyrium_scoparium",
                 "Sorghastrum_nutans","Sporobolus_airoides","Sporobolus_cryptandrus","Ipomopsis_rubra")

## Sample from 10 - 40 species, each sampling number have 100 replication sampling

resample_spp_min = 10
resample_spp_max = length(species_pool) - 2
resample_spp_rpt = 100

resample_spp_list = lapply(resample_spp_min:resample_spp_max, function(x) lapply(1:resample_spp_rpt, function(i) 
  sample(species_pool, x, replace = F)))

resample_spp_list_3 = sapply(1: (resample_spp_max - resample_spp_min + 1), function(i) as.data.frame(resample_spp_list[[i]]))

colnames(resample_spp_list_3) = sapply(1: (resample_spp_max - resample_spp_min + 1), function(i) colnames = paste0("resample_sp_", (resample_spp_min + i - 1)))

rownames(resample_spp_list_3) = paste0("repeat_", 1:resample_spp_rpt)

## Calculate Phylogenetic Signal of each sample pool

res2 = foreach (i = 1:(resample_spp_max - resample_spp_min + 1), .combine = rbind) %:% 
  ## i is the related withdraw number of species, column of resample_spp_list_3
  
  foreach(j=1:resample_spp_rpt, .combine = rbind) %dopar% {
    ## j is the repeat withdraw numbers
    resampling_spp = resample_spp_list_3[[100*(i-1)+j]]
    resample_spp = resample_spp_min + i - 1
    resample_rpt = j
    
    TreeAllMatrix = PruneTree(resampling_spp, treeVascularPlants)
    
    dat_sp = dat %>% filter(dat$Species %in% TreeAllMatrix$tip.label)
    
    PhyloSig_para =  PhyloSig_para(TreeAllMatrix, dat_sp) %>% mutate(resample_spp = resample_spp) %>% mutate(resample_spp_rpt = resample_rpt) %>% mutate(resampling_spp = list(resampling_spp))
    
    as.data.frame(PhyloSig_para)
    
  } 


## explore the data in dplyr

phyloSig1 = res2 

phyloSig1$Blomberg.s.K = as.numeric(phyloSig1$Blomberg.s.K)
phyloSig1$P.value = as.numeric(phyloSig1$P.value)

resam_ttl_var = (phyloSig1 %>% nrow())/4

dir.create("./Output")
setwd("./Output")

foreach(i = seq_along(unique(phyloSig1$Variables)), .combine = rbind) %dopar% {
  
  trait = unique(phyloSig1$Variables) [[i]]
  
  ########## character: all data
  
  ## i is the seed trait in analysis:
  phyloSig_i_full = phyloSig1 %>% dplyr::filter(Variables == trait) 
  
  phyloSig_i = phyloSig1 %>% dplyr::filter(Variables == trait) %>% 
    dplyr::filter(P.value <0.05) 
  #saveRDS(phyloSig_i, file = paste0('phyloSig_i_', trait, '.rds'))
  ## portion distribution along the increase number of sampling species
  ## portion distribution along the increase number of sampling species
  p_dis = phyloSig_i %>% group_by(resample_spp) %>% dplyr::summarise(cnt = n()) %>% mutate(perc = cnt/100) %>% arrange(resample_spp)
  
  pdf(file = paste0("resample_", trait, ".pdf"))
  print( ggplot(p_dis, aes(resample_spp, perc)) 
         + geom_point()
         + theme_bw()
         + ylim(0,1) +
           labs(x = "number of resample species",
                y = "purportion of phylogenetic signalin resamples (%)"))
  
  graphics.off()
  
  ## distribution of K-value
  ## p <= 0.05, black
  ## p > 0.05, grey86
  
  pdf(file = paste0("resample_k_distr_", trait, ".pdf"))
  print(ggplot(phyloSig_i_full, aes(resample_spp, Blomberg.s.K)) +
          geom_point(col = ifelse(phyloSig_i_full$P.value <= 0.05, 'black', 'grey97'), size = 0.05) +
          theme_bw() + ylim(0, 1) +
          labs(x = "number of resample species",
               y = "the value of Blomberg's K"))
  graphics.off()
  
}  

## unbundle the cores

stopImplicitCluster()
