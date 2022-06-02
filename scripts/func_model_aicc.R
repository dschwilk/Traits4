
##############################################################################################################
###################### Function of Model Selection, based on AICc value ######################################
##############################################################################################################

#### Requirement inputs include phylogenetic tree and its comparable trait variable



Model_AICc = function(TreeAllMatrix, dat){
  
  ## TreeAllMatrix are the species tree after PruneTree function
  ## dat is a list of trait value of different species
  
  dist = cophenetic.phylo(TreeAllMatrix)
  phyloposi = MASS::isoMDS(dist, trace = F) %>% as.data.frame()
  
  phyloposi_species = phyloposi %>% mutate(Species = row.names(phyloposi)) %>% 
    mutate(phy1 = round(scale(points.1),digits = 1)) %>% 
    mutate(phy2 = round(scale(points.2), digits = 1))
  
  AllMatrix = merge(dat, phyloposi_species)
  
  AllMatrix_G = merge(AllMatrix, Germination)
  
  AllMatrix_G_scaled = AllMatrix_G %>%
    mutate(scaled_Area = scale(mean_Area)) %>%
    mutate(scaled_Height = scale(mean_Height)) %>%
    mutate(scaled_Mass = scale(mean_Mass)) %>%
    mutate(scaled_Germination = scale(mean_Germination)) %>%
    dplyr::select(scaled_Area, scaled_Height, scaled_Mass, scaled_Germination, phy1, phy2)
  
  
  vars = names(AllMatrix_G_scaled[-4]) 
  models = list()
  for (i in seq_along(vars)){
    vc = combn(vars, i)
    for (j in 1:ncol(vc)){
      model = as.formula(paste0('scaled_Germination ~', paste0(vc[, j], collapse = '+')))
      models = c(models, model)}}
  
  
  glmmodels = lapply(models, function(x) glm(x, data = AllMatrix_G_scaled))
  
  AICc = sapply(models, function(x) AICcmodavg::AICc(glm(x, data = AllMatrix_G_scaled), 
                                         return.K = FALSE))
  
  names(AICc) = sapply(glmmodels, function(x) x$formula)
  
  tab_model = data_frame(names(AICc),AICc)
  colnames(tab_model) = c("Model", "AICc")
  
  #tab_models = tab_model %>% mutate(resample_spp = resample_spp) %>%
   # mutate(resample_rpt = resample_rpt)
  
  tab_models = tab_model[order(AICc),]
  
  tab_models
  
}