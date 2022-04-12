
##############################################################################################################
###################### Function of Pruning and Replace Species on a Tree ######################################
##############################################################################################################

## The functioon help to Prun the tree 
## also place the species according to its general, if the exact species is not on the tree
## limitation: only one species can use replacement method

PruneTree <- function(desir_species, master_tree){
  ## desir_species is a pool of species, can be unique or not, in "genera_species" format
  ## master_tree is a species tree, normally sythesize from large number of sequence data. 
  
  ## prepare the master tree
  tips <-master_tree$tip.label
  Genera<-unique(sapply(strsplit(tips,"_"),function(x) x[1]))
  
  
  ## double check the data
  DesiredSpecies <- unique(desir_species)
  ## Desired Genera
  DesiredGenera <- sapply(strsplit(DesiredSpecies, "_"),function(x) x[1])
  
  ## TreeSpecies Data Frame
  SpeciesGeneraSpecies <- data_frame(TreeSpecies = master_tree$tip.label,
                                     TreeGenera = sapply(strsplit(tips,"_"),function(x) x[1]),
                                     DesiringGenera = TreeGenera %in% intersect(Genera,DesiredGenera) ,
                                     DesiringSpecies = TreeSpecies %in% DesiredSpecies) 
  
  SpeciesListSpecies <- filter(SpeciesGeneraSpecies, DesiringSpecies == "TRUE")
  
  GeneraSpecies <- filter(SpeciesGeneraSpecies, SpeciesGeneraSpecies$TreeGenera %in% 
                            setdiff(DesiredGenera, SpeciesListSpecies$TreeGenera))
  SpeciesListGenera <- group_by(GeneraSpecies, TreeGenera) %>% group_modify(~ head(.x, 1L))
  
  LISTALLSPECIES <- rbind.data.frame(SpeciesListSpecies, SpeciesListGenera)
  
  treeTestedSpecies <- keep.tip(master_tree, LISTALLSPECIES$TreeSpecies)
  
  ## Replacing tip label
  
  aaa <- data_frame(DesiredSpecies = DesiredSpecies, 
                    TreeGenera = sapply(strsplit(DesiredSpecies,"_"),function(x) x[1]))
  bbb <- merge(aaa, SpeciesListGenera)
  
  treeTestedSpecies$tip.label <- mapvalues(treeTestedSpecies$tip.label, c(bbb$TreeSpecies), c(bbb$DesiredSpecies))
  tree_x <- treeTestedSpecies
#  plotTree(tree_x, ftype="i")
  return(tree_x)
}
