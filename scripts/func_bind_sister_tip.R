## function binding tips to their sister branch
## the tip relationship were assumed based on taxonomic relationship
## the branch length were assumed half of the branch length


bind_sister = function(new_tip, sister_tip, tree){
  # sister_tip is the tip were assume phylogenetically closely related based on there taxonomic relationship
  # tree is the target binding phylogeny
  newtree = bind.tip(tree,new_tip,
           edge.length = 0.5*tree$edge.length[which(tree$edge[,2]==
                                                                 which(tree$tip.label==sister_tip))],
           where=which(tree$tip.label==sister_tip),
           position=0.5*tree$edge.length[which(tree$edge[,2]==
                                                            which(tree$tip.label==sister_tip))])
  return(newtree)
}
