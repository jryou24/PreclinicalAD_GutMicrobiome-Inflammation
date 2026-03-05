read_tree <- function(path_to_nwk, path_to_species_list){
  ### This function uses Jian's code for reading in metaphlan trees and most of the comments for it are hers
  # Read in table with metaphlan4 sgb id's and corresponding phylogenies at species level
  
  tree <- ape::read.tree(path_to_nwk)
  tree_meta <- read.delim(path_to_species_list, header = F, sep='\t', check.names = F)

  tree_meta$V1=gsub('SGB|_group','',tree_meta$V1)
  tree_meta=separate(tree_meta,V2,into=c("V2","Other1"),sep=",")
  tree_meta=tree_meta[!duplicated(tree_meta[c('V2')]),]
  tree_meta2=subset(tree_meta,V1 %in% tree$tip.label)
  #subset only those species that are in the m4 tree
  tree=drop.tip(tree, setdiff(tree$tip.label, tree_meta2$V1))
  tree$tip.label=tree_meta2$V2
  return(tree)
}
