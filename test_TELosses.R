path = "~/Project/"
setwd(path) 

source("HTvFunctions.R")
 
# to infer TE losses, we rely on transfer events that involve more than one species
# in a given clade involved. If these species do not compose a monophyletic group,
# this can indicate a loss of TEs after the transfer.

htt = fread("HTThitsAssessed_perClade_retained.txt")
tree = read.tree("datasetTree.nwk")

# For each "independent" transfer, we obtain the species and clade involved
# As transfers, we use "communities" of TE-TE hits (see paper), because hitgroups may frequently
# contain TEs from several transfers, due to the aggressive clustering we applied,
# especially the step where we cluster communities if they don't show homology in TE protein
# sequences (and might *in principle* represent non overlapping parts of the same TE).
sp_clade = htt[independent == T, .(sp = c(species.1, species.2), transfer = c(community, community), clade = rep(1:2, each = .N))]
sp_clade = unique(sp_clade)
sp_clade[,tr_c := paste(transfer, clade)]
nsp = sp_clade[,.N, by = tr_c] # number of species per clade in each transfer

# We select cases involving more than 1 species per clade
sp_clade = sp_clade[tr_c %in% nsp[N > 1, tr_c],]

# we obtain the tip number of the species (in the tree) in a list 
# each element in this list corresponds to a transfer+clade
tipIDs = split(match(sp_clade$sp, tree$tip.label), sp_clade$tr_c)

# This allows us to obtain the MRCAs of each clade in each transfer
mrcas = sapply(tipIDs, function(sp) getMRCA(tree, sp))

# We then obtain all the descendants of each MRCA (tips of the tree)
tips = lapply(mrcas, function(node) tipsForNode(tree, node))

# We can them get the species that are "missing", for each transfer, by applying a set difference
diffSpecies = Map(setdiff, tips, tipIDs)

# and count the number of monophyletic group they represent (for each transfer)
# This is a conservative estimate of the number of TE losses.

cladesForTips = function(tree, tips) {
  # Returns the oldest clades (node numbers) leading to tips of a tree, 
  # ensuring that each clade contains only these tips (and no other)
  if(length(tips) == 0) {
    return(integer(0))
  }
  if(is.character(tips)) {
    tips = match(tips, tree$tip.label)
  }
  edges = data.table(tree$edge)
  test = T
  nodes = NULL
  currentNodes = tips
  while(test) {
    anc = edges[V2 %in% currentNodes]
    descs = tipsForNodes(tree, anc$V1)
    good = descs[,all(tip %in% tips), by = node]
    f = good$V1
    test = any(f)
    if(test) {
      nodes = c(nodes, anc[V1 %in% good[!f, node], V2])
      currentNodes = good[f, node]
    }
  }
  nodes
}


losses = lapply(diffSpecies, function(tips) cladesForTips(tree, tips))
table(lengths(losses) > 0) # only 7 transfers involving possible TE losses. 


