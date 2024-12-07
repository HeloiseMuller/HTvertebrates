##%######################################################%##
#                                                          #
####  This stage counts the number of HTT events that   ####
####              are required to explain               ####
####            the TE-TE hits (hit groups)             ####
#                                                          #
##%######################################################%##


library(data.table)
library(parallel)
library(stringi)
library(ape)
library(dplyr)

path = "~/Project/"
setwd(path) 

source("HTvFunctions.R")

# this script uses
# -the table of HTT hits clustered per clade at stage 10
httHits <- fread("occ200HitGroup_perclade.txt")

#Get coef to use to do copy * coef + hitGroup
coef_hitGroup <- 10^nchar(length(unique(httHits$hitGroup)))


# - the hit group statistics generated at stage 11 (imported later)
# - the "self-blastn" of TE copies against themselves, generated at stage 8

# The output is the table of HTT hits with a logical column indicating
# whether each transfer (hit group) can be "explained" by others (see below)

# PRINCIPLE:
# similar TEs may have been brought in two vertebrate clades A and B by several transfers ("hit groups")
# this may give the impression of a direct transfer between those clades
# hence this focal transfer may be "explained" by others
# if that is the case, the copies involved in the focal transfer are similar to other transfers
# and  the species involved in the focal transfer are the same as, or related to, 
# those of the explanatory transfers

# In the pseudocode section of the paper, a main function "requirement_1" evaluates these conditions
# but a naive implementation of this function would be vastly inefficient
# given how frequently this function will be called (thousands of times)
# In practice, most of the repetitive work underlying the pseudocode function
# is done upstream


################

# STEP ONE, we find TE copies that are similar between transfers --------------------------------------

# For each TE copy from each "focal transfer", we retrieve the highest sequence similarity
# that this copy has with any other copy from every other transfer.
# We do not yet control whether the two transfers involve related species
# Doing so allows changing the criterion for species relatedness independently
# of this first step, should we want to.

# Similarities are only evaluated within TE superfamilies
# since we never looked for homology between copies from different super families.
# So we split the hit table per super family, retaining only the columns we need 
# and replacing slashes with periods to avoid issues with file paths in future files
hitList <- split(
    x = httHits[, .(ID.1, ID.2,             # copy integer IDs
        hitGroup                      # the hitGroup identifiers
    )], 
    f = sub(".*[/]", "", httHits$rep_superF.1)
) 

# Note that all hit groups are taken, even those that we considered
# as unreliable in 11-hitGroupEvaluation.R. Doing so allow changing our
# selection of reliable transfers without having to re-run this step.

blastFiles <- list.files(
    path = "TEs/clustering/selfBlastn_perClade/out_catClades",
    pattern = ".out",
    full.names = T
)

# we split these file names by superfamily 
# for this we extract superfamily names from file names
superF <- sub(".out", "", basename(blastFiles))
blastFiles <- split(blastFiles, superF)

# this is the core function that finds similarities between 
# copies of different hit groups of a super family
findSimilarities <- function(superF) { # superF is the superfamily name (character)

    # we get the HTT hits of the super family being processed
    hits <- hitList[[superF]] 
    
    # this was only to get all the copies involved in these hits
    copies <- hits[, unique(c(ID.1, ID.2))] 

    # this function import blast files of copies against themselves
    selfBlast <- fread(
            input = blastFiles[[superF]],
            sep = "\t",
            header = T,
            select = c(1:3),
        )

    # we remove its not involving copies in the htt hits
    #vue que j'ai qu'un fichier par superF, ça ne devrait rien changer
    selfBlast <- selfBlast[query %in% copies & subject %in% copies]

    # and convert pID to integer for speed and RAM usage
    selfBlast[, pID := as.integer(pID * 1000)]
  
    # we will need to quickly retrieve all the copies that  
    # are homologous to a given one (including itself)
    # for this, we first duplicate the hits by creating reciprocal 
    # ones and add self hits with 100% pIDs
    selfBlast <- rbind(
        selfBlast, # the original hits
        selfBlast[, .(query = subject, subject = query, pID)], # the hits in "reverse"
        data.table(query = copies, subject = copies, pID = 100000L) # and the self hits
    ) 


    # we can now create the following list
    subjectsForQuery <- reList(split(selfBlast$subject, selfBlast$query))
    # subjectsForQuery[[x]] returns all copies that have some homology
    # with copy x (x is an integer)
    #N.B. Split make a list whose name is each query. Eacn list cotains subject(s) that have homology with this query 
    #N.B.B reList return a list of length max(query). subjectsForQuery[x] is null if x not in query
    
    # we will also need to quickly retrieve the pID score of all hits involving a query
    # the principle is the same as above
    idForQuery <- reList(split(selfBlast$pID, selfBlast$query))
    
    #N.B. subjectsForQuery & idForQuery have the same names() but the fist contaon subject, the second pID

    # as well as the hit group(s) in which each copy is found, following the same principle
    # but for this, we first have to obtain all ids of copies involved in each hit group
    copyPerHitGroup <- hits[, .(copy = unique(c(ID.1, ID.2))), by = hitGroup]
    hitGroupsForCopy <- reList(split(copyPerHitGroup$hitGroup, copyPerHitGroup$copy))
    #N.B. hitGroupsForCopy have the same name() as  subjectsForQuery & idForQuery, but contains hitGroup

    copyPerHitGroup <- split(copyPerHitGroup$copy, copyPerHitGroup$hitGroup)

    # we count the number of subjects per query, used later
    nSubjects <- lengths(subjectsForQuery) #la plupart font 0 cat l'indice ne correspond pas a une de mes copies

    # as well as the number of hit groups where each copy is found
    nHitGroups <- lengths(hitGroupsForCopy)

    rm(selfBlast) # to reclaim RAM, as there are lots of hits

    # for the copies of a given hit group ("focal" hit group), the function below
    # retrieves the best similarity this copy has with every other hit group
    findSimilaritiesForHitGroup <- function(copies, focalHitGroup) { 
    # focalHitGroup is just the integer id of the hit group
    
        # we first retrieve all copies that are homologous to 
        # those in the hit group (i.e., 'subject' of the blast)
        subject <- unlist(subjectsForQuery[copies]) #doing so we take only those without NULL

        # with the pIDs associated with this homologies
        pID <- unlist(idForQuery[copies]) #SAME

        # we place them in a data table (which requires replicating
        # copies that are homologous to several "subjects")
        pairs <- data.table(copy = rep(copies, nSubjects[copies]), subject, pID)

        # and we get the hit groups to which these subjects belongs
        explanatoryHitGroup <- unlist(hitGroupsForCopy[pairs$subject])

        # we place all this information in a table
        # for this, we need the number of hitGroups in which the subjects are found
        len <- nHitGroups[pairs$subject]

        # we generate the table (by replacing the previous one we no longer need)
        pairs <- pairs[, data.table(
            copy = rep(copy, len),
            pID = rep(pID, len),
            explanatoryHitGroup
        )]

        # we don't keep rows involving subjects that belong to the focal hit group
        pairs <- pairs[explanatoryHitGroup != focalHitGroup]

        # for each copy from focalHitGroup, we retain only the most similar subject per
        # explanatoryHitGroup
        setorder(pairs, -pID) 
        pairs <- pairs[!duplicated(data.table(copy, explanatoryHitGroup))]

        # we now add the focalHitGroup id to the table. 
        # Doing it earlier would have increased the memory footprint.
        pairs[, focalHitGroup := focalHitGroup]
        cat(".") # progress indicator
        pairs
    }

    # we apply the above function to all hit groups of the super family
    similarities <- Map(
        f = findSimilaritiesForHitGroup,
        copies = copyPerHitGroup, #a list with 1 element per hitGroup. Each element contains ID copies (go through elements one by one)
        focalHitGroup = as.integer(names(copyPerHitGroup)) #All the hitGroup in this superF (go through them one bu one)
    )

    # we write this table for safety, as it is returned by the function
    fwrite(rbindlist(similarities), stri_c("TEs/clustering/networks/", superF, ".similarities.txt"), sep='\t')

    rbindlist(similarities) #jean had hitGroupPairs but this variable does not exist
}

# we apply the function to all super families in parallel
res <- mclapply(
    X = names(hitList),
    FUN = findSimilarities,
    mc.cores = 10,
    mc.preschedule = F
)

fwrite(rbindlist(res),stri_c("TEs/clustering/networks/all.hitGroupSimilarities.txt"), sep='\t')


# We import the results
highestSimilarity <- fread(
  input = "TEs/clustering/networks/all.hitGroupSimilarities.txt",
  col.names = c("copy", "pID", "explanatoryTransfer", "focalTransfer")
  )

# from now on, we will refer to a hit group as a "transfer"
# Using this word helps to understand the concepts we use. 
# This required renaming some columns of the imported table
# in this table, the "copy" belongs to the "focalTransfer".
# Column ("pID") denotes the sequence similarity that the "copy" has with the
# most similar copy of "explanatoryTransfer" 
# (that other copy is not named in the table).
# the explanatoryTransfer may therefore "explain" the focalTransfer.
# But for this, explanatoryTransfer has to fulfil some conditions...




# STEP TWO  ---------------------------------------------------------------------
# we assess whether the explanatoryTransfer involves species 
# that are related to the host species of the "copy" from the focalTransfer

# We determine if the clade involved in the focalTransfer to which the copy belongs 
# is, nested in, or encompasses, either of the clades involved in the explanatoryTransfer. 
# Note that this doesn't require that species are shared between these transfers

# defining the clades requires the species tree
tree <- read.tree("datasetTree.nwk")

# We first retrieve the host species of each copy of the focal transfer,
# we encode species as tip numbers of the timetree
spForCopy <- integer(httHits[,max(ID.1, ID.2)])
spForCopy[httHits[,c(ID.1, ID.2)]] <- httHits[,chmatch(c(species.1, species.2), tree$tip.label)]

# we add the host species in a new column
highestSimilarity[, speciesA := spForCopy[copy]]


# We define two clades involved in each transfer by the mrcas of the species involved
# which are node numbers in the species tree
clades <- httHits[, .(
    cladeA = MRCA(tree, unique(species.1)), # the mrca of the left-clade species (sp1)
    cladeB = MRCA(tree, unique(species.2))  # same for the right-clade species (sp2)
), 
by = .(transfer = hitGroup)
]

# we make integer vectors to quickly retrieve the clades (A or B) of each transfer
cladeA <- cladeB <- integer(max(clades$transfer))
cladeA[clades$transfer] <- clades$cladeA
cladeB[clades$transfer] <- clades$cladeB


# We now attribute each copy involved of the focalTransfer to a clade.
# We start with "query" copies, which belong to clade A
queries <- httHits[,.(copy = unique(ID.1)), by =  hitGroup]
queries[,clade := cladeA[hitGroup]]

# we do the same for "subject" copies, those belonging to clade B
subjects <- httHits[,.(copy = unique(ID.2)), by = hitGroup]
subjects[,clade := cladeB[hitGroup]]

cladeForCopies = rbind(queries, subjects)

# we create a unique integer identifier for each copy in each transfer
cladeForCopies[, copyTransfer := copy * coef_hitGroup + hitGroup] 
 
# we do the same for the highestSimilarity table
highestSimilarity[, copyfocalTransfer := copy * coef_hitGroup + focalTransfer] 

# We now add a column denoting the clade to which the copy belongs in the focalTransfer
# We call this column "clade_A" to match the name we use in the pseudocode
highestSimilarity[,clade_A := cladeForCopies[match(
  copyfocalTransfer,
  copyTransfer
), clade]]

# we add similar columns indicating the clades involved in the explanatoryTransfer
# we use "cladeC" to match the name used in supplementary figure 4
highestSimilarity[, c("clade_C1", "clade_C2") := .(
  cladeA[explanatoryTransfer],
  cladeB[explanatoryTransfer]
)]


#----------
# we create a matrix that is TRUE if a clade (column/row index) is nested in,
# or includes, another clade (row/column index). The diagonal is also TRUE
# see nestedClades() function in HTvFunctions.R
nested <- nestedClades(tree)

# we add a logical column that tells if cladeA is
# nested, or includes, either clade of explanatoryTransfer
highestSimilarity[, nestedOrIncludes := nested[cbind(clade_A, clade_C1)] | 
                        nested[cbind(clade_A, clade_C2)]]




# STEP THREE ------------------------------------------------------------------------
# We prepare our criterion that imposes that the best identity that a copy 
# of the focalTransfer has with any copy of the explanatoryTransfer is higher than
# the lowest identity the copy has within the focalTransfer

# for this, we obtain the min pID that each copy has within each transfer
minCopyIDs <- rbind(
    httHits[, .(lowestSimilarity = min(pID)), by = .(hitGroup, copy = ID.1)], # for "query" copies
    httHits[, .(lowestSimilarity = min(pID)), by = .(hitGroup, copy = ID.2)]  # and for "subject" copies
) 

# as we did previously, we use an integer to identify the copy-transfer pair
minCopyIDs[, copyTransfer := copy * coef_hitGroup + hitGroup] 

# we can now add the column denoting the lowest_similarity
# we convert this pID to integer, to compare it 
# to the pID column
highestSimilarity[, lowestSimilarity := minCopyIDs[
    match(copyfocalTransfer, copyTransfer),
    as.integer(lowestSimilarity * 1000)
]] 



# we now select pairs of homologous copies fulfilling our conditions -----------------------------------

# We will also discard the transfers that were considered unreliable in the previous stage.
# Doing it only now allows changing the retained transfers 
# without having to re-run the whole script

# We thus import the statistics we generated in the previous script
hitGroupStats <- fread("TEs/clustering/hitGroupStats_perClade.txt")

retainedTransfers <- hitGroupStats[retained == T, hitGroup]

selectedCopies <- highestSimilarity[, 
                                   
    # the best identity it has with a copy of the explanatory transfer
    # must be higher than the the lowest similarity it has within the focal transfer
    pID > lowestSimilarity & 
  
    # its host clade must be nested in, or encompass, either clade of
    # the explanatory transfer
    nestedOrIncludes == TRUE &

    # the transfers must be among those considered as reliable
    explanatoryTransfer %in% retainedTransfers & focalTransfer %in% retainedTransfers
    ]


# STEP FIVE ----------------------------------------------------------------------------
# we evaluate whether explanatory transfers could have brought TE copies
# in all the species composing the focalTransfer


# We make several objects (lists) to optimize the procedure.
# We first list the species that carry copies that are
# similar to those of other transfers harboring related species

# at this stage, we no longer care about TE copies, 
# so we simplify the table with unique()
explainedSpecies <- unique(highestSimilarity[selectedCopies, .(
  focalTransfer, 
  explanatoryTransfer, 
  speciesA
  )])

# we split speciesA of this table by explanatoryTransfer then by focalTransfer
# we use a modification of the split() function that adds recursiveness
explainedSpecies <- Split(
    x = explainedSpecies$speciesA,
    f = explainedSpecies[, list(explanatoryTransfer, focalTransfer)],
    drop = T,
    recursive = T
)

# explainedSpecies[[x]][[y]] returns the species that carry the TEs 
# involved in focal transfer "y" that may have been brought 
# by explanatory transfer "x"
# names(explainedSpecies[[x]]) gives the ids of transfers that are  
# partly explained by transfer x

explainedSpecies <- reList(explainedSpecies, max(retainedTransfers))
# now x can be an integer, which speeds-up access to first-level elements


# As we will check whether all species of the focal transfers have TEs that
# may have been brought by others transfers, we need to list all species per transfer
# we encode species as tree tip numbers
spForTransfer <- httHits[, unique(chmatch(c(species.1, species.2), tree$tip.label)), 
    by = hitGroup]

# We again make a list for quick access. spForTransfer[[x]] will return all species
# involved in transfer "x", x being an integer
spForTransfer <- reList(split(
    x = spForTransfer$V1,
    f = spForTransfer$hitGroup
))


# These lists were generated to optimize the speed of the function below
# This is the function that tells if a focal transfer can be explained by others
# i.e., if "requirement 1" as defined in the methods and pseudocode is passed

requirement1_passed <- function(transfer) {

    # we retrieve all species whose TE copies may have been brought by explanatory transfers
    explSpecies <- unlist(
      explainedSpecies[explanatoryTransfers], 
      recursive = F)
    
    # we of course only consider the species from the transfer we investigate
    explSpecies <- explSpecies[names(explSpecies) == transfer]
    # this selection is not very efficient, but a faster version of 
    # the function was more complex to understand and distant from the 
    # pseudocode
    
    # we return whether these species constitute all species of the focal transfer
    # and if there are at least 2 contributing transfers
    all(spForTransfer[[as.integer(transfer)]] %in% 
          unlist(explSpecies, use.names = F)) & 
      length(explSpecies) >= 2L
}



# we evaluate requirement 1 on transfers iteratively, sorted by "reliability" score ---------------------------------

# to determine if a transfer can be explained by others, 
# we will inspect the less "reliable" transfers first
# these are considered less likely to represent a 
# "direct" transfer event between clades

# we put the best htt hits on top, for each transfer
setorder(httHits, hitGroup, -pID)

# the reliability score of a transfer is based on 
# the pID of best hits of copies involved (see Methods)
# the sum of the best pIDs for over "query" copies (cladeA)
pIDQ <- httHits[!duplicated(data.table(hitGroup, ID.1)), .(sumID = sum(pID)), by = hitGroup]

# and for the subject copies (cladeB)
pIDS <- httHits[!duplicated(data.table(hitGroup, ID.2)), .(sumID = sum(pID)), by = hitGroup]

# we add the "score" as a new column (the lowest of the two sum of pIDs)
hitGroupStats[, score := pmin(pIDQ$sumID, pIDS$sumID)]

# we extract retained transfers, ordered by reliability score, as an integer vector
orderedTransfers <- hitGroupStats[retained == T, hitGroup[order(score)]]


# we are now ready to iterate over transfers -----------------------
# we initialise two vectors:
# this one will contain the identifiers of transfers that are explained by others
explained <- NULL

# and this one contains the transfers than may explain others
# initially, all transfers are allowed to explain others
# but if one is explained, hence "indirect", it does not
# correspond to a "movement" of TEs, and therefore it cannot 
# explain other transfers
explanatoryTransfers <- retainedTransfers

#prend qq heures
for (focalTransfer in orderedTransfers) {
    if (requirement1_passed(focalTransfer)) {
        # if the transfer can be explained by others,
        # we may remove it from explanatory transfers
        explanatoryTransfers <- setdiff(explanatoryTransfers, focalTransfer)
      
        # we investigate if transfers that were explained by the focal one
        # can still be explained without it
      
        # we define the transfers to investigate
        toInvestigate <- intersect(
          explained, 
          names(explainedSpecies[[focalTransfer]])
          )
        
        # we check if all these transfers can 
        # still be explained without the focal transfer
        stillExplained <- sapply(toInvestigate, requirement1_passed)
        
        if (any(!stillExplained)) {
          
          # if any of these transfers can no longer be explained 
          # after removal of the focal transfer we restore the 
          # focal transfer into the list of explanatory transfers
          explanatoryTransfers <- c(explanatoryTransfers, focalTransfer)
          
          cat(".") # progress indicator
          
        } else {
          # else we consider the focal transfer as explained
          explained = c(explained, focalTransfer)
        }
    }
}

# WE ARE NOW FINISHED WITH THE COUNT OF HTT EVENTS -------------------------------------



# we only retain the columns we need in the htt hit table
httHits <- httHits[, .(
  copie1, copie2, ID.1, ID.2,
  TEconsensus.1, TEconsensus.2, #I added these 2 at the version greedy
  species.1, species.2, rep_superF.1,
  mrca, divTime, pID, length, qStart, qEnd, 
  sStart, sEnd, dN, dS, length.aa, com, hitGroup
  )] 

# we add two logical columns:
# - "retained" is TRUE for transfers we retained (based on the contamination filter , sufficient divergent time and low dS)
# - "independent" is TRUE is a transfer is not explained by other 
# i.e. a htt event that may be seen as "independent" or "direct"
httHits[, c("retained", "independent") := .(hitGroup %in% hitGroupStats[retained==T,]$hitGroup, ! hitGroup %in% explained)] 

# we save this whole table, including transfers that are not retained
fwrite(httHits, "HTThitsAssessed_perClade.txt", sep='\t') 


# we remove transfers and columns we don't retain
retainedHits <- httHits[retained == T, -"retained"]

# we replace column names with more user-friendly ones
setnames(
    x = retainedHits,
    old = c("copie1", "copie2", "com", "rep_superF.1"),
    new = c("copy1", "copy2", "community", "superfamily")
)

# we add more common super family names to the hits
retainedHits[, c("subClass", "newName") := .(splitToColumns(superfamily, "/", column=1), splitToColumns(superfamily, "/", column=2))] 


fwrite(retainedHits, "HTThitsAssessed_perClade_retained.txt", sep='\t')