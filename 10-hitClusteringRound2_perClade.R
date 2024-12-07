## %######################################################%##
#                                                          #
####   This stage applies the 2nd round of clustering   ####
#                                                          #
## %######################################################%##

library(data.table)
library(dplyr)
library(stringr)
library(stringi)
library(ape)
library(Biostrings)

path = "~/Project/"
setwd(path) 

source("HTvFunctions.R")

# We cluster hits from different communities into "hit groups" (see the paper methods for details)

# this script uses
# - the tabular file of selected TE-TE hit ("HTT hits") from stage 07-TEKsAndHTTfilter.R
# with community IDs added at stage 09
httHits <- fread("occ200ds05_comm_perClade.txt") 
# - the self-blastn output generated at stage 8
# - the timetree of the species (both imported later)

# output is
# - a tabular file of HTT hits with a new column attributing a hit-group number to each, "ocCD00_comm.txt"
# - a tabular file of statistics for each pair of hit community (much too long to describe in a comment)


# STEP ONE, we apply criterion 1 to hits from different communities---------------------------------------------------------------
# this is done for communities of hits between TEs of the same super family and MRCA
# Indeed, two hits involving pairs of species that do not have the
# same MRCA cannot results from the same HTT (see paper for details)
# The procedure is similar to what we did in the previous stage

# here a "group" of hits to cluster is defined by the TE super family and MRCA
# we therefore add a column indicating the "group" of each hit

httHits[, group := stri_c(gsub("/", ".", rep_superF.1, fixed = T), "_", mrca)]

# like in the previous stage, we split the hits by group and select the column we need for the clustering
hitList <- split(
    x = httHits[, .(
        #NOTE Because of the following, hit numbering is different than in table httHits.
        #It is not a problem because hitList is only for iterativeSecondClustering, 
            #that output all.txt, which does not contain the column hit
            #and that does not use httHits
        hit = 1:.N, 
        query = ID.1,
        subject = ID.2,
        pID = as.integer(pID * 1000L), # pID is again converted to integer for speed
        group,
        com
    )],
    f = httHits$group
)

# we determine the number of hits per group. 
# No clustering need be done if there is just one community of hits per group
nCom <- sapply(
    X = hitList,
    FUN = function(x) length(unique(x$com))
)

hitList <- hitList[nCom > 1L] 

# and we sort the hit list by decreasing number of hits (for better CPU usage)
hitList <- hitList[order(sapply(
    X = hitList,
    FUN = nrow
),
decreasing = T
)]

# results of this round of clustering will go there:
dir.create("TEs/clustering/round2") 

# we save the list of hits that will be clustered by a dedicated script
saveRDS(
    object = hitList,
    file = "TEs/clustering/round2/hitsFor2ndClustering.RDS"
)

system("Rscript iterativeSecondClustering.R 20") # applies criterion 1 with 20 CPUs. (Job=3544208)

# we collect results. 
# these are statistics for every relevant pair of hit communities (See iterativeSecondClustering.R)
# among groups with >1 com, stats shows all pair of com
stats <- fread("TEs/clustering/round2/all.txt")

# We now evaluate the connections between hit communities that may represent the same HTTs
# the number of hit pairs that pass criterion 1 between two communities
stats[, links := tot - allDiff] 

# criterion 1 needs to be passed by more than 5% of hit pairs
# to consider the communities as representing the same HTT.
stats[, crit1 := links / tot > 0.05]

# But we are not done with criterion 1 yet (see STEP THREE)



# STEP TWO, applying "criterion 2" -------------------------------------------------------------------------------------
## we evaluate whether 2 communities can represent the same HTT, based on
# "criterion 2", which relies on the inferred age of the transfer. To reflect the
# same HTT, the transfers (represented by the 2 communities of hits) should not be
# more recent that both clades. We infer the maximal times of transfers from
# the mean Ks between copies at hits of the same community. We compare it do the
# Ks of buscos of corresponding clades. Se we must determine the clades involved. For 2
# clusters of hits between young clades A-B (community 1) and C-D (community 2)
# both coalescing to the same MRCA (mrca of A and B is the same that of C and D),
# 2 clades are involved. One is composed of A-C and the order of B-D, if
# A and C belong to one of the 2 subclades diverging from the MRCA. Remember that
# we took care to place species of one subclade in column sp1 and species of the
# other subclade in column sp2, for each mrca. (done at stage 08-prepareClustering.R)

# for each community, we compute the mean dS of TE-TE hits
dS_stats <- httHits[, .(meandS = mean(dS)), by = .(com, mrca)]

# we get the core gene Ks threshold we used to filter hits between species of a given MRCA
# we thus import BUSCO gene Ks (200 AA alignments and one score per BUSCO gene per 
# pair of clades) generated in stage 05-coreGeneKs.R
dS <- fread("Busco/dS/dNdS/dS100AA_filterBusco_median.txt") 

# for each pair of subclades diverging from an MRCA ("clade"), we get the dS threshold mentioned above
dS_thresholds <- dS[, .(q05 = quantile(dSmedian, 0.005)), by = mrca]

# we add the dS threshold as a new column
dS_stats = left_join(dS_stats, dS_thresholds, by="mrca") %>% dplyr::rename(., "threshold"=q05)

# we identify the 2 clades (A-C and B-D in the explanation above) that are 
# involved in every pair of communities ----------------
# for this we need the timetree
tree <- read.tree("datasetTree.nwk")

# the matrix of MRCA for all species of the tree
mrca <- mrca(tree) 

# We identify the clades involved in the HTT of the communities.
# we first retrieve two species for each community of a pair, those involved in the first hit
species <- cbind(
    httHits[match(stats$com1, com), cbind(species.1, species.2)],
    httHits[match(stats$com2, com), cbind(species.1, species.2)]
)

# since all sp1 and sp2 from a hit community are from different "young" clades of <40 My,
# their MRCA define the clades we want (mrca of clades A-B and of clades C-D)

# we add these MRCA to the table
stats[, c("AB", "CD") := .(
    mrca[species[, c(1, 3)]], # AB defines the MRCA of clades A and B
    mrca[species[, c(2, 4)]]
)]

# we add a last row for a dummy clade with Ks 0
# this will be useful later
dS_thresholds <- rbind(
    dS_thresholds,
    data.table(mrca = 0L, q05 = 0)
)

# the "nomatch" argument below below is set to the last row, hence the q05 value retrieved would be 0
# Ks would thus be 0 for an MRCA for which we did not compute Ks (between its children clades)
# this is the case for lineage that diverged within the last 40 My
stats[, c("dSAB", "dS") := .(
    dS_thresholds[match(AB, mrca, nomatch = .N), q05],
    dS_thresholds[match(CD, mrca, nomatch = .N), q05]
)]

# and we retrieve the mean Ks of the hit communities involved
stats[, c("dS1", "dS2") :=
    data.table(
        dS_stats[match(com1, com), meandS],
        dS_stats[match(com2, com), meandS]
    )]

# we are now ready to apply criterion 2 by adding a logical column to this effect
stats[, crit2 := (pmin(dS1, dS2) >= pmin(dSAB, dSCD) &
    pmax(dS1, dS2) >= pmax(dSAB, dSCD))]



# STEP THREE ------------------------------------------------------------------------------------------
# to check whether different clusters of hits could represent the
# retention of different parts of a TE transferred once, we compare
# their protein regions.
# the approach is that used in Peccoud et al. 2017 PNAS

# For this we reuse the blastx results of copies against repeat proteins (done in stage 06-filterTEhits.R)
blastx <- fread("TEs/blastx/all.copies.successiveBlastx.out", 
      col.names = c("query", "subject", "pID", "length", "mismatches", "gapOpen",
        "qStart", "qEnd", "sStart", "sEnd", "evalue", "score"))


# we retrieve copy integer IDs, to speed up some functions
copyIDs <-fread("TEs/clustering/dcMegablast/IDcopies_full.txt")

# we use this id to refer to queries of the blastx results
blastx[, query := copyIDs[chmatch(query, copy), ID]]

# to reduce the workload, we discard hits that do not  
# involve copies that we retained in the pipeline
blastx <- blastx[!is.na(query)]

#I think we could even reduce workload more doing this:
# blastx[query %in% unique(c(httHits$ID.1,httHits$ID.2)),]
#Just make sure we don't need the other copies

# we blastp proteins that are similar to copies against themselves. -----------
# This will be used to determine whether two copies have some homologies
# that is, is they have homology to protein regions that are themselves homologous. 
# This should be quite sensitive.

# we get the protein that are similar to retained copies
prot <- blastx[, unique(subject)] 

# we these proteins are from the database used by repeat modeler
# we therefore import it
protSeqs <- readAAStringSet("DB_clustered/findDubious/RepeatPeps.lib") 

# we have to discard protein descriptions in protein names
names(protSeqs) <- splitToColumns(names(protSeqs), " ", "",  1) 

# results of the blastp search will go there
dir.create("TEs/blastp") 

# we create a fasta of all proteins that are similar to TEs
writeXStringSet(
    x = protSeqs[prot], #save only prot that have hits with our TEs
    filepath = "TEs/blastp/involvedProtOCC200dS05.fas"
)
#NOTE: we continue script with protSeqs, but only a subset was saved

# we make a blastp databse
system("makeblastdb -in TEs/blastp/involvedProtOCC200dS05.fas -dbtype prot")

# and launch the blastp search
system(paste(
    "blastp -query TEs/blastp/involvedProtOCC200dS05.fas -db TEs/blastp/involvedProtOCC200dS05.fas -max_target_seqs",
    
    # we set max_target_seqs as the query size to ensure that all possible hits are reported
    length(prot), 
    "-outfmt 6 -num_threads 5 -evalue 1e-4 -out TEs/blastp/involvedProtOCC200dS05.self.out"
))

# imports results of the self blastp launched above
blastp <- fread(
    input = "TEs/blastp/involvedProtOCC200dS05.self.out",
    header = F,
    col.names = c(
        "query", "subject", "pID", "length", "mismatches", "gapOpen",
        "qStart", "qEnd", "sStart", "sEnd", "evalue", "score"
    )
)

# since we blasted protein against themselves, we may remove reciprocal hits.
# Among two reciprocal hit, we will retain the ones with best score 
# (score are not always identical, for some reason)

# we select the best HSP for each subject-query pair (including reciprocal hits)
blastp <- removeReciprocal(blastp, removeSelf = F)

# we will also need to make all hits reciprocal for our procedure.
# we thus make an equivalent table with query and subject reversed
blastpR <- copy(blastp[query != subject])
blastpR[, c("query", "subject", "qStart", "qEnd", "sStart", "sEnd") :=
    .(subject, query, sStart, sEnd, qStart, qEnd)]

# we rbind this table with the original hits
blastp <- rbind(blastp, blastpR)
rm(blastpR)

# we creae a unique identifier for the query-subject pair
blastp[, pair := paste(query, subject)]


# --------------
# TE copies from the same community may have several protein-coding regions that are
# overlapping or adjacent on a given protein. We can combine these regions, so we
# get larger, non-overlapping, regions that should more likely represent longer
# ancestral elements


# for the procedure, we generate these tables:
# all copies involved in each community, we also retrieve the mrca of species involved
dt <- httHits[, .(copy = unique(c(ID.1, ID.2))), by = .(com, mrca)]
# this table is only needed to generate the following one:

# blastx hits for copies in each community. (This command takes a bit of time)
protHits <- dt[, blastx[query %in% copy], by = .(com, mrca)]
    #FOr each copy of dt, add coloumn with result about blastx
    #if several lines for a same copy in dt (when a copy is invovled in several com ad/or mrca),
        #the results of blastx will be duplicated
    #eg: 851298 invovled in 2 com (so 2 lines in dt)
        # 851298 has hits on 6 prot of blastx
        #--> 851298 has 2*6 lines in protHits

# we count the number of different communities that have copies
# matching proteins for each MRCA (within which hits were clustered).
comPerMRCA <- protHits[, length(unique(com)), by = mrca]

# if there is just one such community, there is nothing to do so we can ignore these
protHits <- protHits[mrca %in% comPerMRCA[V1 > 1L, mrca]]

# we combine protein regions that are distant by 10 amino
# acids or less on the same protein in a given community
protRegions <- combineRegions(protHits[, .(
    
    # as we combine protein regions only within communities,
    # we attach the community number to the protein name
    stri_c(com, subject, sep = " "), 
    sStart, sEnd
)],
distance = 10L
)

# see function combineHomologous() in HTvFunctions.R for details
protRegions <- combineHomologous(protRegions, blastp, protSeqs) 
#NOTE: that we use protSeqs or protSeqs[prot], it is the same
#NORE2: the com column is now called commID

fwrite(protRegions, "TEs/clustering/occ200ds05_protRegionsInCommunities_perClade.txt", sep='\t')

# We now compare communities of TEs in respect at protein regions -----------------------------------

# we will generate every possible pair of rows from the protRegions table in order
# to compare protein regions with each others and assess how they overlap
# protein regions are identified by row indices in the protRegion table

setorder(protRegions, commID)

# In Zhang et al. (2020), we generated all possible pairs of protein regions
# It takes too much RAM with all the hits we have here
# --> I change that part of the code to do only part we need (don't need to compare regions of different superF)

protHits = mutate(protHits, superF = sub(".*#", "", subject))
com_superF = select(protHits, c(com, superF)) %>% distinct()

# we will retain only the region pairs for community pairs that are relevant
# we thus create an identifier for community pairs
#Such as for pair of copies, identifier is com1*coef+com2
coefCom = 10^nchar(max(protRegions$commID))
# Note that com1 < com2 
stats[, comPair := com1 * coefCom + com2]

fwrite(stats, "TEs/clustering/occ200ds05_stats_perClade_save", sep='\t')
#NOTE If one continue with fread(stats_save),
#comPair will be read as an interger64. 
#Whereas it was a double before we saved it
#And we need a double, so make it a double again 
# (eiter by renning again the line that make comPair or by using as.double)

protRegions = mutate(protRegions, region = 1:nrow(protRegions))

dir.create("TEs/clustering/round2/protRegions_tmp/")
pairsOfRegions_perSuperF <- function(superFi){
    #Get only part of the table of superF
    protRegions_tmp = filter(protRegions, commID %in% filter(com_superF, superF==superFi)$com)

    setorder(protRegions_tmp, commID)
    nr <- nrow(protRegions_tmp)
    rows <- 1:(nr - 1L)
    
    #Function to add info to pairsOfRegions_tmp
    f_addInfo <- function(pairsOfRegions_tmp){
        # we retrieve the communities of these protein regions
        pairsOfRegions_tmp[, c("com1", "com2") := data.table(
            protRegions_tmp[region1_tmp, commID],
            protRegions_tmp[region2_tmp, commID]
        )]


        #we retrieve the number of row in the full protRegions
        pairsOfRegions_tmp[, c("region1", "region2") := data.table(
            protRegions_tmp[region1_tmp, region],
            protRegions_tmp[region2_tmp, region]
        )]
    
        #We don"t need the numerotation of regions specific to this tmp anymore
        pairsOfRegions_tmp = select(pairsOfRegions_tmp, -c(region1_tmp, region2_tmp))

        #We add compair
        # note that com1 is always <= com2 since we sorted the protRegions table by commID
        pairsOfRegions_tmp[, comPair := com1 * coefCom + com2] 

        # and remove pairs of regions we no longer need
        pairsOfRegions_tmp <- pairsOfRegions_tmp[comPair %in% stats$comPair]
    }
    
    # if protRegions_tmp is too big (too big to be contained in an integer),
    # data.table() will convert region1 & region2 to integer, which will fail (even if I use double for region1 & region2).
    # So I have to split the table. Here I split it in 3.
    #NOTE tmp1 & tmp2 are shorter than tmp3, but their pairsOfRegions will be bigger anyway
    if((nr*(nr-1)/2)>2^31){ 

        tmp1_rows = rows[1:floor(length(rows)/4)]
        tmp2_rows = rows[((floor(length(rows)/4))+1):(floor(length(rows)/2))]
        tmp3_rows = rows[((floor(length(rows)/2)+1)): length(rows)]

        f_pairsOfRegions_tmp <- function(tmp_rows){
            dt = data.table(region1_tmp = rep(tmp_rows, nr - tmp_rows), region2_tmp = unlist(lapply(X = tmp_rows[2]:(tmp_rows[length(tmp_rows)]+1), FUN = function(x) x:nr)))
        }
        pairsOfRegions_tmp1 <- f_pairsOfRegions_tmp(tmp1_rows)
        pairsOfRegions_tmp2 <- f_pairsOfRegions_tmp(tmp2_rows)
        pairsOfRegions_tmp3 <- f_pairsOfRegions_tmp(tmp3_rows)

        pairsOfRegions_tmp1 <- f_addInfo(pairsOfRegions_tmp1)
        pairsOfRegions_tmp2 <- f_addInfo(pairsOfRegions_tmp2)
        pairsOfRegions_tmp3 <- f_addInfo(pairsOfRegions_tmp3)
        pairsOfRegions_tmp = rbind(rbind(pairsOfRegions_tmp1, pairsOfRegions_tmp2), pairsOfRegions_tmp3)

    #if protRegions_tmp is not that big, no need to split it
    } else { 

        pairsOfRegions_tmp <- data.table( 
            region1_tmp = rep(rows, nr - rows),  
            region2_tmp = unlist(lapply( 
                X = 2:nr,
                FUN = function(x) x:nr
        )))

       pairsOfRegions_tmp <- f_addInfo(pairsOfRegions_tmp)
    }
    
    print(paste("done", superFi))

    #save in case
    fwrite(pairsOfRegions_tmp, paste0("TEs/clustering/round2/protRegions_tmp/", sub(".*[/]", "", superFi), ".txt"), sep='\t')

    return(pairsOfRegions_tmp)
}  

pairsOfRegions <- rbindlist(lapply(unique(com_superF$superF), pairsOfRegions_perSuperF))


gc() # reclaims some memory

# we create a table containing the coordinate of the protein
# regions between which we want to determine the homology
regionPairCoordinates <- pairsOfRegions[, data.table(
    comPair,
    protRegions[region1, .(prot, start, end)],
    protRegions[region2, .(prot2 = prot, start2 = start, end2 = end)]
)]

# and we determine the homology at each pair of protein regions
# see proteinOverlap() function in HTvFunctions.R
regionPairCoordinates <- proteinOverlap(
    regionPairCoordinates,
    blastp,
    protSeqs 
) 

fwrite(regionPairCoordinates, "regionPairCoordinates_perClade_save", sep='\t')


# for each pair of communities, we obtain statistics about protein homology
comPairStats <- regionPairCoordinates[, data.table(
    
    # length of the longest hsp between protein regions of the 2 communities
    maxOverlap = max(interAll), 
    
    # the total length of hsp between protein regions of the 2 communities
    sumOverlap = sum(interAll), 
   
    # the max pID of the HSP 
    maxID = max(pID, na.rm = T)
), 
by = comPair
]

fwrite(comPairStats, "comPairStats_perClade_save", sep='\t')


# we merge these statistics with our previous statistics on community pairs
mStats <- merge(stats, comPairStats, by = "comPair", all = T)

# if a community was not even present in the protRegions, the values are NAs.
# We put zeros (we leave NAs for com1 == com2 as such pairs are irrelevant)
mStats[
    com1 != com2 & is.na(maxOverlap),
    c("maxOverlap", "sumOverlap", "maxID") := 0L
]

#I should save mStats here
fwrite(mStats, "TEs/clustering/occ200ds05_mStats_perClade_save", sep='\t')

# we also obtain information per hit community regarding the proteins it covers
# (these pieces of information are not used afterwards, but it can be useful in other contextes)

# this first requires retrieving total protein lengths from the protein sequences
protRegions[, protLength := nchar(protSeqs[prot])]

comStats <- protRegions[, .(
    
    # number of different proteins (remaining after combineHomologous()),
    nProt = length(unique(prot)), 
    
    # length of longest protein, 
    protLength = max(protLength), 
    
    # length of longest protein region that is covered by TE copies
    longestRegion = max(end - start + 1)
    ), 
    
by = commID]

# we add these statistic to our table of community pairs
mStats[, c(
    "nProt1", "longestRegion1", "protLength1", # these columns will apply to com1
    "nProt2", "longestRegion2", "protLength2"  # same for com2
) 
:= data.table(
        comStats[match(com1, commID), .(nProt, longestRegion, protLength)],
        comStats[match(com2, commID), .(nProt, longestRegion, protLength)]
    )]

# we replace NAs by zeros
mStats[is.na(nProt1), c("nProt1", "longestRegion1", "protLength1") := 0L]
mStats[is.na(nProt2), c("nProt2", "longestRegion2", "protLength2") := 0L]



# we can finally evaluate criterion 1-----------
# In addition the the condition we already evaluated, criterion 1 is passed...
mStats[
    #if the longest aligned protein part ("maxOverlap") 
    #between copies of the 2 transfers is less than 100 aa.
    maxOverlap < 100 & 
    
    # and if here is no nucleotide identity between copies 
    #within clade (i.e, no self blastn hit),  
    valid == 0,    

    #Crit1 becomes TRUE
    crit1 := T
] 
#NOTE: rows that already had TRUE for crit1, keep TRUE

#save too
fwrite(mStats, "TEs/clustering/occ200ds05_mStats_perClade_save2", sep='\t')

# this leaves the possibility that communities with no nucleotide homology 
# and insufficient protein homology be different parts of the same TEs 
# that were retained by different clades after a single HTT.




# STEP FOUR, we apply complete linkage clustering to cluster communities in "hit groups"----------------------
# the approach is explained in Peccoud et al. 2017 PNAS
# we link communities that pass criteria 1 and 2 as being possibly the result of
# the same HTT. An HTT can only include communities that are all linked with each
# other (= complete-linkage clustering). Because linkage doesn't mean 2
# communities MUST reflect the same HTT, it could just mean a lack of information
# to distinguish HTTs (while two communities not passing the criteria cannot
# possibly be included in the same HTT). So we have to take decisions as to which
# group to from. Within each clade pair and super family, we want to group first:
# communities that cover very similar (or the same) proteins (high maxID)

# we thus place these pairs of communities on top, since we will regroup from top to bottom
mStats <- mStats[order(-maxID)]

# and apply our complete-linkage clustering
hitGroups <- mStats[, cLinkFromPairs(
    V1 = com1,
    V2 = com2,
    linked = crit2 & crit1
)]


# and we add the hitGroup indentifiers (integers) to our community pair table
mStats[, c("hitGroup1", "hitGroup2") :=
    data.table(
        hitGroups[as.character(com1)],
        hitGroups[as.character(com2)]
    )]


# we check that all pairs of communities belonging to the same group passed both criteria
mStats[hitGroup1 == hitGroup2 & com1!=com2, all(crit2 & crit1)] # this should return TRUE

# we save these statistics to disk
fwrite(mStats, "TEs/clustering/occ200ds05_comPairStats_perClade.txt")

# we assign each HTT hit to a hit group (new column added)
httHits[, hitGroup := hitGroups[as.character(com)]]

# any hit that is not in a hit group constitutes its own hit group
#I want hits of a same com in the same hitGroup
corres <- data.table(a = unique(httHits[is.na(hitGroup),]$com),
                     b = 1:length(unique(httHits[is.na(hitGroup),]$com)) + max(hitGroups))
httHits[is.na(hitGroup), hitGroup := corres[match(com, a), b]]

# we write results to disk
fwrite(httHits,"occ200HitGroup_perclade.txt", sep='\t')