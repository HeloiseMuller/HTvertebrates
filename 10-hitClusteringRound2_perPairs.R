library(data.table)
library(dplyr)
library(Biostrings)
library(stringi)
library(parallel)


path = "~/Project/"
setwd(path) 

source("HTvFunctions.R")

#Check how many copies we have
#This number will determine how to calculate unique identifiers when will compares hits of copies
copyIDs <-fread("TEs/clustering/dcMegablast/IDcopies_full.txt")
nCopies <- nrow(copyIDs)
print(nCopies)

coef = 10^nchar(nCopies)
#So unique identifiers for a pair of copies will be copy1*coef + copy2
#eg: ID of the hit of copy 1 vs copy 4598 will be 
    # 100004598 if nCopies = 9999 (coef = 4)
    # 100000004598 if nCopies = 2271890 (coef = 7)


#NOTE will use the same method latter in this script to make ID of pair of communities (comPair) 

#STEP ONE: compare com 2 by 2 to see whether they could result from the same HTT event
httHits <- fread("occ200ds05_comm_perPairs.txt") #Output script 9
httHits[, group := stri_c(sub(".*/", "",  rep_superF.1), ":", assembly.1, ":", assembly.2)]

# like in the previous stage, we split the hits by group and select the column we need for the clustering
hitList <- split(
    x = httHits[, .(
        hit = 1:.N,
        query = q,
        subject = s,
        pID = as.integer(pID * 1000L), # pID is again converted to integer for speed
        group,
        ucomm
    )],
    f = httHits$group
)

# we determine the number of hits per group. 
# No clustering need be done if there is just one community of hits per group
nCom <- sapply(
    X = hitList,
    FUN = function(x) length(unique(x$ucomm))
)

hitList <- hitList[nCom > 1L] 

# and we sort the the hit list by decreasing number of hits (for better CPU usage)
hitList <- hitList[order(sapply(
    X = hitList,
    FUN = nrow
),
decreasing = T
)]

## From here, we technically enter the script iterativeSecondClustering.R

# we retrieve super family names from the list of hits
superF <- splitToColumns(names(hitList), ":", columns=1)
nHits <- sapply(hitList, nrow)
# to process the super families with more hits first,
# we sum number of hits per superF (no matter group) and we order it
nHits <- sort(tapply(nHits, superF, sum), T)  

#Get outputs of the selfblast
blastFiles <- list.files(
    path = "TEs/clustering/selfBlastn/out",
    full.names = T, recursive=T
)

# we name the blast files elements with superfamily names (extracted from file names)
names(blastFiles) <- sub("/",":", sub(".IDs_selfBlastn.out", "", sub("TEs/clustering/selfBlastn/out/", "", blastFiles))) 

#  we "confront" every two communities of hits, we do this per group ------------------------------------------------------------
communityPairStats <- function(group) { 
    # we identify the self-blastn outputs we need, and import them
    superF = sub(":.*", "", group)
    assembly1 = splitToColumns(group, ":",columns=2)
    assembly2 = splitToColumns(group, ":",columns=3)
    
    #Read self blast of assembly 1
    if(file.size(blastFiles[paste0(assembly1,":",superF)])>0){
        blast1 <- fread(
            input = blastFiles[paste0(assembly1,":",superF)],
            header = F,
            select = 1:3, 
            col.names = c("query", "subject", "pID")
        )
    } else {
        blast1 = data.table("query"=numeric(), "subject"=numeric(), "pID"=numeric())
    }
    
    #Read self blast of assembly 2
    if(file.size(blastFiles[paste0(assembly2,":",superF)])>0){
        blast2 <- fread(
            input = blastFiles[paste0(assembly2,":",superF)],
            header = F,
            select = 1:3,
            col.names = c("query", "subject", "pID")
        )
    }  else {
        blast2 = data.table("query"=numeric(), "subject"=numeric(), "pID"=numeric())
    }
    
    blast = rbind(blast1, blast2)    
    
    
    if(nrow(blast)>0){
        blast[, pID := as.integer(pID * 1000L)] # we convert pIDs to integers as we did for the htt hits (to save memory)

        # but here we don't create a matrix of pIDs where subject and queries are rows and columns numbers,
        # such matrix would be too big for this round.
        # We create a unique subject-query identifier based on the fact that these are integers
        #The coefficient by which multiply query was calculated at the begginging of the script. 
        blast[, copyPair := query * coef + subject]  
    } else {
        blast$copyPair=numeric()
    }

    # we select the group of hits (one per group) for this superfamily
    hits <- hitList[[group]]
    

    # we sort these hits per community to ensure that the
    # community pairs are formed in the same "direction"
    # (com1 vs com2 and never the opposite)
    setorder(hits, ucomm)
 
    # the ids of copies involved in the hits to cluster
    uCopies <- hits[, unique(c(query, subject))]

    # we selected the TE hits involving these copies. We only
    # need the subject-query pair identifier and the blast pID
    selectedHits <- blast[query %in% uCopies & subject %in% uCopies, .(copyPair, pID)]

    # we add a copyPair that does not actually exist (with 0 pID)
    # at the end of this table. This trick will help us later
    #Important even if selectedHits is empty
    selectedHits <- rbind(selectedHits, data.table(copyPair = 0, pID = 0L))

    #We create a table that will analyse each pair of hits at once
    uHits = hits$hit
    pairs = data.table(
        hit1 = rep(uHits,(length(uHits)-1):0), 
        hit2 =  unlist(lapply(X = 2:length(uHits),  FUN = function(x) uHits[x:length(uHits)]))
    )
    

    # we add columns for the ids of copies involved in the 2 hits
    pairs[, c(
        "q1",                        # query ids (copies from clade A) for hit1
        "s1",                        # subject ids (copies from clade B)
        "com1",                      # community id
        "inter1",                    # the inter-clade identity (representing the HTT)
        "q2", "s2", "com2", "inter2" # and the equivalent for hit2
        ) 
        := data.table(
            hits[match(hit1,hit), .(query, subject, ucomm, pID)],
            hits[match(hit2,hit), .(query, subject, ucomm, pID)]
        )]

        # we create query-subject pair identifiers as we did for the self-blastn output
        # we must ensure that the left copy id is always lower than right one
        # (as it is the case for the copyPair identifiers in the self-blast files)
        #ATTENTION Here we don't swap the other columns. Make sure we don't care to know what pID belong to who
        pairs[q1 > q2, c("q1", "q2") := .(q2, q1)] # we swap these ids if it is not the case
        pairs[s1 > s2, c("s1", "s2") := .(s2, s1)]

        # we create the same pair identifiers we created for the blast results
        pairs[, c("qPair", "sPair") := .(
            q1 * coef + q2,
            s1 * coef + s2
        )]

        # so we can retrieve within-clade sequence identities.
        # The nomatch argument below makes it so that the pID of the last row
        # in selectedHits (which has pID 0) is returned where there is no match (no hit).
        # This speeds things up a little
        pairs[, c("intra1", "intra2") := .(
            selectedHits[match(qPair, copyPair, nomatch = .N), pID], #if no match, give last value of selectedHits, which is 0 since we added this line 
            selectedHits[match(sPair, copyPair, nomatch = .N), pID]
        )]

        # for hits of a copy with itself (not in the blast outputs), we define 100% identity
        # this is equivalent to what we did for the matrix' diagonal in the other script)
        pairs[q1 == q2, intra1 := 100000L]
        pairs[s1 == s2, intra2 := 100000L]

        pairs[, maxIntra := pmax(intra1, intra2)]
        #Here corrected mistake in Zhang et al. 2020: "<="," not "<"
        pairs[, connected := inter1 <= maxIntra |
                inter2 <= maxIntra]
         cat("*")
            

        # we return statistics for each community pair
        # (note that a community pair may involve the same community twice,
        # which can help to assess the degree to which hits within a community are connected)

       statsPerCom <- pairs[, .(
            allDiff = sum(!connected), # number of pairs of hits not passing the criterion
            diff = sum(maxIntra > 0L & !connected), # the same, but for hits whose copies have some degree of intra-clade identity
            valid = sum(maxIntra > 0L), # the number of pairs of hits whose copies have some degree of intra-clade identity,
            tot = .N # the total number of pairs of hits evaluated
        ), 
        by = .(com1, com2) # all this is done for each community pair
        ] 
        
    print(paste("done", group))

    statsPerCom
}

# we apply the function in parallel for the different batches
#much faster than in script 09 (took me less than 2 hours)
stats <- rbindlist(mclapply(
    X = names(hitList),
    FUN = communityPairStats,
    mc.cores = 20,
    mc.preschedule = F
))


# We now evaluate the connections between hit communities that may represent the same HTTs
# the number of hit pairs that pass criterion 1 between two communities
stats[, links := tot - allDiff] 

# criterion 1 needs to be passed by more than 5% of hit pairs
# to consider the communities as representing the same HTT.
stats[, crit1 := links / tot > 0.05]

fwrite(stats, "TEs/clustering/occ200ds05_comm_valid_perPairs.txt", sep='\t')
#If valid = 64845, it means that there are 64845 pairs of hits whose copies have some degree of intra-clade identity 
    #(out of the total number of hits evaluated = tot)
    #We will evaluate whether 2 communities invovle different part of a same TE only if valid == 0


# But we are not done with criterion 1 yet (see STEP THREE)

#MAKE figure S5 from Peccoud et al 2017
library(ggplot2)
#more precisely, draw in grey those of same com, and in black those of different com 
par(mfrow=c(1,2))
hist(stats[com1==com2,]$links/stats[com1==com2,]$tot, col=rgb(0.5,0.5,0.5), border=F, breaks=100, xlab = "Connectivity between hits of a same community", main="")
hist(stats[com1!=com2,]$links/stats[com1!=com2,]$tot, col=rgb(0,0,0), border=F, breaks=100, xlab = "Connectivity between hits of different communities", main="")
#I can't overlap both figures as in Peccoud et al. 2017 because y is not at the same scale

#########################################################################################

#STEP TWO: when clustering per pair of species, we don't need it (criterion 2)

#Because criterion 2 is about the age of the clade, but here we do clustering per pair of species

#########################################################################################

#STEP THREE

# The first lines being in common with 10-hitClusteringRound2_perClade.R, we removed them from this script. 
# Instead, we directly read the output of interest :
blastp <- fread(
    input = "TEs/blastp/involvedProtOCC200dS05.self.out",
    header = F,
    col.names = c(
        "query", "subject", "pID", "length", "mismatches", "gapOpen",
        "qStart", "qEnd", "sStart", "sEnd", "evalue", "score"
    )
)

#NOTE the folowing lines are also in common with 10-hitClusteringRound2_perClade.R :
blastp <- removeReciprocal(blastp, removeSelf = F)
blastpR <- copy(blastp[query != subject])
blastpR[, c("query", "subject", "qStart", "qEnd", "sStart", "sEnd") :=
    .(subject, query, sStart, sEnd, qStart, qEnd)]
blastp <- rbind(blastp, blastpR)
rm(blastpR)
blastp[, pair := paste(query, subject)]


#We still need to read blastx (generated in script 06-filterTEhits.R):
blastx <- fread(
    input = "TEs/blastx/all.copies.successiveBlastx.out",
    header = T,
    sep = "\t",
    col.names = c(
        "query", "subject", "pID", "length", "mismatch", "gapopen",
        "qStart", "qEnd", "sStart", "sEnd", "evalue", "score"
    )
)
blastx[, query := copyIDs[chmatch(query, copy), ID]]
blastx <- blastx[!is.na(query)]


#and also need to read protSeqs
protSeqs <- readAAStringSet(paste0(path, "/DB_clustered/findDubious/RepeatPeps.lib"))
names(protSeqs) <- splitToColumns(names(protSeqs), " ", "",  1) 


# --------------
# TE copies from the same community may have several protein-coding regions that are
# overlapping or adjacent on a given protein. We can combine these regions, so we
# get larger, non-overlapping, regions that should more likely represent longer
# ancestral elements

httHits <- fread("occ200ds05_comm_perPairs.txt")

# for the procedure, we generate these tables:
httHits = mutate(httHits, pair = ifelse(assembly.1<assembly.2, paste0(assembly.1, "-", assembly.2), paste0(assembly.2, "-", assembly.1))) 
dt <- httHits[, .(copy = unique(c(q, s))), by = .(ucomm, pair)]
# this table is only needed to generate the following one:

# blastx hits for copies in each community. (This command takes a bit of time)
       
protHits <- dt[, blastx[query %in% copy], by = .(ucomm, pair)] 
    #for each line of dt, we add info about blastx
    #if this copy hit several prot during blastx (we did a max of 5 successive blastx), one line per hit --> nrow(protHits)>nrow(dt)

# we count the number of different communities that have copies
# matching proteins for each MRCA (within which hits were clustered).
comPerMRCA <- protHits[, length(unique(ucomm)), by = pair] 
    #e.g. the pair GCA_905171765.1-GCA_905171775.1 has 154 communities(column V1)

# if there is just one such community, there is nothing to do so we can ignore these (can't be conected to another community if alone)
# so we keep only pairs that have more than 1 community
protHits <- protHits[pair %in% comPerMRCA[V1 > 1L, pair]] 

#Save what communities correspond at what pair
comm_pair = select(protHits, c(ucomm, pair)) %>% distinct() 
fwrite(comm_pair, "occ200ds05_comm_perPairs.txt", sep='\t')
#NOTE: Pairs that have only one comm where not saved in this table

# we combine protein regions that are distant by 10 amino
# acids or less on the same protein in a given community
protRegions <- combineRegions(protHits[, .(
    
    # as we combine protein regions only within communities,
    # we attach the community number to the protein name
    stri_c(ucomm, subject, sep = " "), 
    sStart, sEnd
)],
distance = 10L
)

# see function combineHomologous() in HTvFunctions.R for details
protRegions <- combineHomologous(protRegions, blastp, protSeqs) 
fwrite(protRegions,"TEs/clustering/occ200ds05_protRegionsInCommunities.txt", sep='\t')


# We now compare communities of TEs in respect at protein regions -----------------------------------

# we will generate every possible pair of rows from the protRegions table in order
# to compare protein regions with each others and assess how they overlap
# protein regions are identified by row indices in the protRegion table

setorder(protRegions, commID)
nr <- nrow(protRegions)
rows <- 1:(nr - 1L)

setorder(protRegions, commID)
protRegions = mutate(protRegions, region = 1:nrow(protRegions))
pairsOfRegions_perPair <- function(v.pair){
    #Get only part of the table of that pair
    protRegions_tmp = filter(protRegions, commID %in% filter(comm_pair, pair==v.pair)$ucomm )

    setorder(protRegions_tmp, commID)
    nr <- nrow(protRegions_tmp)
    rows <- 1:(nr - 1L)
    
    pairsOfRegions_tmp <- data.table( 
    region1_tmp = rep(rows, nr - rows), 
    region2_tmp = unlist(lapply( 
        X = 2:nr,
        FUN = function(x) x:nr
    )))
    
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
    
}  

pairsOfRegions <- rbindlist(lapply(unique(comm_pair$pair), pairsOfRegions_perPair))


################################"

# we will retain only the region pairs for community pairs that are relevant
# we thus create an identifier for community pairs
# note that com1 is always <= com2 since we sorted the protRegions table by commID
#Like for coef on the number of copies, calculate the coef to use for comPair
coefCom = 10^nchar(max(pairsOfRegions$com2))
pairsOfRegions[, comPair := com1 * coefCom + com2] 

# and similar identifiers for our table containing statistics on community pairs.
# Note that com1 < com2 is both cases
stats[, comPair :=  com1 * coefCom + com2] 

# and remove pairs of regions we no longer need
#This is important, oterhwise I will have NA in column valid in mStats.
#NOTE: or stats, we did not compare com from different superF. So if NA here, it means they belong to different superF
pairsOfRegions <- pairsOfRegions[comPair %in% stats$comPair]


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
    #if no blastp hit between the pair of prot --> NA for coords and 0 for pID, length, interAll
regionPairCoordinates <- proteinOverlap(
    regionPairCoordinates,
    blastp,
    protSeqs
)

#Can save this intermediate file
fwrite(regionPairCoordinates, "regionPairCoordinates_perPairs_save", sep='\t')

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

#Can also save this intermediate
fwrite(comPairStats, "comPairStats_perPairs_save", sep='\t')

# we merge these statistics with our previous statistics on community pairs
#NOTE: comPair is so big that it will read it as integer64.
# But I want double since comPair of stats is a double. I need the same type to merge.
comPairStats <- fread("comPairStats_perPairs_save", integer64 = "double") 

# we merge these statistics with our previous statistics on community pairs
mStats <- merge(stats, comPairStats, by = "comPair", all = T)

#NOTE I have plenty of NA. 
    #NA due to absence in table pairsOfregion (so NA in column maxOverlap): 
    #this is when there is only one region for this pair of com, so there was no region to compare homolgy with

# We replaced these NA by zeros 
# Expect when com1 == com2 as such pairs are irrelevant
mStats[
    com1 != com2 & is.na(maxOverlap),
    c("maxOverlap", "sumOverlap", "maxID") := 0L
]
 
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
   
   #crit1 take value TRUE
   crit1 := T
] 
#To note, other lines of mStats already had the value TRUE because of the criteria 5%
#crit1 in thus a OR. Either pass the 5% criteria and/or pass the protein criteria

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
# DO I NEED THAT STEP?

# we thus place these pairs of communities on top, since we will regroup from top to bottom
mStats <- mStats[order(-maxID)]

# and apply our complete-linkage clustering
hitGroups <- mStats[, cLinkFromPairs(
    V1 = com1,
    V2 = com2,
    linked = crit1
)]


# and we add the hitGroup indentifiers (integers) to our community pair table
mStats[, c("hitGroup1", "hitGroup2") :=
    data.table(
        hitGroups[as.character(com1)],
        hitGroups[as.character(com2)]
    )]


# we check that all pairs of communities belonging to the same group passed both criteria
mStats[hitGroup1 == hitGroup2 & com1!=com2, all( crit1)] # this should return TRUE

#com1=com2 are useless, remove them:
mStats = mStats[com1 != com2,]

# we save these statistics to disk
fwrite(mStats, "TEs/clustering/occ200ds05_comPairStats_perPairs.txt", sep='\t')

# we assign each HTT hit to a hit group (new column added)
httHits[, hitGroup := hitGroups[as.character(ucomm)]]

# any hit that is not in a hit group constitutes its own hit group
# (these small hit groups will be removed by our checkpoints in script 11)

#I want hits of a same ucomm in the same hitGroup
corres <- data.table(a = unique(httHits[is.na(hitGroup),]$ucomm),
    b = 1:length(unique(httHits[is.na(hitGroup),]$ucomm)) + max(hitGroups))
httHits[is.na(hitGroup), hitGroup := corres[match(ucomm, a), b]]

# we write results to disk
fwrite(httHits, "occ200HitGroup_perPairs.txt", sep='\t')