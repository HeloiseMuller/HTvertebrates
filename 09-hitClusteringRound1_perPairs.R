library(data.table)
library(Biostrings)
library(dplyr)
library(stringi)
library(stringr)
library(igraph) #for function graph_from_data_frame
library(parallel)


path = "~/Project/"
setwd(path)

source("HTvFunctions.R")

httHits <- fread("TEs/clustering/dcMegablast/occ200ds05_dcMegablast.txt")
copyIDs <-fread("TEs/clustering/dcMegablast/IDcopies_full.txt")

nCPUs <- 30


########################

#We add ID to copy name (each copy has a numeric ID)
httHits[, c("q", "s") := .(copyIDs[chmatch(copie1, copy), ID], copyIDs[chmatch(copie2, copy), ID])]

#We keep just the colums we need for clustering 
conv = select(httHits, c(hit, pID, q, s, assembly.1, assembly.2, rep_superF.1, divTime, mrca)) %>% 
    mutate(., pID = as.integer(pID*1000L))

#Make sure that for each pair of species, it is always the same one on the right
#NOTE I looked at superF1 only because at this step superF1=superF2
conv = mutate(conv, group = ifelse(assembly.1<assembly.2, 
    paste0(assembly.1,"-",assembly.2, "-", sub(".*/", "", rep_superF.1)), 
    paste0(assembly.2,"-",assembly.1, "-",  sub(".*/", "", rep_superF.1))))

# we split these hits according to "groups"
hitList <- split(conv, conv$group)
hitList <- hitList[order(-sapply(hitList, nrow))]

# if there is just one hit in a group, no clustering needs to be done so we can remove the group
hitList <- hitList[sapply(hitList, nrow) > 1L]

groups <- names(hitList)
length(groups) #The function hitCommunities will happens on that many groups

groupList <- data.table(
    group = names(hitList),

    # number of hits to cluster per "group"
    nHits = sapply(hitList, nrow),
    superF = sapply(hitList, function(x){unique(x$rep_superF1)})
)

#Intra species and intra superF:
files_intra <- list.files("TEs/clustering/selfBlastn_perPairs/out/", recursive = T, pattern="_selfBlastn.out")
files_intra_dt <- data.frame(assembly = sub("assembly_", "", sub("/.*", "", files_intra)), 
    superF = str_match(files_intra, "/\\s*(.*?)\\s*.IDs")[,2], 
    file = paste0("TEs/clustering/selfBlastn_perPairs/out/", files_intra))

blast <- rbindlist(apply(X = files_intra_dt, 1, FUN = function(x) {
    b = fread(
        input = as.character(x['file']),
        header = F,
        select = c(1,2,3), 
        col.names = c("query", "subject", "pID"),
        sep = "\t"
    )
    #add info about assembly and superfamilly
   b = mutate(b, assembly = unique(x['assembly']), superF = unique(x['superF']))
   if(nrow(b)>0){
    return(b)
   }
}))

    #same as for httHits
blast = mutate(blast, pID = as.integer(pID*1000L))

### Function from script iterativeFirstClustering.R
hitCommunities <- function(group) {

    # we first retrieve the hits to cluster. We copy the data table to
    # avoid some side effect related to how data.table functions work
    hits <- copy(hitList[[group]])

    # we obtain the different copies (integer ids) in these hits
    ucopies <- hits[, sort(unique(c(q, s)))]

    # we convert these to smaller integers so that we can make matrices of blast pID
    # between copies, where a copy id will be a row/column index in the matrix
    hits[, c("qid", "sid") := .(
        match(q, ucopies),
        match(s, ucopies)
    )]

    # we thus need to use the same ids for the copies in the blast output (will be NA for copies of other groups)
    blast[, c("qid", "sid") := .(
        match(query, ucopies),
        match(subject, ucopies)
    )]

    # we extract the hits involving these copies
    selectedHits <- blast[!is.na(qid) & !is.na(sid)]

    # we make a matrix of pID for each copy pair (with pID zero by default, when there is no hit)
    pIDmatrix <- matrix(0L, length(ucopies), length(ucopies))

    # the diagonal is set to 100% pID (*1000), to represent perfect pID for a copy with itself
    diag(pIDmatrix) <- 100000L #I initially wrote 1000L...

    # we fill the matrix, (both semi matrices)
    pIDmatrix[cbind(selectedHits$qid, selectedHits$sid)] <- selectedHits$pID
    pIDmatrix[cbind(selectedHits$sid, selectedHits$qid)] <- selectedHits$pID

    # we cannot always cluster all the hits at once due to the max size of vectors in R (and RAM required)
    # so we get the number of hits to cluster
    nHits <- nrow(hits)

    # and determine the number of batches of hits, since we will have to make all possible pairs of hits
    nBatches <- ceiling(nHits^2 / 2^28)

    # in the following, a hit corresponds to a row index in the "hits" table
    # we split the hits into several batches 
    hitBatches <- splitEqual(1:(nHits - 1L), n = nBatches)


    # we "connect" hits 2 by 2 according to criterion 1, with this function:
    criterion_1 <- function(batch) {

        # for a batch of hits, we make all possible pairs of hits.
        # The left-hand hit (hit1) is from the batch, and the right-hand
        # hit (hit2) includes all other hits (including other batches)
        # this ensures that, over all batches, all possible pairs will be made
        pairs <- data.table(
            hit1 = rep(batch, nHits - batch),
            hit2 = unlist(x = lapply(
                X = batch[1]:max(batch) + 1L,
                FUN = function(hit) hit:nHits
            ))
        )

        # we retrieve the ids of copies involved in the 2 hits (2 per clade)
        pairs[, c("q1", "s1", "q2", "s2") := data.table(
            hits[hit1, .(qid, sid)],
            hits[hit2, .(qid, sid)]
        )]

        # and retrieve the blast pIDs (percentage identity) within each clade (intra) of the 2 hits (inter)
        pairs[, c(
            "inter1", # between-clade pID, representing the HTT (for hit1)
            "inter2", # same for the right-hand hit1
            "intra1", # pID of copies within the left-clade (clade A)
            "intra2"
        ) # and for the right-clade
        := data.table(
                hits[hit1, pID],
                hits[hit2, pID],
                pIDmatrix[cbind(q1, q2)],
                pIDmatrix[cbind(s1, s2)]
            )]

        # to apply criterion 1, we get the highest within clade identity dans
        pairs[, maxIntra := pmax(intra1, intra2)]
        cat("*") # progress indicator (this can be long)

        # and we finally apply criterion 1 to "connect" the hits
        # in effect, we return pairs of hits where the best intra-clade
        # identity is higher than at least inter-clade identity (that of hits)
        pairs[inter1 <= maxIntra | inter2 <= maxIntra, maxIntra, .(hit1, hit2)] #return only those that pass criteria
    }

    # we apply the function for batches of hit pairs and concatenate the results
    pairs <- rbindlist(lapply(hitBatches, criterion_1))

    #Check how many hits are returned. A hit returned has at least 1 connection with another hit
    nbHitsReturned <- length(unique(c(pairs$hit1, pairs$hit2)))

    
    # we now do the clustering if there are "connected hits"
    #NOTE: the algo cluster_fast_greedy has a weird behavious:
        #if everything is connected with everything, it  will always put the last node in another com 
        #(whereas everything should be together)
   
   #If pairs not empty, it means at least one connection
   if (nrow(pairs) > 0 ) {
   
        #and if not everything connected to everything
        if ((nbHitsReturned*nbHitsReturned-nbHitsReturned)/2 != nrow(pairs) ) { 
            cls <- graph_from_data_frame(pairs, directed = F)

            # we cluster into "communities" done with clauset et al. algorithm
            cls <- data.tableFromCommunities(cluster_fast_greedy(cls))
        
        #if everything is connected (in pairs), do not use cluster_fast_greedy
        #instead, put everything in the same com
        } else {
            cls <- data.table(member = unique(c(pairs$hit1, pairs$hit2)), community = 1)
        }

        cat("-") # to monitor progress

        # we attribute hits to "communities" named by integer numbers
        com <- integer(nrow(hits))
        com[cls$member] <- cls$community
        # so com[x] will give the community of hit "x" (which is an integer)

        # any hit that is not in a community now forms its own community
        # those are the hits not returned at all by pairs, because they don't have a single connection
        com[com == 0L] <- 1:(sum(com == 0L)) + max(com)

        # we add the community column to the original table of htt hits
        hits[, comm := com]
        
    #if pairs is empty, it means that nothing is connected. So give one com per hit
    } else {
        hits[, comm := 1:.N]
    }

    # we write results to disk for safety (the function also returns the results)
    fwrite(hits[, .(hit, comm, group)], stri_c("TEs/clustering/round_perPairs/", group, ".groups.txt"), sep='\t')

    hits[, .(hit, comm, group)]
}

mrca = fread("Busco/pairsALL.txt")

# we apply the above function to the different groups of the super family, in parallel
#NOTE: the following can take a couple days, so one might want to run it in background 
res <- mclapply(
    X = groups,
    FUN = hitCommunities,
    mc.cores = nCPUs,
    mc.preschedule = F
)

res <- rbindlist(res)


# we generate unique integer community ids (accross all groups)
res[, ucomm := toInteger(string = stri_c(comm,
    group,
    sep = " "
))]

#Add info on divTime & mrca
res <- select(conv, c(group, mrca, divTime)) %>% distinct() %>% right_join(., res)
 
# we add these identifiers to the HTT hit table
httHits <- left_join(httHits, select(res, c(hit, ucomm))) 

#When only 1 hit in a superfamily and pair, we didn't included it in the clustering step
#So com = NA for these hits
#We now give them a unique community numbers:
httHits[is.na(ucomm), ucomm := 1:.N + max(httHits$ucomm, na.rm=T)]

# write the new HTT hit table to disk
fwrite(httHits, "occ200ds05_comm_perPairs.txt", sep='\t')