## %######################################################%##
#                                                          #
####          This script shuffles species to           ####
####       generate random species pairs involved       ####
####            in transfers, for each super            ####
####       family, avoiding generating "illegal         ####
####   pairs" for which no HTT can be infered, due to   ####
####         our requirement that species must          ####
####             be divergent by 2*120 My.              ####
#                                                          #
## %######################################################%##

require(ape)
require(data.table)
require(dplyr)
library(stringr)
library(parallel)
library(stringi)
require(RColorBrewer)
library(matrixStats)

## Variables to set :
nCPUs <- 20

#Pick one, comment the others
included <- "all"
#included <- "Chordata"
#included <- "Insects"


path = "~/Project/"
setwd(path) 

source("HTvFunctions.R")

# a second optional argument is a node number to permute only species within this clade. 
# Defaults to the basal node to permute all species
tree <- read.tree("datasetTree.nwk")
node <- which.max(node.depth(tree))

# a third optional argument is the number of replicates we want, defaults to 1000
n <-1000

# we compute the number of legal permutation we need per job launched in parallel
maxn <- ceiling(n / nCPUs)

#######################################

###  Metadata ##
meta <- fread("metadata.tbl")
meta <- mutate(meta, species = gsub(" ", "_", species))
meta$group <- str_to_title(meta$group)
meta[, sp :=  chmatch(species, tree$tip.label)]

#Make groups at same taxonomical scale as other studies
#N.B. We call this column "order" because it corresponds to taxo order for insects (but not for others)
meta = mutate(meta, order = 
    ifelse(group == "Archosauria" & habitat == "aquatic", "Crocodylidae",
    ifelse(group == "Archosauria" & habitat == "terrestrial", "Aves",
    ifelse(group == "Chelicerata" & habitat == "terrestrial", "Arachnida",
    ifelse(species=="Nymphon_striatum", "Pycnogonida",
    ifelse(group == "Chelicerata" & habitat == "aquatic" & species!="Nymphon_striatum", "Limulidae",
    ifelse(dbBusco=="lepidoptera_odb10", "Lepidoptera",
    ifelse(group=="Amphiesmenoptera" & dbBusco=="endopterygota_odb10", "Trichoptera",
    ifelse(group %in% c("Syrphoidea", "Nematocera", "Ephydroidea"), "Diptera",
    ifelse(species %in% c("Baetis_rhodani", "Ephemera_danica", "Cloeon_dipterum"), "Ephemeroptera",
    ifelse(species %in% c("Pantala_flavescens", "Rhinocypha_anisoptera", "Calopteryx_splendens", "Ladona_fulva", "Ischnura_elegans"), "Odonata",
    ifelse(species %in% c("Clitarchus_hookeri", "Timema_cristinae"), "Phasmatodea",
    ifelse(species %in% c("Vandiemenella_viatica", "Gryllus_bimaculatus", "Laupala_kohalensis", "Apteronemobius_asahinai"), "Orthoptera",
    ifelse(species %in% c("Coptotermes_formosanus", "Cryptotermes_secundus", "Blattella_germanica"), "Blattodea",
    ifelse(species=="Anisolabis_maritima", "Dermaptera", group)))))))))))))))

#######################################

# STEP ONE, we select the hits and delineate the scope of permutations -------------------------------------------

# we import the hit representing htt
retainedHits <- fread("HTThitsAssessed_perClade_retained")

#If not all species are included, get the subset of meta of interest
if(included=="Chordata"){
    meta <- meta[phylum=="Chordata",]
} else if(included=="Insects") {
    meta <- meta[phylum=="Arthropoda" & !group %in% c("Chelicerata", "Crustacea"),]
}

#The following will be a subset only if meta is a subset:
retainedHits <- filter(retainedHits, species.1 %in% meta$species & species.2 %in%  meta$species) %>% data.table
tree <- keep.tip(tree, tip = meta$sp) 
node <- which.max(node.depth(tree))

#Update number sp in case changed when took subset of tree
meta[, sp :=  chmatch(species, tree$tip.label)]

setorder(meta, sp)
taxon <- meta$order

# as several species are involved in a hit group even though this would just represent one transfer, 
# we only retain the two species for the best hit (best pID) per hit group. We also replace their 
# names with tip numbers of the tree, as integer numbers speed things up and save memory considerably

# we put the best hits on top 
setorder(retainedHits, -pID)

# we extract the ones we want (note that we only consider "independent" transfers)
# and we replace species names with tip numbers at the same time
HTThits <- retainedHits[independent == T, .(
    sp1 = chmatch(species.1[1], tree$tip.label),
    sp2 = chmatch(species.2[1], tree$tip.label),
    repl = 0L
    # repl is a replicate number (= 0 for original, non-permuted species)
), by = .(hitGroup, superfamily)]


#I regroup some superfamilies
HTThits = mutate(HTThits, superF = ifelse(str_detect(superfamily, "TcMar"),"DNA/Mariner",
    ifelse(str_detect(superfamily, "hAT"), "DNA/hAT",
    ifelse(str_detect(superfamily, "Sola"), "DNA/Sola",
    ifelse(str_detect(superfamily, "CMC"), "DNA/CMC",
    ifelse(str_detect(superfamily, "Crypton"), "DNA/Crypton",
    ifelse(str_detect(superfamily, "MULE"), "DNA/MULE",
    ifelse(str_detect(superfamily, "PIF"), "DNA/PIF-Harbinger",
    ifelse(str_detect(superfamily, "Ginger"), "DNA/Ginger",
    ifelse(str_detect(superfamily, "LTR/ERV"), "LTR/ERV",
    ifelse(str_detect(superfamily, "LINE/RTE") | str_detect(superfamily, "Proto2"), "LINE/RTE",
    ifelse(str_detect(superfamily, "LINE/L1"), "LINE/L1",
    ifelse(str_detect(superfamily, "LINE/R2") | str_detect(superfamily, "LINE/CRE") | str_detect(superfamily, "R4"), "LINE/R2",
    ifelse(str_detect(superfamily, "LINE/L2") | str_detect(superfamily, "LINE/CR1"), "LINE/Jockey",
    ifelse(str_detect(superfamily, "LINE/R1"), "LINE/I",
    superfamily
)))))))))))))))


# we determine the scope of permutations. Ideally, we permute species involved in htt within a TE super family, 
# but some contain too few HTT for this to be meaningful, 
# and others contain too many, which makes impossible to obtain only "legal" permutations
# to make our choices, we compute the number of independent transfers per super family, 
nTr <- HTThits[, .N, by = superF]

# we add a column to denote the TE class, as we pool super families 
# that are involved in less than 20 transfers within classes
nTr[, class := ifelse(str_detect(superF, "DNA/") | superF=="RC/Helitron",
    "DNA",
    "RNA"
)]

nTr[, combined := ifelse(test = N < 20, 
                         yes = paste("other", class, sep = " "),  #for pooled superfamilies, we use "other" + the TE class
                         no = superF)]

# we replace names of underrepresented super families with the combined names
HTThits[, superfamily2 := nTr[chmatch(HTThits$superF, superF), combined]]

# we recompute the numbers of transfers with the pooled super families
nTr <- nTr[, .(N = sum(N)), by = combined]

# as there may be too many hits per superfamily to obtain only "legal" permutations, 
# we will split certain super families in batches of hits (when they encompass more than 120 transfers)
# We add a column for the number of time a superfamily has been seen in transfers
HTThits[, occ := occurrences(superF)]

# we compute the max number of transfers we allow per superfamily batch
nTr[, maxi := N / ceiling(N / 121)]

# which we transfer to the hits table
#HTThits[, maxi := nTr[match(HTThits$superF, combined), maxi]]
HTThits[, maxi := nTr[match(HTThits$superfamily2, combined), maxi]]

# for mariners, we need even smaller batches of ≤ 61 transfers 
# (for some reason, it is harder to obtain legal permutations in these TEs)
HTThits[superfamily2 == "DNA/Mariner", maxi := 61L]

#I do the same thing with DNA/hAT since they are about as many as Mariner
HTThits[superfamily2 == "DNA/hAT", maxi := 61L]

HTThits[, batch := ceiling(occ / maxi)]

# we split the hits by superfamily and batch
# a "hit" will be a pair of species associated with a permutation number
hitList <- split(x = HTThits[, .(sp1, sp2, repl)],
                 f = list(HTThits$superfamily2, HTThits$batch),
                 drop = T) 

# we get the total number of species
nSpecies <- length(tree$tip.label)

# this matrix will contain the permuted "replacement species" during the work, one set of permuted species per column.
newsp <- matrix(rep(1:nSpecies, 10^5), ncol = 10^5)
# Note that we anticipate 10^5 permutations per batch, although we keep much fewer. 
# This is because many permutations may be "illegal".
# The row indices of this matrix = the original species (integer ids)

# we determine which simulated transfers (those with permuted species) are "legal"-------------------------------
# we create the matrix of divergence times
divTimeMat <- cophenetic(tree)

# this vector will be TRUE for species pairs that are "too close" to be involved in HTT 
# (a species here is a row/column index of the logical matrix)
#tooClose <- divTimeMat < 80

#Other illiegal pairs are the ones for which dS busco <0.81
mrca_fitlerB = readRDS("mrca_FilterB.Rds")
pairsALL = fread("Busco_new/pairsALL.txt")


pairsALL <- mutate(pairsALL, tooClose = ifelse(mrca %in% mrca_fitlerB | divTime<80, TRUE, FALSE ))

#Only keep pairs that invovled species included in this analysis
pairsALL <- filter(pairsALL, sp1 %in% meta$species & sp2 %in% meta$species)

#write sp1 & sp2 in both sense
pairsALL2 = select(pairsALL, c(sp2, sp1, tooClose)) %>%
    rename(., sp1=sp2, sp2=sp1) %>% 
    rbind(., select(pairsALL, c(sp1, sp2, tooClose))) %>% distinct()

tooClose <- tapply(pairsALL2$tooClose, pairsALL2[,1:2], FUN=print)
names(dimnames(tooClose)) <- NULL

# we retrieve the species we will permute (tip numbers of the tree, for the clade/node we focus on)
toShuffle <- tipsForNode(tree, node) 

#######################################

# STEP TWO, we permute species for htts (hits) of a batch ----------------------------------------------------------
shuffleSpecies <- function(hits, superfamily) {
    # we print progress, which is the only use of the superfamily argument
    print(superfamily)
    
    # to speeds things up, we generate 10^5 permutations in a row as the vast
    # majority lead to illegal transfers. This allows taking advantage of
    # vectorised functions after that.
    # for that, we need the number of hits
    nHits <- nrow(hits)
 
    # we replicate the hits 10^5 times, but incrementing species ids at each replicate (by the total number of species = l)
    sp1 <- rep(0:(10^5 - 1), each = nHits) * nSpecies + hits$sp1
    sp2 <- rep(0:(10^5 - 1), each = nHits) * nSpecies + hits$sp2
    
    # this function performs the permutations in parallel. job is a simple integer identifier
    shuffleWork <- function(job) {

        # we prepare a matrix of "legal permutations" 
        # it is the same format as the newsp matrix, where row numbers = species ids, 
        # columns are different permutations and values are replacement species
        legalPermutations <- NULL

        # indicator to tell when to stop
        g <- 0
         
        # until we obtain maxn legal permutations:
        repeat {
           
            # this creates a matrix of permuted species, with 10^5 columns (each is a vector of shuffled species)
            # this is the workhorse function, all the rest is result handling
            newsp[toShuffle, ] <- replicate(10^5, sample(toShuffle))

            # we replace original species by the sampled ones in the transfers
            # the left-hand column is the original species, the right-hand one the replacement one
            newPairs <- cbind(newsp[sp1], newsp[sp2])

            # we compute the number of illegal transfers per permutation (we use a matrix to quickly count them via colSums)
            illegal <- colSums(matrix(tooClose[newPairs], nrow = nHits))
            
            if (any(illegal == 0L)) {

                # we extract the columns corresponding to permutations with no illegal transfer and add them to the retained permutations
                legalPermutations <- cbind(legalPermutations, newsp[, illegal == 0])
                g <- ncol(legalPermutations)

                # progress indicator
                cat(".")
            }
            
            # we exit once we have the number of legal permutations we want
            if (g >= maxn) {
                break
            }
        }

    
        # we unrolls the matrix of legal permutations (useful to replace original species with the sampled ones),
        newSp <- as.vector(legalPermutations[, 1:maxn])
        #  the index in this vector will be the original species identifier, and its value is the replacement species.
        # We do not retain more than maxn legal permutations (there may actually be up to 10^5 if there were few hits).

        # we generate a permutation id number
        repl <- rep(1:maxn, each = nHits)
        
        # and we make a table of hits involving the permuted species
        # we replicate the transfers maxn times but incrementing species numbers, 
            
        randomHits <- data.table(
            sp1 = rep(hits$sp1, maxn) + (repl - 1L) * nSpecies,
            sp2 = rep(hits$sp2, maxn) + (repl - 1) * nSpecies,
            repl

        )
        
        # we can now easily replace original species with the sampled ones.
        # This is similar to what we did to create the newPairs matrix, except this time we use a data table
        # we also use the job id to generate final permutation identifiers, which will have to differ between jobs
        randomHits[, c("sp1", "sp2", "repl") := .(newSp[sp1], newSp[sp2], repl + job * maxn)]
        
        randomHits
    }
    
    # we apply the above function to batches of hits in parallel
    randomHits <- mclapply(1:nCPUs - 1L,
        shuffleWork,
        mc.cores = nCPUs,
        mc.preschedule = F
    )
    
    #and we stack the results in a single table
    randomHits = rbindlist(randomHits)
    
    # we replace replication numbers with smaller number that do not exceed the number
    # of permutation that was specified
    randomHits[, repl := match(repl, unique(repl))]
    
    randomHits
}

# we apply the permutations to superfamilies successively
randomHTTs <- Map(shuffleSpecies, hitList, names(hitList))

# we stack randomized hits with read ones (that we can differentiate since their replication number is 0)
randomHTTs = Map(rbind, hitList, randomHTTs)

dir.create("permutations")
saveRDS(randomHTTs, file = stri_c("permutations/allPermutations_", included, "_Node.", node, ".RDS"))

#######################################

# STEP THREE, contrasting the number of HTT per taxa in real and randomised hits.---------------------------------------

# we import the permutations generated above (= pairs of species representing real and simulated HTT)
randomHTTs <- readRDS("permutations/allPermutations_", included, "_Node.", node, ".RDS")

# this function counts htt per taxon in a given TE super family (or batch within a larger superfamily) 
# to obtain statistics (means, quantiles) for number of simulated transfers for each taxon of interest
countHTTs <- function(superfamily) {

    # we "unfold" the species pairs in htts, separating species 1 and 2, duplicating HTTs. 
    # Effectively, if a hit involves 2 species of the same taxon, two HTTs will be counted for the taxon
    unFolded <- randomHTTs[[superfamily]][, .(sp = c(sp1, sp2), repl = c(repl, repl))]
   
    # we add a column denoting the taxon of each species
    unFolded[, taxon := taxon[sp]]

    # we count the numbers of transfers involving each taxon in each permutation (note that repl == 0 for real HTT events)
    counts <- unFolded[, .(N = .N / 2), by = .(taxon, repl)]
    
    # we retrieve the number of permutations
    nRepl <- max(randomHTTs[[superfamily]]$repl)
    
    # we need to add 0 counts for taxa not involved in HTT in each replicate
    # since they are simply missing (if we didn't do this, the stats would be biased upwards)
    # for this we need to know the taxa involved
    uTaxa <- unique(taxon)
    
    # ad we create a table indicating 0 transfer for each taxon in each replicate
    missingTaxa <- data.table(
        taxon = rep(uTaxa, nRepl),
        repl = rep(1:nRepl, each = length(uTaxa)),
        N = 0
    )
    
    # we place these rows at the end of our previous table
    counts <- rbind(counts, missingTaxa)

    # and we remove rows of 0 counts for taxa that were already present
    counts <- counts[!duplicated(data.table(taxon, repl))]
    
    # we now obtain the statistics we need, per taxon
    stats <- counts[repl > 0L, .(
        simulated = mean(N),       # mean number of transfer over permutations
        minNumber = min(N),        # min number '''
        maxNumber = max(N),        # max number '''
        qLc = quantile(N, 0.005),  # and the various quantiles to test for significance
        qRc = quantile(N, 0.995),
        qL5 = quantile(N, 0.025),
        qR5 = quantile(N, 0.975),
        superfamily
    ), by = taxon]

    # and we add number of real number of HTT events to these stats
    merge(stats, counts[repl == 0L, .(taxon, observed = N)], by = "taxon", all = T)
}


# we apply the function to all super families (batches)
permutationStats <- lapply(names(randomHTTs), countHTTs)

# we stack results in a single table
permutationStats <- rbindlist(permutationStats)

# NAs must be replaced with zeros for taxa not observed in the real transfers
permutationStats[is.na(observed), observed := 0]


# We make one plot per TE class ---------------------------------------------

# we thus add a column for TE classes
permutationStats = mutate(permutationStats, class = ifelse(str_detect(superfamily, "DNA/") | str_detect(superfamily, "RC/"), "DNA", "RNA"))

# we compute stats for randomised and observed HTT numbers of the different taxa we selected and for each TE class
# this means that we sum stats obtained over super families within classes
statsPerClass <- permutationStats[, .(
    simulated = sum(simulated),
    minNumber = sum(minNumber),
    maxNumber = sum(maxNumber),
    qLc = sum(qLc),
    qRc = sum(qRc),
    qL5 = sum(qL5),
    qR5 = sum(qR5),
    observed = sum(observed)
), by = .(class, taxon)]


# we count the number of simulated HTT per taxon, to put the taxa the most involved on top
sumSimulated <- statsPerClass[, sum(simulated), by = .(node = taxon)]
statsPerClass <- statsPerClass[order(sumSimulated[match(taxon, node), -V1], class), ]

# building the connected barplots 

# we attribute taxa-specific colour used for the plots  (the same as on previous figures)
colTaxa <- fread("colorsTaxa.txt")
colTaxa <- left_join(meta[, c("group", "order")], colTaxa)
statsPerClass[, col := colTaxa[match(taxon, order), col]]


########################################

#Final figure showing results of statistical test

pdf(
    file = paste0("Figure3_", included, ".pdf"),
    width = 8,
    height = 4.5
)


# in this layout, the rightmost section is there to leave room for the legend 
layout(cbind(1, 2, 3), widths = c(1, 1, 0.4))

par(
    lwd = 0.5,
    xpd = NA,
    las = 1,
    mai = c(0.6, 0.55, 0.55, 0.6)
)

# the barplots for DNA transposons
linkedPlots(statsPerClass[class == "DNA"], ylab = "number of transfers")

#and for retrotranspons
linkedPlots(statsPerClass[class == "RNA"], legend = T)

dev.off()