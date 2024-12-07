## %######################################################%##
#                                                          #
####    This stage establishes dN dS rate for pairs     ####
####                 of TEs copies that                 ####
####      diverged by transposition within genomes      ####
#                                                          #
## %######################################################%##


# the principle is to compute dN and dS values for related copies within a genomes
# These copies belong to the same "community" of hits, so are unlikely to represent
# different htt events. They would have diverged by transposition within the genome

# This script uses 
# - the table of HTT hits between TEs, to obtain the communities to which these hits belong (clustring per clade)
# - the self-blastn output of TE copies against themselves, generated at stage 8


library(stringi)
library(data.table)
library(stringr)
library(parallel)
library(ape)

path = "~/Project/"
setwd(path) 

source("HTvFunctions.R")

####################

# STEP ONE, we select relevant blast hits of TEs to compute dN dS rates.--------------------------------

# We do it per superfamily (for big super families, several blast searches were performed in
# parallel so there are multiple output files)

# we retrieve blastn hits of copies against themselves (generated by
# TEselfBlastn.R at stage 08-prepareClustering.R). 
# As these files are very big, we will select only the hits we need. 
# We filter the hits separately from dN dS computation because it takes some time
# Should be need to redo the dN dS computations, we won't have to filter hits again


blastFiles <- list.files(
    path = "TEs/clustering/selfBlastn_perClade/out", 
    pattern = ".out",
    recursive = TRUE,
    full.names = T
)

# some output file may be empty (no hits), which could cause an issue during import
blastFiles <- blastFiles[file.size(blastFiles) > 0]

# we extract superfamily names from file names
superF <- splitToColumns(basename(blastFiles), "_", columns=1)

# we split files names per superfamily
blastFiles <- split(blastFiles, superF)

# we import the copy integer identifiers as these were 
# used in the blasts, instead of very long copy names
ids <- fread("TEs/clustering/dcMegablast/IDcopies.txt", col.names=c("copy", "id"))

# we retrieve the host species of copies, as we will only compare copies within the same genome
ids[, assembly:= splitToColumns(ids$copy, split="-", columns=1)]

meta = fread("metadata.tbl")
meta$species  = gsub(" ","_", meta$species)

ids[, species:= meta[chmatch(ids$assembly, meta$assembly), species]]


# We also import the tree to give integer ids (tip numbers) to species, for speed
tree <- read.tree("datasetTree.nwk")

# we make the correspondance between the copy id (index of
# the vector below, as the id is an integer) and species id
spForCopy <- ids[match(1:max(id), id), chmatch(species, tree$tip.label)]

# we import the HTT hits (including hitgroups that were not retained by
# our criteria in stage 11, just in case we change our filters afterwards)
httHits <- fread("HTThitsAssessed_perClade.txt")

# we extract the information we need. 
# We need to know the "communities", which are clusters of related copies
# that should have diverged only by transposition events, avoiding potential HTT

copiesInCom <- httHits[, union(ID.1, ID.2), by = com]

# we also make a vector for the correspondance between a copy (its index) and a community (its value)
#To note, this will only save the first com of each hits. So we will miss things. But that's fine, we have have enough hits intra
comForCopy <- copiesInCom[match(1:max(ids$id), V1), com] 

dir.create("TEs/TEevolution/selectedHits/", recursive = T) # where selected hits will go


# we will now select the hits we need process self blastn output --------------------------------------------------------------
# the hits must involve copies from the same community and species (genome)
# we do it per super family as there are no hit between super families (no blast)

#I want the outputs of TEs/clustering/selfBlastn_perClade/out_catClades, that were generated by the script iterativeSecondClustering.R
    #However, I didn't save some columns I need here. So I have to generate these output again, with these columns
    #For that, I copy and past the function used in the script iterativeSecondClustering.R, and slightly modify it to keep the columns I need, and change the location of the ouput

getHits <- function(superFi) {
    files <- blastFiles[[superFi]]
    
    import <- function(file) {
        hits <- fread(
            input = file,
            sep = "\t",
            header = F,

            # the function will be used in parallel, so we don't need more than 1 CPU here
            nThread = 1,
            
            #we only need these columns for the clustering (we don't use the score, but we import it just in case)
            drop = c(3, 9), 
            col.names = c("q", "s", "length", "qStart", "qEnd", "sStart", "sEnd")
        )
        cat(".")
        return(hits)
    }

      # we apply the function for each blast file of the super family, in parallel with 20 CPUs
    hits <- mclapply(
        X = files,
        FUN = import,
        mc.cores = min(20, length(files)),
        mc.preschedule = F
    )
    hits = rbindlist(hits)

    # some q-s may have been selected several times, if their species are included in several mrca
    #--> we keep only one q-s hits, they should all be identical anyway
    # N.B. We always have query < subjet
    hits <- hits[!duplicated(data.table(q, s)), ] 

    #here, we don't need hits <300 (since we will have to compute dN & dS)
    hits <- hits[length >= 300L]

    #here, we only want hits between the same species
    hits <- hits[spForCopy[q] == spForCopy[s] & 
                     # and same community
                     #N.B. Here we look only at one one com of each copy, so we will miss many hits
                        #That's fine, we don't need all of them. It will be enough that this.
                        #The improtant is to be sure that those hits diverged by vertical transfer
                     comForCopy[q] == comForCopy[s], ] 

    cat("-")
    return(hits[,-"length"])
}


intraGenomeHits <- rbindlist(lapply(unique(superF), getHits))


fwrite(intraGenomeHits, "TEs/TEevolution/selectedHits/all.out", sep='\t')


# STEP TWO, we filter TE hits that cover a protein region of at least 300 bp-----------------------------
# this is to avoid unnecessary dN/dS computations on codon alignments that will be too short

# we import the a blastx file generated earlier (stage 6), listing
# alignment coordinates of copies against repeat proteins
# this one doe not have overlapping HSPs
blastx <- fread("TEdS/blastxOCC2000cov1.out")

#Normal that many NA since many hits did not reach this step

# we replace copy names by their integer ids, since ids are used in the self blastn files
blastx[, id := ids[chmatch(blastx$copy, copy), id]]

# we only retain blastx results for copies involved in our selected hits
# this will speed things up a bit

blastx <- blastx[id %in% intraGenomeHits[, c(q, s)]]

# as there are several HSPs (alignments) per TE, we get the minimal start
# and maximal end of these, to approximate the protein region covered by a TE
blastxCoords <- blastx[, .(start = min(start), end = max(end)), by = id]

# we add columns of TE hit coordinates where starts always < ends
intraGenomeHits[, c("qSt", "qEn", "sSt", "sEn") := data.table(
  pmin(qStart, qEnd),
  pmax(qStart, qEnd),
  pmin(sStart, sEnd),
  pmax(sStart, sEnd)
)]

# we add protein ranges for query and subject in blastn hits
intraGenomeHits[, c("qSprot", "qEprot") := blastxCoords[match(q, id), .(start, end)]]
intraGenomeHits[, c("sSprot", "sEprot") := blastxCoords[match(s, id), .(start, end)]]

# we obtain the region of the TE query that aligns on the TE subject and also aligns on proteins
intraGenomeHits[, c("qSprot", "qEprot") := data.table(intersection(qSt, qEn, qSprot, qEprot, T))]

# same for the subject TE
intraGenomeHits[, c("sSprot", "sEprot") := data.table(intersection(sSt, sEn, sSprot, sEprot, T))]

# we convert the above into query coordinates
intraGenomeHits[, c("qsS", "qsE") := .(qSt + sSprot - sSt, qEn + sEprot - sEn)]

# we estimate the length of the protein part of copies that align with each other
intraGenomeHits[, inter := pmin(qEprot, qsE) - pmax(qSprot, qsS) + 1L]

# selects hits for which this region is at least 300 bp, and the columns we need
# we also retrieve full TE names as these will be needed for dN dS computation, 
# because these names are those in the fasta file of copy sequences
intraGenomeHits <- intraGenomeHits[inter >= 300L, data.table(q, s,
    query = ids[match(q, id), copy],
    subject = ids[match(s, id), copy],
    qStart, qEnd, sStart, sEnd
)]

selectedHits <- data.table(intraGenomeHits, htt = F)

# we prepare the files needed for the TEdNdS.R script
TEhitFile <- "TEs/TEevolution/selectedHits/selected.out"
blastxFile <- "TEs/TEevolution/selectedHits/blastx.out"

#I did some modification in the script TEdNdS.R, so I adapt my table
colnames(selectedHits)[3:4] = c("copie1", "copie2")


# for these tables, we remove columns that are not used in TEdNdS.R
fwrite(selectedHits[, -c("q", "s", "htt")], TEhitFile,  sep='\t')
fwrite(blastx[, -"id"], blastxFile, sep='\t')


# STEP THREE, we compute dN and dS rates -------------------------------------------------------------

outFolder <- "TEs/TEevolution/dS/"
dir.create(outFolder)
nCPUs <- 30

system(paste(
    "Rscript 07bis-TEdNdS.R",
    TEhitFile, blastxFile,
    "TEdS/selectedCopiesAA300cl2000.fas",  #way more copies than needed, but at least I am sure it contains all I need
    outFolder, nCPUs
))


# We imports results generated by TEdNdS.R (dN dS from VT)
dNdS <- fread(stri_c(outFolder, "/alldNdS.txt"))

selectedHits[, hit := 1:.N] # required for the merged below, as the dNdS results have a hit identifier
dNdS <- merge(selectedHits, dNdS, by = "hit")

#We also want dN dS from HT
httHits <- httHits[retained==T,]

#Need to add distance (this info was removed at a previous stage).
    #For this, read the file genrated at stage 07-08-Precision
httHits_q05 <- fread("TEs/clustering/dcMegablast/occ2000dS05_dcMegablast.txt")

#Keep only hits we kept as HT until the end
httHits_q05[, ID:=ifelse(copie1<copie2, paste0(copie1, "-", copie2), paste0(copie2, "-", copie1))]
httHits[, ID:=ifelse(copie1<copie2, paste0(copie1, "-", copie2), paste0(copie2, "-", copie1))]
httHits_q05 <- httHits_q05[ID %chin% httHits$ID,]

#Merge dNdS of HT & VT
dNdS <- rbind(
  data.table(dNdS[,-"hit"], htt = F), # we add a logical column to differentiate the two types of hits
  data.table(httHits_q05[, colnames(dNdS)[-1], with=F], htt = T) #with=F so variable as column names work
)

#for some reason I don"t have q in my version, so I firstly have to add this
ids <- fread("TEs/clustering/dcMegablast/IDcopies_full.txt")
dNdS[, "q" := .( ids[match(copie1, copy), ID]   )]

#I can now do like Jean:
# we add columns for superfamily names and community ids based on copy ids
dNdS[, c("superF", "com") := .(
    ids[match(q, ID), rep_superF],
    comForCopy[q]
)]

fwrite(dNdS, "TEs/TEevolution/TEdNdSndDistance.txt", sep='\t') 