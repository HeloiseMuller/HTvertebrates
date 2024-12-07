# ## %######################################################%##
# #                                                          #
# ####  This stage computes Ks between TEs for filtered   ####
# ####   TE-TE hits, to select hits that represent HTT.   ####
# #                                                          #
# ## %######################################################%##
# 
library(data.table)
library(stringi)
library(parallel)
library(dplyr)
library(stringr)
 
path = "~/Project/"
setwd(path) 
 
# # this script uses :

# # - the TE-TE hits we selected in the previous stage
selectedHits <- fread("allOCC2000.txt")

# - the blastx hits between TE copies and proteins also generated in the previous stage
# (imported later)

# the final outputs are :
# - a fasta file of TE copies involved in the hits (as it will be used in further stages)
# - the table of filtered TE-TE hits with new columns for dN and dS values between copies

########################################

# STEP ONE, we extract TE copy sequences --------------------------------------------
# they were actually not extracted in the blastx we preformed prior, as we piped directly from seqtk to diamond

copies <- rbind(select(selectedHits, c(copie1, assembly1)), setnames(select(selectedHits, c(copie2, assembly2)), c("copie1", "assembly1")))
copies <- distinct(copies)

# we will write copy names for each species so we can use seqtk to extract them from fastas -------
# so we split copies per species as 
copiesPerSpecies <- split(copies$copie1, copies$assembly1)

# where the copy sequences will go
dir.create("TEKs/selectedCopies", recursive = T)

# the file of TE copy names that seqtk will use
fileNames <- stri_c("TEKs/selectedCopies/", names(copiesPerSpecies), ".txt")

# we write these to disk
m <- Map(
    f = write,
    x = copiesPerSpecies,
    file = fileNames
)
# ----------


# we import copy sequences -------------
# we list gzip fasta files of TE copies (generated in stage 2)
fasNames <- stri_c("RepeatMasker/copies/", names(copiesPerSpecies), ".TEs.fasta")

# we import selected copies from TE fastas in parallel
seqtk <- function(fas,
                  bed,
                  out,
                  ex = "seqtk subseq",
                  formated = F) {
  # calls seqtk to return sequence from a fasta fas, based on bedfile bed.
  # Frite to fasta file out (if specified) or return to R (if not). 
  # if formatted is TRUE, return sequences as a DNAStringset (if not, returns a character vector 
  # corresponding to the fasta) "ex" specifies how seqtk subseq should be executed
  
  if (missing(out)) {
    seqs <- system(paste(ex, fas, bed), intern = T)
    
    if (formated) {
      f <- stri_sub(seqs, 1, 1) == ">"
      dt <- data.table(content = seqs[!f], id = cumsum(f)[!f])
      concat <- DNAStringSet(dt[, stri_flatten(content), by = id]$V1)
      names(concat) <- stri_sub(seqs[f], 2L, nchar(seqs[f]))
      return(concat)
      
    } else {
      return(seqs)
    }
    
  } else {
    system(paste(ex, fas, bed, ">", out))
  }
}

seqs <- mcMap(
    f = seqtk,
    fas = fasNames,
    bed = fileNames,
    #no out so will not write any output. Fasta returned in seqs
    mc.cores = 20,
    mc.preschedule = F
)

# we stack these sequence in a single DNAStringset
allSequences <- unlist(seqs, use.names = F)

# and write it to a single fasta
write(allSequences, file = "TEKs/selectedCopiesAA300cl2000.fas") 

########################################

# STEP TWO, we prepare and launch computations of dS value of TE hits. -------------------------------------------
# This is based on blastx results of copies against repeat proteins (used by repeat modeler for TE classification)


# for this, we import "raw" results from the blastx performed prior
# which differ from blastx results used to select hits, as those had successive HSPs combined.
blastx <- fread(
    input = "TEs/blastx/allRaw.out", 
    header = F, 
    sep = "\t",
    col.names = c(
        "query", "prot", "pID", "length", 
        "mismatches", "gapopen", "qStart", "qEnd", 
        "sStart", "sEnd", "eValue", "score", "qlen"
        ),
)

#If query is a subpart of the copy, its name is: copyName:start-end --> we want to get copyName only
copyFields <- data.frame(blastx[, stri_split(query, fixed = ":", simplify = T)])

copyFields = data.table(X1 = copyFields$X1, X2 = copyFields$X2, X3 = copyFields$X3)

# we reconstruct copy names as they were in the selectedHits table
#If copy name under form rnd or ltr, copy name = X1 + X2 & coordinates in X3 if subpart
#If copy name without form rnd nor ltr (=else of ifelse), copy name = X1  & coordinates in X2 if subpart
copyFields = mutate(copyFields, copy = ifelse(str_detect(X2, "rnd-") | str_detect(X2, "ltr-"), stri_c(X1,X2, sep=":"), "")) #when else I would prefer to replace by X1, but there is a bug (put a number Instead)
copyFields[which(copyFields$copy==""),]$copy = copyFields[which(copyFields$copy==""),]$X1 #Pour contourner le bug
copyFields = mutate(copyFields, starts = ifelse(str_detect(X2, "rnd-") | str_detect(X2, "ltr-"), as.integer(splitToColumns(X3, "-", columns=1)), as.integer(splitToColumns(X2, "-", columns=1)))) 
copyFields[is.na(copyFields$starts),]$starts <- 1L

# we select hits between copies present in the blastn hits
f <-  copyFields$copy %chin% selectedHits[, c(copie1, copie2)]
blastx <- blastx[f]
copyFields <- copyFields[f,]

# we also replace copy names with shorter names matching the TE-TE hits
blastx[, copy := copyFields$copy]

# so we can convert blastx coordinates in original copy coordinates
blastx[, c("qStart", "qEnd") := .(qStart + copyFields$starts - 1L, qEnd + copyFields$starts - 1L)]

# we write the file for safety, but it is not used afterwards
write.table(blastx, "TEKs/blastxOCC2000.out", col.names = F, na = "NA",  sep = "\t", row.names=F, quote=F)

# we select hits between TEs with sufficient evalue 
blastx <- blastx[eValue <= 0.001]

#we create new start and end coordinates so that the former is always lower
blastx[, c("start", "end") := .(pmin(qStart, qEnd), pmax(qStart, qEnd))]

# we select regions of TE copies covered by just one hsp, hence for which there is no visible conflict between reading frames.
cov <- blastx[, genomeCov(data.table(copy, start, end), successive = F)]
noOverlappedHSP <- cov[cov == 1L]

# we retrieve the hit covering each of these regions (the row index of the blastx table)
# the new column we add is the row index of the hit in the blastx table
#e.g. if hit==13747, it means that it corresponds to blastx[13747,]
noOverlappedHSP[, hit := assignToRegion(blastx[, data.table(copy, start, end)], 
                                        data.table(copy, start))]

# we retrieve the start position of the hit on the protein, and whether the hit was in
# reverse orientation. This is necessary to determine positions in codons later
noOverlappedHSP[, c("protStart", "rev") := blastx[hit, .(qStart, qStart > qEnd)]]

# we remove columns we no longer need
noOverlappedHSP[, c("cov", "hit") := NULL]

write.table(noOverlappedHSP, "TEdS/blastxOCC2000cov1.out", col.names = T, na = "NA",  sep = "\t", row.names=F, quote=F)

# dN dS computations performed via this script :
system("Rscript 07bis-TEdNdS.R allOCC2000.txt TEdS/blastxOCC2000cov1.out TEdS/selectedCopiesAA300cl2000.fas TEdS/ 30")

########################################

# STEP THREE, we remove hits that may result from VT based on dS --------------------------------------------

# we gather results and merge them with the table of selected hits
TEdS <- fread("TEdS/alldNdS.txt")

# we add an integer identifier to the hits, as it was used in TEKaKs.R 
selectedHits[, hit := 1:.N]

#this allows merging the dN dS results with our table of hits
selectedHits <- merge(
    x = selectedHits,
    y = TEdS,
    by = "hit",
    all = T,
    suffixes = c("", ".aa")
)

#May NA because TEdS has less lines than selectedHits.
selectedHits = filter(selectedHits, ! is.na(vks))

# I would like to save selectedHots, but it so big it would not be readable.
# So I firstly remove some columns I don't need anymore
selectedHits = select(selectedHits, -c(superF_blastx.2, nMut, K80distance, rawDistance))
    #Since I already filtered superF_blastx.1==superF_blastx.2, I can keep only one of those colomns
    #I don't need the stats that I remove

# we write the table for safety, this file is not used afterwards
write.table(selectedHits, "allOCC2000dS.txt", quote=F, row.names=F)


# file of filtered BUSCO dS (200 AA alignments and one score per BUSCO gene per pair of clades) generated in stage 05-coreGenedS.R
dS <- fread("Busco/dS/dNdS/dS100AA_filterBusco_median.txt") 

# 0.5% quantiles of dS, per pair of sister clades
dSQuantile <- dS[, .(q05 = quantile(dSmedian, 0.005)), by = mrca]

# gets the relevant quantile for every hit, that corresponding to the clade pair
q05dS <- dSQuantile[match(selectedHits$mrca, mrca), q05]

# we select hits that should constitute HTT, according to our criteria 
# the dS value must be lower than 99.5% of the core gene dS
# this value must be computed on at least 100 codons and be lower than 0.5
httHits <- selectedHits[ks < q05dS & length.aa >= 100L & ks < 0.5, ]

write.table(httHits, "occ2000dS05.txt", quote=F, row.names=F)
