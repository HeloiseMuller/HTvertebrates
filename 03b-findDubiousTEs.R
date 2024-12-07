library(data.table)
library(stringi)


path = "~/Project/"
setwd(paste0(path), "DB_clustered/findDubious/")

source("../../HTvFunctions.R")

######################
##### Process nr #####
######################

blast_nr <- fread(input = "TE_cat_all.long.clustered_rep_seq_on_nr.out", header = F,  sep = "\t",
    col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) 

# remove hits of low evalue
blast_nr <- blast_nr[evalue <= 0.001]
setorder(blast_nr, qseqid, qstart, -qend)

# remove hits nested within others (in the query), as these hits are of lower quality (shorter)
blast_nr[, set := regionSets(data.frame(qseqid, qstart, qend))]
blast_nr <- blast_nr[!duplicated(set)]

# we combine relatively adjacent hits (with a lot of margin). This function also renames columns
combined <- combineHits(
    blast = blast_nr,
    maxDist = 300,
    maxOverlap = 10000,
    maxDiff = 300,
    blastX = T
)


# puts best hits on top, per query
setorder(combined, query, -score)

# attributes a number to every hit related to the same query (starting at 1)
combined[, hit := occurrences(query)]

#############################
##### Process RepeatPep #####
#############################

blast_repBase <- fread(input = "TE_cat_all.long.clustered_rep_seq_on_RepeatPeps.out", header = F, sep = "\t",
    col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "a")
)

blast_repBase <- blast_repBase[evalue <= 0.001]
setorder(blast_repBase, qseqid, qstart, -qend)
blast_repBase[, set := regionSets(data.frame(qseqid, qstart, qend))]
blast_repBase <- blast_repBase[!duplicated(set)]

# replaces subject name by this generic term
blast_repBase[, sseqid := "repProtein"]

# here we combine all hsps of each consensus, regardless of the protein 
combinedRep <- combineHits(
    blast = blast_repBase,
    maxDist = 100000,
    maxOverlap = 100000,
    maxDiff = 100000,
    blastX = T
)

##################################################
# compare hits on repbase proteins to hits on nr #
##################################################


# the function below confronts the hit on repbase to the hit on nr for each consensus:
dubious <- function(occ) {
    merged <- merge(
        x = combined[hit == occ],
        y = combinedRep,
        by = "query",
        suffixes = c("", ".r"),
        all = T
    )

    # we compute the length of the part that aligns on the nr protein before it aligns on the rep protein
    merged[, leading := pmin(qStart.r, qEnd.r) - pmin(qStart, qEnd)]

    # the length pf the part that aligns on the nr protein after it aligns on the rep protein
    merged[, trailing := pmax(qStart, qEnd) - pmax(qStart.r, qEnd.r)]

    # we retain pairs of hits where there is no alignment on repbase,
    m <- merged[is.na(subject.r) |
        
        # or where there are parts aligning on nr protein not covered by alignment on rep proteins (of at least 90 bp)
        leading > 90 | trailing > 90]
    m
}

# we apply the above fonctions to hits involving every consensus
dubiousConsensusHits <- lapply(unique(combined$hit), dubious)
dubiousConsensusHits <- rbindlist(dubiousConsensusHits)

# we only retain hits on nr proteins that are not annotated as repeat proteins :
dubiousConsensusHits <- dubiousConsensusHits[!grepl(
    pattern = "transpo|retro|reverse|integrase|gag|pol-|pol |rna-dependent|polyprot|mobile|jockey|jerky|setmar|copia|recombinase|crypton|mariner|tcmar|tc1|gypsy|helitron|harbi|piggy|maveri|polinton|academ|ltr|cmc|envelop",
    # the keyword list above has been determined after carefull inspection of protein names in the hits
    x = subject,
    ignore.case = T
)
                                            # and with alignment length >100 
                                            & abs(qStart - qEnd + 1) > 100, ]

# we put accession numbers of nr proteins in this new column
dubiousConsensusHits[, acc := splitToColumns(subject, " ", columns=1)]



##################################################
#  blastp nr proteins similar to TE consensuses #
############ against repeat proteins ############
##################################################

# we put their names in a file for seqtk 
# (we use all proteins, not just those in the dubiousConsensusHits table. We can afford it)
writeLines(unique(blast_nr$sseqid), "hitted_nrProteins.bed")

# we extract the protein sequences on bash:
# seqtk subseq nr.gz hitted_nrProteins.bed > hitted_nrProteins.fasta

# we launch the search on bash:
# nohup diamond blastp -p 20 --more-sensitive -q hitted_nrProteins.fasta  --more-sensitive -d repeatPeps.dmnd --quiet -k 20 -f 6 "qseqid" "sseqid" "pident" "length" "mismatch" "gapopen" "qstart" "qend" "sstart" "send" "evalue" "bitscore" -o nrProtsOnRepeatPeps.out &> nohup.nrProtsOnRepeatPeps &

# processing the hits --------------------------------------
blastp <- fread("nrProtsOnRepeatPeps.out",
    header = F,
    sep = "\t",
    col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
)

blastp <- blastp[evalue < 0.001]

# again, combining HSPs on the same query
combinedP <- combineHits(
    blast = blastp,
    maxDist = 1000,
    maxOverlap = 100000,
    maxDiff = 1000
)
    

# we get accession number of nr proteins that have an acceptable hit with a rebpase proteins.
# We will consider them as legit TE proteins
repProtAccessions <- combinedP[abs(qEnd - qStart) + 1 >= 100 & identity >= 35, 
                               unique(splitToColumns(query, " ", columns=1))]


# so we can extract consensuses that hit to other proteins than these
dubiousFamilies <- dubiousConsensusHits[!acc %in% repProtAccessions, unique(query)]

# we write these family names to disk. 
write(dubiousFamilies, "familiesToIgnore.txt")

# We still blast these TEs to find HTT.
# This leaves us the possibility to remove these TEs afterwards
