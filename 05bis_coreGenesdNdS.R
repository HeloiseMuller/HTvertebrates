## %######################################################%##
#                                                           #
####        This script computes pair-wise dN and        ####
####         dS between species at BUSCO genes           ####
#                                                           #
## %######################################################%##

library(seqinr)

# this script is run at stage 05-coreGenedS.R

# the only user input is the number of CPUs to use
args <- commandArgs(trailingOnly = TRUE)
nCPUs <- as.integer(args[1])

dir.create("dS/dNdS")

# we import all protein and CDS sequences at once
allAA <- readAAStringSet(list.files("CDS", pattern = ".translated.fas", full.names = T))
allCDS <- readDNAStringSet(list.files("CDS", pattern = ".CDS.fas", full.names = T))

# we list results files of mmseq2 easy-rbh searches (tabular)
mmseqRBH_Files <- list.files(
    "dS/mmseqRBH",
    pattern = ".out",
    recursive = T,
    full.names = T
) 

#We keep only the pairs for which we have an mmseq output
pairs = fread("pairsToMMseqRBH.txt")
pairs = filter(pairs, out %in% mmseqRBH_Files)

pairs = split(pairs, pairs$out)

# our function to compute dN and dS of pairs of genes based on reciprocal blastp results (one species vs another)
compute_dS <- function(pairs) {
    
   
    # the fields present in the blast outputs
    fields <- c("query", "subject", "identity", "length", "mismatches", "indels",
                "qStart", "qEnd",  "sStart","sEnd",  "eValue", "score")

    RBH <- suppressWarnings(fread(pairs$out, header = F, col.names = fields, sep = "\t"))
    
    # retain only the best alignment (HSP) per hit (because there can be overlapping alignments)
    setorder(RBH, -score)
    RBH <- RBH[!duplicated(data.table(query, subject)), ]
    
    # because there appears to be a bug with diamond where rare hits have coordinates longer than the query or subject
    RBH <- RBH[qEnd <= nchar(allAA[query]) & sEnd <= nchar(allAA[subject])] 
    
    #returns empty results if there are no hits (almost impossible)
    if (nrow(RBH) == 0) return (data.table()) 
    
    # we extract protein regions involved in hits. 
    # If the alignment doesn't work (it appears that Biostring sometimes has issues with aligning AAStringSets), 
    # it could be useful to replace subseq() below by stri_sub(), which returns character vectors
    #Mon com: fonctionne pas avec RBH post diamond mais fonctionne avec mmseq easy-rbh (enfin avec stri_sub)
    sp1Seqs <- RBH[, stri_sub(allAA[query], qStart, qEnd)] #Jean: sp1Seqs <- RBH[, subseq(allAA[query], qStart, qEnd)]
    sp2Seqs <- RBH[, stri_sub(allAA[subject], sStart, sEnd)] #Jean: sp2Seqs <- RBH[, subseq(allAA[subject], sStart, sEnd)]

  
    # aligns these regions
    aln <- alignWithEndGaps(sp1Seqs, sp2Seqs)
    
    # extracts CDS regions corresponding to aligned protein regions
    sp1Nuc <- RBH[, stri_sub(allCDS[query], qStart * 3 - 2, qEnd * 3)]  #Jean: sp1Nuc <- RBH[, subseq(allCDS[query], qStart * 3 - 2, qEnd * 3)] 
    sp2Nuc <- RBH[, stri_sub(allCDS[subject], sStart * 3 - 2, sEnd * 3)] #Jean: sp2Nuc <- RBH[, subseq(allCDS[subject], sStart * 3 - 2, sEnd * 3)]
    
    # converts protein alignments into nucleotide alignments
    sp1Nuc <- aaToCDS(aln$pattern, sp1Nuc)
    sp2Nuc <- aaToCDS(aln$subject, sp2Nuc)
    
    # makes a list of seqinrAlignment objects for Ka/Ks computation
    seqinrAlns <- apply(cbind(sp1Nuc, sp2Nuc), 1, seqinrAlignment)
    dNdS <- lapply(seqinrAlns, kaks) #kaks is a function from seqinr
    
    # we stack the results into a data.table
    dNdS <- rbindlist(lapply(
        X = dNdS,
        FUN = as.data.table
    ))
    
   res <- data.table(RBH[, c(1, 2, 7:10), with = F], 
                      dNdS[,.(dS = ks, dN = ka)],
                      alnLength = nchar(gsub("-", "", aln$pattern)))
    
   res
}



# we process batches in parallel
res <- mclapply(pairs,
    FUN = compute_dS,
    mc.cores = nCPUs,
    mc.preschedule = F
)

write.table(rbindlist(res), "dS/dNdS/all.dNdS.txt", quote=F, row.names=F, na="NA", sep='\t')