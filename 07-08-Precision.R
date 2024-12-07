#############################################################
#                                                           #
####     This stage prepares is to get a precise pID     ####
#                   and to add some hits                    #
#############################################################


# This script is an addition of the initial pipeline of Zhang et al. 2020
# which is run between the scripts 7 and 8.
# We added this step for two reasons :
#   - In this new pipeline we are using MMseq2
#   instead of blast, which give unprecise pID. 
#   We replaced blast by MMseq2 because of its speed.
#   Now that we identittfied copies resulting of HTT, 
#   out table is much smaller, so blast is runnable.
#   - Increase the number of hits.
#   A small amount of hits resulting of a single HTT event
#   could actually result in several artificial events
#   after clustering.
#   Additional hits will help to form more bridges
#   to cluster them together
#   Note: additional hits are for cluster with less than 200 occurences

# Once the blast is done, the next steps of this scipt is very redundant with what was done in script 06

# Thist script also need the script 07-08bis-changeID.sh

library(data.table)
library(Biostrings)
library(dplyr)
library(stringi)
library(strex)

path = "~/Project/"
setwd(path) 

source("HTvFunctions.R")

# where output file will go
dir.create("TEs/clustering")

# selected TE-TE hits from the previous stage
httHits <- fread("occ2000dS05.txt")


################################################
################## Functions ##################
################################################

extract_fasta <- function(bed, out){
  write.table(select(bed, copy), paste0(out, ".bed"), quote=F, row.names=F, col.names=F)
  system(paste0("seqtk subseq RepeatMasker/copies/", unique(bed[,'assembly']), ".TEs.fasta ", out, ".bed > ", out, ".fasta"))
  system(paste0("rm ", out, ".bed"))
  system(paste0("bash TEs/clustering/changeID.sh /home/muller/Project_Aquatic/TEs/clustering/dcMegablast/IDcopies.txt ", out, ".fasta ", out, ".IDs.fasta"))
  system(paste0("rm ", out, ".fasta"))
}

################################################
################################################
################################################

# STEP ONE, we blast copies involved in HTT 

httHits <- mutate(httHits, pair = paste0(assembly1, "-", assembly2))
length(unique(httHits$pair)) #Gives the number of pairs involved in HTT

summary_httHits <- httHits %>% group_by(pair,superF_blastx.1) %>% summarize(n=n())

dir.create("TEs/clustering/dcMegablast/copies", recursive = T)

#Get all copies (.1 or .2) of the table, with its consensus name, its assembly and its superfamily
copies <- get_copies(httHits)

#we give a distinct number for each copies to reduce the size of blastn files, and improve speed in further stages
copies <- mutate(copies, ID = 1:nrow(copies))

#We save this table
write.table(copies, "TEs/clustering/dcMegablast/IDcopies_full.txt", row.names=F, col.names=T, quote=F, sep='\t')
#We also save a version with only 2 columns. Will be read by 07-08bis.
write.table(select(copies, c(copy, ID)), "TEs/clustering/dcMegablast/IDcopies.txt", row.names=F, col.names=F, quote=F, sep='\t')


#Look at what is the biggest db
sizeDB <- copies %>% group_by(assembly, rep_superF) %>% summarize(n=n())
max(sizeDB$n) 
#IMPORTANT When we will do self blast, max-seqs will have to be at least this size

#We split per assembly
copies <- mutate(copies, sp_family = paste0(assembly, "-", rep_superF))
copies_list <- split(copies, copies$sp_family)

# Extracr fasta of these copies, using the ID numbers
lapply(copies_list, function(x){
    family = sub(".*[/]","", unique(x$rep_superF))
    out = paste0("TEs/clustering/dcMegablast/copies/", unique(x$assembly), "_", family)
    if(!file.exists(paste0(out, ".IDs.fasta"))){ 
      extract_fasta(x, out)
    }
})


files <- list.files("TEs/clustering/dcMegablast/copies", pattern=".IDs.fasta$", full=T) 
f <- !file.exists(paste(files, ".nin", sep = ""))

makedb <- function(fasta) {
    system(paste("makeblastdb -in ",  fasta, " -dbtype nucl "),  ignore.stderr = T)
    return(NULL)
}

# we apply the above function in parallel on 10 CPUs
res <- mcMap( 
    f = makedb,
    fasta = files[f],
    mc.cores = 30,
    mc.preschedule = F
)

pairs <- select(httHits, c(assembly1, assembly2, superF_blastx.1)) %>% distinct() %>% mutate(., family=sub(".*/", "", superF_blastx.1))
files_comp <- sub(".IDs.fasta", "", sub(".*/", "",files)) 
files_tb <- data.frame(file = files, assembly = str_before_nth(files_comp, "_", 2), family = sub(".*_","", files_comp))
pairs <- left_join(pairs, files_tb, by=c("assembly1" = "assembly", "family" = "family")) %>%
    left_join(., files_tb, by=c("assembly2" = "assembly", "family" = "family"), 
    suffix = c(".1",".2")
    ) 

#Make file to run blastn on all the pairs
# Here is an example to run on a server NOT using slurm
dir.create("TEs/clustering/dcMegablast/out")
sink("TEs/clustering/dcMegablast/commandblast")
paste("cd", getwd())
paste("blastn -task dc-megablast -query ",  pairs$file.1, " -db ",  pairs$file.2,  " -max_target_seqs 100 -max_hsps 5 -outfmt 6 -num_threads 20 | awk '{if ($3>=75 && $4>=300 && $12>=200) print $0}' > TEs/clustering/dcMegablast/out/dcMegablast_", pairs$assembly1, "_vs_", pairs$assembly2, "_", pairs$family, " & PID=$!", sep="") #$1=querry, $2=target, $3=threads
paste("blastn -task dc-megablast -query ",  pairs$file.2, " -db ",  pairs$file.1,  " -max_target_seqs  100 -max_hsps 5 -outfmt 6 -num_threads 20 | awk '{if ($3>=75 && $4>=300 && $12>=200) print $0}' > TEs/clustering/dcMegablast/out/dcMegablast_", pairs$assembly2, "_vs_", pairs$assembly1,  "_", pairs$family, " & PID=$!", sep="") #$1=querry, $2=target, $3=threads
sink()


#ATTENTION: do not forget to add the following before running blastn (moslty WAIT !)
#sed -i 's/.*"cd/"cd/' TEs/clustering/dcMegablast/commandblast
#sed -i 's/.*"blastn/"blastn/' TEs/clustering/dcMegablast/commandblast
#sed -i 's/"//g' TEs/clustering/dcMegablast/commandblast
#sed -e 's/$/\nwait $PID/' -i TEs/clustering/dcMegablast/commandblast

# One can split this file with bash to run several blastn in parallel
# In this example, we split the 157,538 lines into 4 files
#sed -n '1,39384p' TEs/clustering/dcMegablast/commandblast >TEs/clustering/dcMegablast/commandblast_1 
#sed -n ' 39385,78769p' TEs/clustering/dcMegablast/commandblast >TEs/clustering/dcMegablast/commandblast_2 
#sed -n '78770,118154p' TEs/clustering/dcMegablast/commandblast > TEs/clustering/dcMegablast/commandblast_3
#sed -n '118155,157538p' TEs/clustering/dcMegablast/commandblast > TEs/clustering/dcMegablast/commandblast_4 

################################################

# STEP TWO we filter hits

dir.create("TEs/clustering/dcMegablast/filtered")
pairs <- mutate(pairs,
    out = paste0("TEs/clustering/dcMegablast/out/dcMegablast_", assembly1, "_vs_", assembly2, "_", family),
    rev = paste0("TEs/clustering/dcMegablast/out/dcMegablast_", assembly2, "_vs_", assembly1, "_", family),
    filtered =  paste0("TEs/clustering/dcMegablast/filtered/dcMegablast_", assembly1, "_vs_", assembly2, "_", family, "_filtered")
)

files <- list.files(path="TEs/clustering/dcMegablast/out", pattern = "dcMegablast_", full.names = T, recursive = T)

#Normaly, this line should not change anything
pairs <- filter(pairs, out %chin% files & rev %chin% files) 

hitNumbers <- pairs[, mcMap(
    f = filterHits,
    out,
    rev,
    filtered,
    mc.cores = 20,
    mc.preschedule = F
)] 
#NOTE: If out & rev are empty, it will not output anything for that pair

hitNumbers <- data.table(do.call(rbind, hitNumbers))
hitNumbers <- mutate(hitNumbers,
    pair = str_before_nth(sub("_vs_", "-", sub("dcMegablast_", "", V1)), "_", 3), 
    superF_blastx.1 = str_after_nth(sub("_filtered", "", V1), "_", 6)
    )
colnames(hitNumbers) = c("file", "pair", "nHit_dcMegablast", "nbHit_filtered")

write.table(summary_httHits, "TEs/clustering/dcMegablast/summary_nbHits.txt", sep='\t', col.names=T, row.names=F, quote=F) 

# on bash:
# cat TEs/clustering/dcMegablast/filtered/*_filtered > TEs/clustering/dcMegablast/filtered/all.score200.dcMegablast.out 
# If there are too many files, one has to split the job.
# Example to split 39176 into two jobs:
    # ls -lh | grep "_filtered" | head -n19588 > cat1-19588
    # ls -lh | grep "_filtered" | tail -n19588 > cat19589-39176

    #Then on each file, do :
    # awk '{print $9}' cat1-19588 > a 
    # mv a cat1-19588
    # cat cat1-19588 |  tr '\n' ' '  > cat1-19588.sh #put file names on the same line
    # sed 's/^/cat /'  cat1-19588.sh -i #add cat at the begining of the line
    # sed 's/$/ >  cat1-19588.quantile005score200/'  cat1-19588.sh -i #add name of the concatenated file at the end of the line
    # bash  cat1-19588.sh
        
    #once done on both files:  cat cat*200 > all.score200.dcMegablast.out


################################################

#STEP THREE: we keep only hits with at least 300bp of coding (already the same superF since I did blast between same superF)

TEhits <- fread("TEs/clustering/dcMegablast/filtered/all.score200.dcMegablast.out",  
    col.names = c("ID.1", "ID.2", "pID", "length", "qStart", "qEnd", "sStart", "sEnd", "score")
    )
copies <- fread("TEs/clustering/dcMegablast/IDcopies_full.txt")

#add name of copy and info about it (assembly, consensus, superF)
TEhits <- left_join(TEhits, copies, by=c("ID.1" = "ID")) %>% 
    left_join(., copies, by=c("ID.2" = "ID"), suffix = c(".1",".2")) 

## From here, the next steps are just a copy and paste of what was run in script 06
# Many hits will obsviouly be considered as resulting of HT
# Yet among the additional hits, there will be some vertical tansfers.
# This is why we have to go through all these steps again

# we import the combined blastx results generated in 06tris 
# (no need to generate that file again)
blastx <- fread(
    input = "TEs/blastx/all.copies.successiveBlastx.out",
    header = T,
    sep = "\t",
    col.names = c(
        "query", "subject", "pID", "length", "mismatch", "gapopen",
        "qStart", "qEnd", "sStart", "sEnd", "evalue", "score", "ex"
    )
)

#Combine coordinates by copy
#Here we don't really care if the copy has hits on different proteins. 
#This step is just to select coordinates of interest for the next step (better to select too many than not enough)
perCopy <- blastx[ ,
    .( start = min(c(qStart, qEnd)),
        end = max(c(qStart, qEnd))
        ),
    by = query
]

# to select TE-TE hits that cover a protein region that is long enough,
# we need to get the blastn coordinates where starts always < ends
TEhits[, c("qSt", "qEn", "sSt", "sEn") := data.table(
    pmin(qStart, qEnd),
    pmax(qStart, qEnd),
    pmin(sStart, sEnd),
    pmax(sStart, sEnd)
    )]

# we add protein ranges for query and subject in blastn hits
TEhits[, c("qS", "qE") := perCopy[chmatch(copy.1, query), .(start, end)]]
TEhits[, c("sS", "sE") := perCopy[chmatch(copy.2, query), .(start, end)]]

# we prevent bugs that are caused by NAs in intersection(). Happens when the copy does not hit a prot
TEhits[is.na(qS), c("qS", "qE") := 0L]
TEhits[is.na(sS), c("sS", "sE") := 0L]

# we get the region of the query that aligns on the subjects and also aligns on proteins
TEhits[, c("qS", "qE") := data.table(intersection(qSt, qEn, qS, qE, T))]
TEhits[, c("sS", "sE") := data.table(intersection(sSt, sEn, sS, sE, T))] # same, for the subject

# converts the above into query coordinates
TEhits[, c("qsS", "qsE") := .(qSt + sS - sSt, qEn + sE - sEn)]

# so we can compute the lenght of the TE alignment that also aligns on proteins
TEhits[, inter := pmin(qE, qsE) - pmax(qS, qsS) + 1L]
TEhits[is.na(inter) | inter < 0L, inter := 0L]

# we remove columns we no longer need
TEhits[, c("qSt", "qEn", "sSt", "sEn", "qS", "qE", "sS", "sE", "qsS", "qsE") := NULL]

# we select hits  with sufficient protein region (300 bp)
TEhitsCoveringProteins <- TEhits[inter >= 300L] 

#Save
write.table(TEhitsCoveringProteins, "TEs/clustering/dcMegablast/all300aa_dcMegablast.txt", row.names = F, na = "NA",  sep = "\t") 

################################################

#STEP FOUR: single linkage clustering: don't keep more than 2000 copies per cluster and pair of species

# we place the best hits on top (longest protein region then higher score)
setorder(TEhitsCoveringProteins, -inter, -score)

# we add an integer column denoting the single-linkage cluster of copies, based on hits
TEhitsCoveringProteins[, cl := clusterFromPairs(copy.1, copy.2)]

# we add an integer column indicating number of times we have seen a cluster for each a species pair.
# Within a species pair, hits of the same cluster would represent the same transfer and are thus redundant.
TEhitsCoveringProteins[, occ := occurrences(data.table(cl, assembly.1, assembly.2))]

summary = mutate(TEhitsCoveringProteins, pairs = paste0(assembly.1, "-", assembly.2)) %>% group_by(cl, pairs) %>% summarize(nbHit=n())


# we select at most 2000 hits between copies from a given cluster for any species pair
sel <- TEhitsCoveringProteins[occ <= 2000L, ]

# we also add an integer column for the mrca of each species pair, for later use
mrca <- fread("Busco/pairsALL.txt")
#assemblies don't have version so I have to remove them also in sel:
sel$assembly.1bis <- sub("[.].*","", sel$assembly.1)
sel$assembly.2bis <- sub("[.].*","", sel$assembly.2)

# We add the column mrca, but assembly.1 - assembly.2 are not necessary in the same order as in sel
# This is why we have those additional lines
sel <- left_join(sel, mrca, by = c("assembly.1bis" = "assembly.1", "assembly.2bis" = "assembly.2")) %>% 
    left_join(., mrca, by = c("assembly.1bis" = "assembly.2", "assembly.2bis" = "assembly.1")) %>% 
    mutate(., mrca = ifelse(is.na(mrca.x), mrca.y, mrca.x), 
            divTime = ifelse(is.na(divTime.x), divTime.y, divTime.x),
            species.1 = ifelse(is.na(sp1.x), sp2.y, sp1.x), species.2 = ifelse(is.na(sp1.x), sp1.y, sp2.x)
           ) %>%
    # Remove columns we don't need anymore
    select(., -c(mrca.y, mrca.x, divTime.y, divTime.x, sp1.x, sp1.y, sp2.x, sp2.y, assembly.1bis, assembly.2bis)
    )

# One may check that we wrote the right copy with the write assembly with :
sel[assembly.1 != gsub("-.*", "", copie1),] #should be empty
    
# We need to rename column of copies because of script 07bis
sel = rename(sel, copie1 = "copy.1", copie2 = "copy.2") 

# and write the filtered hits to disk
fwrite(sel, "TEs/clustering/dcMegablast/allOCC2000_dcMegablast.txt", col.names = T, row.names = F, na = "NA",  sep = "\t", quote=F) 

################################################

# STEP FIVE, we extract TE copy sequences --------------------------------------------
# they were actually not extracted in the blastx we preformed prior, as we piped directly from seqtk to diamond

copies <- rbind(select(sel, c(copie1, assembly.1, ID.1)), setnames(select(sel, c(copie2, assembly.2, ID.2)), c("copie1", "assembly.1", "ID.1")))
copies <- distinct(copies)

# we will write copy names for each species so we can use seqtk to extract them from fastas -------
# so we split copies per species as 
copiesPerSpecies <- split(copies$copie1, copies$assembly.1)

# where the copy sequences will go
dir.create("TEs/clustering/dcMegablast/TEdS/selectedCopies", recursive = T)

# the file of TE copy names that seqtk will use
fileNames <- stri_c("TEs/clustering/dcMegablast/TEdS/selectedCopies/", names(copiesPerSpecies), ".txt")

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
write(allSequences, file = "TEs/clustering/dcMegablast/TEdS/selectedCopiesAA300cl2000_dcMegablast.fas")

################################################

# STEP SIX we compute dS

system("Rscript 07bis-TEKaKs.R TEs/clustering/dcMegablast/allOCC2000_dcMegablast.txt TEdS/blastxOCC2000EXcov1.out TEs/clustering/dcMegablast/TEdS/selectedCopiesAA300cl2000_dcMegablast.fas TEs/clustering/dcMegablast/TEdS/ 30")
    #$1 are the hits that passed the filtered after blastn max_target_seqs 100
    #$2 are info about coding region (not new since 07)
    #$3 is the fasta of all TE copies in $1
    #$4 is where to save output
    
#change name output
system("mv TEs/clustering/dcMegablast/TEdS/alldNdS.txt TEs/clustering/dcMegablast/TEdS/alldNdS_dcMegablast.txt")
    
################################################

# STEP SEVEN, we remove hits that may result from VT based on dS

# we gather results and merge them with the table of selected hits
TEdS <- fread("TEs/clustering/dcMegablast/TEdS/alldNdS_dcMegablast.txt") 
selectedHits <- fread("TEs/clustering/dcMegablast/allOCC2000_dcMegablast.txt")

# we add an integer identifier to the hits, as it was used in TEKaKs.R 
selectedHits[, hit := 1:.N]

#this allows merging the Ka Ks results with our table of hits
selectedHits <- merge(
    x = selectedHits,
    y = TEdS,
    by = "hit",
    all = T,
    suffixes = c("", ".aa")
)

selectedHits = filter(selectedHits, ! is.na(vks))

# we write the table for safety, this file is not used afterwards
write.table(selectedHits, "TEs/clustering/dcMegablast/allOCC2000dS_dcMegablast.txt", quote=F, row.names=F)


# file of filtered BUSCO dS generated in stage 05-coreGenedS.R
dS <- fread("Busco/dS/dNdS/dS100AA_filterBusco_median.txt") 

# 0.5% quantiles of dS, per pair of sister clades
dSQuantile <- dS[, .(q05 = quantile(Ksmedian, 0.005)), by = mrca]

# gets the relevant quantile for every hit, that corresponding to the clade pair
q05dS <- KsQuantile[match(selectedHits$mrca, mrca), q05] 

# we select hits that should constitute HTT, according to our criteria 
# the dS value must be lower than 99.5% of the core gene dS
# this value must be computed on at least 100 codons and be lower than 0.5
httHits <- selectedHits[ks < q05dS & length.aa >= 100L & ks < 0.5, ]

write.table(httHits, "TEs/clustering/dcMegablast/occ2000dS05_dcMegablast.txt", quote=F, row.names=F) 

################################################

# STEP height, we remove dS<0 & do single linkage clustering again

#dS<0 means couldn't compute dS
httHits <- httHits[ks>=0,] 

# per single-linkage cluster, as hits from the same cluster are very likely
# to represent the same HTT
# For more details see Zhang et al. 2020: 
    # We applied single-linkage clustering to connect any two hits sharing a TE copy (e.g hit A vs B and hit A vs C share the copy A)
    # and we retained no more than 200 hits per resulting hit cluster per species pair. 

# We have to calculate again position in clustering, because it is not impossible that
# occ 190 did not pass the filter but 201 did. In this case, we do want to keep occ 201 

# We favor hits that involve the longest protein coding regions
setorder(httHits, -inter, -score)
httHits[, occ2 := occurrences(data.table(cl, assembly.1, assembly.2))]
httHits <- httHits[occ2 <= 200L, ]
    
# and write them to file
write.table(httHits, "TEs/clustering/dcMegablast/occ200ds05_dcMegablast.txt", quote=F, row.names=F)