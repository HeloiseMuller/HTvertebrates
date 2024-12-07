# %######################################################%##
#                                                          #
####           This stage computes dS between           ####
####        species for core genes extracted by         ####
####             BUSCO, in order to compare             ####
####           dS between TEs and core genes            ####
#                                                          #
## %######################################################%##

# the principle of the approach is described in Peccoud et al. 2017 PNAS


# this script uses the fasta files of complete BUSCO genes in each species in a "Busco" folder 
# (The BUSCO pipeline was not automated in scripts)

# example of BUSCO run (launched from bash with version >= 5.4.0:
# busco -i genomes/GCA_xxx.fa -o Busco/GCA_xxx_busco -l database/busco/yyy/ -m genome -c 1
# where yyy is the busco database associated with the genome of GCA_xxx. It should be the database which is as precise as possible (eg aves instead of metazoa)
# ATTENTION after busco was ran one 1 genome, make sure the outputs of interest are there:
  # Busco/GCA_xxx_busco/run_yyy/busco_sequences/single_copy_busco_sequence/zzz.fna
  # Busco/GCA_xxx_busco/run_yyy/busco_sequences/single_copy_busco_sequence/zzz.faa
  # where zzz is the name of the busco gene
# In my case there was no .fna because of a bug that has been fixed in v5.4.0 

# The steps of this script are
  # Step 1: retrieve BUSCO sequences and do similirity seach
  # Step 2: compute dN & dS

# the final outputs of this script 5 are:
# - dS/dNdS/dS100AA_filterBusco_median.txt :
  # a tabular file of dN and dS median values for pairs of BUSCO genes
# - dSdistribution_of_BUSCO_overTime.pdf: figure of the correlation between dS and divergence time of species 

# An important intermediate file that will be used a lot latter in the piepleine is pairsALL.txt
# This table contains all the possible pairs of species (for which we looked for HTT or not) including species against itself
# Its nrow = ((nSpecies * nSpecies - nSpecies) / 2) + nSpecies, where nSpecies is the number of species in the dataset


library(Biostrings)
library(stringi)
library(dplyr)
library(data.table)
library(parallel)
library(ape)
library(ggplot2)

path = "~/Project/"
setwd(paste0(path, "/Busco"))

source("../HTvFunctions.R")


# STEP ONE, we perform diamond blastp search of BUSCO proteins to find orthologs between which we compute dS, and to align proteins afterwards ---------------------------------------
# for this, we do reciprocal blastp search for each pair of species, only retaining the best hit per query

# we will also use the blastp hps coordinates to align proteins 
# at stage 2 since dN dS computations are are based on codons


# we first retrieve the BUSCO CDS sequences --------------------------------------------

folders <- list.dirs(path=".")
folders <- folders[grep("single_copy_busco_sequences", folders, fixed = T)]

getBUSCOs <- function(folder, type) {
    # concatenates busco genes of a species (having a specific folder) in 
    # single fasta and renames sequences by adding the species name to each
    
    #Get the name of the assembly
    assembly <- sub('_busco.*', "", folder) %>% sub("./", "", .)
    
    if(type=="cds"){
        files <- list.files(folder, pattern = ".fna", full.names = T)
    } else if(type=="aa"){
        files <- list.files(folder, pattern = ".faa", full.names = T)
    }

    # readDNAStringSet takes a vector of file paths, but lapply was 
    # necessary to avoid a bug due to too many open connections
    seqs <- lapply(files, readDNAStringSet)

    # concatenating sequences (unlist doesn't reliably work on DNAStringSets)
    concat <- do.call("c", lapply(seqs, as.character))

    # because somehow the c() function above would add the object number in the list to all sequence names
    nams <- unlist(lapply(seqs, names))

    # renaming sequences
    names(concat) <- stri_c(nams, ":", assembly) #Code Jean supprime coor donc on se retrouve avec plein de seq avec le meme nom: names(concat) <- stri_c(stri_extract_first(nams, regex = "[^:]+:"), assembly) 
    
    # and writing them to disk
    if(type=="cds"){
        writeXStringSet(DNAStringSet(concat), stri_c("CDS/", assembly, ".CDS.fas"))
    } else if(type=="aa"){
        writeXStringSet(DNAStringSet(concat), stri_c("CDS/", assembly, ".aa.fas"))
    }
    print(paste("done", assembly))
    
    if(type=='aa'){
        return(data.table(assembly=assembly, nbBusco = length(nams)))
    }
}

translateCDS <- function(cds, aa) {
    # translates CDS into proteins, cds is the single path to CDS fasta, and aa is the path to output protein fasta
    # and returns some reports regarding number of sequences with internal stop codons
    
    seq <- readDNAStringSet(cds)

    # checks for final stop codons
    f <- lastChars(seq, 3) %chin% c("TAA", "TAG", "TGA")

    # if present, removes them
    seq[f] <- stri_sub(seq[f], 1, nchar(seq[f]) - 3)

    # does the translation
    tr <- suppressWarnings(translate(seq, if.fuzzy.codon = "solve"))

    # find if stop codons are present
    stops <- stri_count(tr, fixed = "*")

    # we only write the sequences without stop codons
    writeXStringSet(tr[stops == 0L], aa)
    
    # we return a table on the presence of stop codons for each CDS
    data.table(
        cds = basename(cds),
        nseq = length(seq),
        stops = sum(stops > 0)
    )
}


# where the BUSCO CDS fasta files will go (one file per species), this is done by the function below
dir.create("CDS")

# we obtain the BUSCO CDS sequence for all species with 10 parallel jobs
m <- mclapply(
    X = folders,
    FUN = getBUSCOs,  
    type="cds",
    mc.cores = 10,
    mc.preschedule = F
)

# we obtain the BUSCO proteins sequence for all species with 10 parallel jobs. Also run return paths aa
#Actually, I am not going to use this, Jean said it's saffer to translate CDS, instead of taking aa generated by busco
#aa <- mclapply(X = folders, FUN = getBUSCOs, type="aa",  mc.cores = 10, mc.preschedule = F)


# we list fastas of BUSCO CDS generated by the above command
CDSfiles <- list.files("CDS", pattern = ".CDS.fas", full.names = T)
#we write names aa
aaFiles <- gsub(".CDS.fas", ".translated.fas", CDSfiles)

# to avoid redoing translations that were already done if the job needs to be relaunched
f <- !file.exists(aaFiles)


# we translate the BUSO CDS into proteins in parallel jobs is done by the function below
# this functions writes a single fasta for all the CDS of a given species
# and returns a summary about the presence of internal stop codons (which may happen)
translationSummary <- mcMap(
    f = translateCDS,  #ATTENTION: prot are different that what is generated by BUSCO!!!
    cds = CDSfiles[f],
    aa = aaFiles[f],
    mc.cores = 10,
    mc.preschedule = F
)

# we stack the reports about stop codons generated in a single table
translationSummary <- rbindlist(translationSummary)


# we add a column for asembly names, which we extract from sequence names
translationSummary$assembly = sub(".CDS.*", "", translationSummary$cds)

write.table(translationSummary, "translationSummary.txt", col.names=T, row.names=F, quote=F, sep='\t')

#########################################
######## Step 1 bis : Make pairs ########
#########################################
# To reduce the workload, two tips can be used:
  # - one might decide not to look for HTT between very related species (eg <40Myrs)
  # In this case, one shouldn't do TE similarity searches for these pair of species eitehr (cf script 4)
  # - there is no need to measure dS in all possible species pairs between two large clades (only a subset of species is used)

# For this, we need the tree of the dataset to have information about time of divergence
tree <- read.tree("../datasetTree.nwk") 

# the matrix of divergence time between every tip (as a row and column index)
distMat <- cophenetic(tree) #ATTENTION cophenetic calculate additive distance (so twice more the divergence time since their common ancestor)

# the matrix listing the MRCA of every tip (as a row and column index)
mrca <- mrca(tree)
mrca[upper.tri(mrca)] <- NA #Thanks to mmseq rbh, we do not need reciprocal search, so put NA on half the matrix

# we turn these matrices into a data table for every species pair
pairs = setNames(data.table(as.table(distMat), as.vector(mrca)), 
                 c("sp1","sp2","divTime","mrca"))
pairs = filter(pairs, !is.na(mrca)) #Remove NA

# Although the tips have species name, our files have the assembly names. Read metadata to get information:
meta <- fread("../metadata.tbl")
meta$species = gsub(' ', '_', meta$species) #write species like in the tree

#add name assemblies in the table of pair of species
pairs <- left_join(pairs, select(meta, c("assembly", "species")), by=c("sp1" = "species")) %>% left_join(.,  select(meta, c("assembly", "species")), by=c("sp2" = "species"),  suffix=c(".1",".2"))

#remove version of assemblies because not written in file names
pairs$assembly.1 = sub("[.].*","", pairs$assembly.1)
pairs$assembly.2 = sub("[.].*","", pairs$assembly.2)

# we list species for which we have translated BUSCO (at least one did not have annotated genes)
aaFiles <- list.files("CDS", pattern = "translated.fas", full.names = T)
aaAssemblies <- sub("[.].*","", basename(aaFiles))

# we only retain pairs of species that have BUSCO genes
pairs <- pairs[assembly.1 %chin% aaAssemblies & assembly.2 %chin% aaAssemblies] #ATTENTION make sure that still as many row

#We save this because it will be very usefull afterward 
pairsALL = pairs 
fwrite(pairsALL, "pairsALL.txt", sep='\t')
# This table contains all the possible pairs of species, for which we looked for HTT or not, including species against itself
# Its nrow = ((nSpecies * nSpecies - nSpecies) / 2) + nSpecies, where nSpecies is the number of species in the dataset

# then make output names
pairs[, out := stri_c("dS/mmseqRBH/", assembly.1, "_vs_", assembly.2, ".out")]

# then the input names
pairs[, aa1 := stri_c("CDS/", assembly.1, ".translated.fas")]
pairs[, aa2 := stri_c("CDS/", assembly.2, ".translated.fas")]

# for clades older than 250 My (=500 Myrs with cophenetic), which are very large, we do not compute dS between all species as it is overkill.
# We will select one species per smaller "young" clades of <=30 My, based on the number of BUSCO detected in its genome

# we retrieve all these young clades with their species (see function in HTvFunctions.R)
youngClades <- cladesOfAge(tree, 30, withTips = T) 

#Add name assemblies for after
youngClades = left_join(youngClades, select(meta, c(assembly, species)), by=c("tip"="species"))
#remove version of assemblies because not written in file names
youngClades$assembly = sub("[.].*","", youngClades$assembly)

# number of translated BUSCO cds for these species
youngClades = left_join(youngClades, translationSummary)
youngClades = mutate(youngClades, nCDS = nseq-stops) %>% select(., -c(nseq, stops))

# we put species that have more BUSCO genes on top
setorder(youngClades, -nCDS)

# so we may ignore all the other species for each clade
toIgnore <- youngClades[duplicated(node), tip]


#Retained only these pairs for the similarity search
pairs <- pairs[
  # For clades older then 500, keep just one one species per young clades
  (divTime >= 500 &
    !sp1 %chin% toIgnore &
    !sp2 %chin% toIgnore) |
  # For other clades, do not keep those which are too related  
    (divTime > 80 & divTime < 500)]


# If some time of divergences were arbitrary set to <40Myrs, but might biologicaly be >40Myrs, we want to keep these pairs (keepPai==T in pairsDivTimeUnknown)
  # Read file in which are the pair of species for which the time of divergence was arbitrary chosen :
pairsDivTimeInconnue = fread("../pairsDivTimeUnknown", sep = " ", header = T, select= c("sp1_name", "sp2_name", "keepPair"), col.names = c("sp1", "sp2", "keepPair"))
keepPairs = filter(pairsDivTimeUnknown, keepPair=="T") %>%
  select(., c("sp1", "sp2"))
keepPairs = inner_join(keepPairs, pairs)

# Add those pairs to our table
pairs = rbind(pairs, keepPairs)

write.table(pairs, "pairsToMMseqRBH.txt", row.names = F, quote=F, sep='\t')

# Give a batch numvber to all pairs we need to run
# We do 10 batch
pairs[!file.exists(out), batch := rep(1:10, length.out = .N)] 

# This file contains all the pairs of species for which we want to do a similarity search of their BUSCO
write.table(pairs, "pairsToMMseqRBH.txt", row.names = F, quote=F, sep='\t')


# Function that create a file per batch, containing several similarity search
# Here, htis file is made to run on a server that does not use slurm
mmseqRBH <- function(job, pairs) {
    pairs = filter(pairs, batch==job)
    tb = data.frame(command = rep("mmseqs easy-rbh", nrow(pairs)), aa1=pairs$aa1, aa2=pairs$aa2,  out=pairs$out, para= rep("tmp --threads 20 & PID=$!", nrow(pairs)))
    write.table(tb, paste0("commandMMseqsRBH_job", job, ".txt"), row.names=F, col.names=F, quote=F, sep=' ')
}

for (i in 1:10){
    mmseqRBH(i, pairs)
}

#Then add wait between each line, so it wait that one line is done before running the next
#sed -e 's/$/\nwait $PID/' -i commandMMseqsRBH_job*.txt

#Lauch each job manually:
#nohup bash commandMMseqsRBH_job1.txt &> nohup_job1 &

#Since the jobs are split into 10 files, 10 jobs can therocally be run at the same time


#########################################
###### Step TWO : Compute dS & dN ######
#########################################

# Computing and processing of dS & dN values are done by the script 05bis
system("Rscript 05bis_coreGenesdNdS.R 30")

# establishing dS distribution between sister clades from dS obtained 

# we import all results off dnDs computations
dS <- fread("dS/dNdS/all.dNdS.txt", header = T)
  #Quey and subject under the form bedBusco:assembly
  #ALso columns for coordiantes hit, length of alignement 
  #And of for columns for dS and dN

# extracts species names from sequence names
dS <- mutate(dS, assembly.1 =  sub(".*:", "", query), assembly.2 = sub(".*:", "", subject))

# extracts bed busco
dS <- mutate(dS, bedBusco.1 =  sub(":GCA.*", "", query), bedBusco.2 = sub(":GCA.*", "", subject))

# We need a correspondance between the name of the busco file and the id of the fasta
# Indeed, busco name is the name of the file
  # e.g. GCA_017591415_busco/run_actinopterygii_odb10/busco_sequences/single_copy_busco_sequences/100028at7898.fna: 
  #100028at7898 is the name of the busco while the id of the fasta is scfAssembly:start-end

# To obtain this correspondance, use the script 05tris-FindNamesBusco.sh, whose output is buscoNames-ids.txt2

buscoNamesIds <- fread("buscoNames-ids.txt2", header=F, col.names=c("assembly", "busco", "bedBusco")) 
buscoNamesIds <- select(buscoNamesIds, -assembly)

#ATTENTION sometimes 2 busco end up with the same bed coordinates ...
viewCorrection <- filter(buscoNamesIds, bedBusco %in%  buscoNamesIds[duplicated(buscoNamesIds$bedBusco),]$bedBusco) 

#If several busco have the same coordinates, I keep only one
correction <- data.frame(buscoOld =  viewCorrection[seq(1, nrow(viewCorrection), by=2), ]$busco, 
  buscoNew = viewCorrection[seq(2, nrow(viewCorrection), by=2), ]$busco, 
  stringsAsFactors=FALSE) %>% 
  distinct()

for(i in 1:nrow(correction)){
  buscoNamesIds <- mutate(buscoNamesIds, busco = ifelse(busco==correction[i,]$buscoOld, correction[i,]$buscoNew, busco))
}
buscoNamesIds <- distinct(buscoNamesIds)

dS <- left_join(dS, buscoNamesIds, by=c("bedBusco.1"="bedBusco")) %>% 
  left_join(., buscoNamesIds, by=c("bedBusco.2"="bedBusco"), suffix=c(".1", ".2"))

# determining the MRCA (node/clade number) for all pairs of species, 
# which we will use to establish the distribution of dS values for each clade
# so we add the MRCA (clade) of each species pair as a new column.

# Info about mrca of each pair is in this file
pairs = fread("pairsToMMseqRBH.txt")

# Join it to add mrca in the table of dS
dS = left_join(dS, pairs, by = join_by(assembly.1, assembly.2))

#Note: it is normal that different buso names have hits, 
# because a same gene does not have the same name if it comes from a different busco database
dS[,mean(busco.1==busco.2)] #for example, here I had 14% of the hits

#Check that when same db for busco, busco have the same names
#299 = "Hydrodamalis_gigas" "Dugong_dugon" vs "Trichechus_manatus_latirostris"
filter(dS, mrca==299 & busco.1==busco.2) %>% nrow #yes (most of the hits: 8247)

# We are going to need a single value of dS per mrca and
# For this, we will calcualte a median

# Yet, we want the median to be calculated on good alignements,
#  so we keep only those >=100aa
dS100 = filter(dS, alnLength>=100)

# We also want to remove busco that look weird
# A weird "busco" would be a "busco" which is not found in many genomes 

# For this, we have to count the number of occurence for each busco per busco db
# Get info about what busco db was used on what genome:
db_busco = fread("metadata.tbl", select=c("assembly","dbBusco"))
db_busco$assembly = sub("\\..*", "", db_busco$assembly)
dS100 = left_join(dS100, db_busco, by=c("assembly.1"="assembly")) %>% 
  left_join(., db_busco, by=c("assembly.2"="assembly"), suffix=c(".1", ".2"))

#Add a column with a unique number per hit
dS100 = mutate(dS100, hit = seq(1:nrow(dS100)))

#To count correclty, make a new table (dS100_bind) in which 
# we write busco names of subject and of query in the same colomn:
dS100.q = select(dS100, c(hit, mrca,busco.1, dbBusco.1))
colnames(dS100.q) = c("hit", "mrca", "busco", "dbBusco")
dS100.s = select(dS100, c(hit, mrca,busco.2, dbBusco.2))
colnames(dS100.s) = c("hit", "mrca", "busco", "dbBusco")
dS100_bind = rbind(dS100.q, dS100.s)

#We don't keep duplicated lines, 
# to not count as 2 occurences a hit that involves the same busco in query & subject
#eg if hit 1 is buscoA in species1 vs buscoA in species2, counts as 1 occurence for buscoA
#eg if hit 2 is buscoA in species1 vs buscoB in species2, also counts as 1 opccurence
dS100_bind_uniq = distinct(dS100bind)

#For each busco db, look how many times each of its busco occurs in the table of hits
dS100_bind_uniq_summarize = dS100_bind_uniq %>% group_by(dbBusco, busco) %>% summarize(n=n())
  
#Get the 5% quantile & 95% quantile of these occurences
quantile_filter = dS100_bind_uniq_summarize %>%
  group_by(dbBusco) %>%
  summarize(q5 = quantile(n, 0.05), q95 = quantile(n, 0.95))
dS100_bind_uniq_summarize = left_join(dS100_bind_uniq_summarize, quantile_filter)

#We keep busco between 5% and 95% of the quantiles 
# + we set an absolute threashold with at least 30 occurences 
keep_busco = filter(dS100_bind_uniq_summarize, n>q5 & n<q95 & n>=30)
dS100_filter = filter(dS100, busco.1 %chin% keep_busco$busco &  busco.2 %chin% keep_busco$busco)

write.table(dS100_filter, "dS/dNdS/dS100AA_filterBusco.txt", quote=F, row.names=F, na="NA", sep='\t') 

# One can look at the distribution of what was kept or removed:
# (some adjustment might be necessary depending on the data, mostly regarding the absolute threashold of 30)
pdf("hist_distributionNbBusco_perDBbusco.pdf")
par(mfrow=c(3,2)) 
# One panel per db busco used
for(i in unique(dS100_bind_uniq_summarize$dbBusco)){
  dt = filter(as.data.frame(dS100_bind_uniq_summarize), dbBusco==i)
  hist(dt$n, breaks=1000, main=i, cex.main=0.8, xlab="# of occurence")
  abline(v=quantile(dt$n, 0.05)[[1]],col="red",lwd=0.5, lty=3)
  abline(v=quantile(dt$n, 0.95)[[1]],col="red",lwd=0.5, lty =3)
}
dev.off()


# Now that we filtered out the weird busco and the alignements too short, 
# we can calcualte the median dS with the hits that left

# For this, we need to have all species at the right of the node in the same column 
# and all species at the left of the node in the other column

#For this, obtain children of nodes:
edges <- as.data.table(tree$edge)

dS100_filter = mutate(dS100_filter, inv = NA, inv2=NA)

#Go trough each mrca
for(i in unique(dS100_filter$mrca)){
    print(paste("mrca:", i))

    #Get children of this mrca
    edge_sub = filter(edges, V1 == i)

    # NOTE: we could use mrca in dS100_filter, 
    # but the problem is there are not all there
    # since we didn't do similiarty searches for all pair of species.
    # This would be a problem; we would end up with some empty children_R and children_L 
    # This why we are using the table pairsALL, that contains all possible pair of species

    children_all = filter(pairsALL, mrca==i)
    children_all = unique(c(children_all$sp1, children_all$sp2))

    #species on the right of the node:
    children_R = filter(pairsALL, mrca == edge_sub[1,]$V2)
    children_R = unique(c(children_R$sp1, children_R$sp2))

    #species on the left of the node:
    children_L = filter(pairsALL, mrca == edge_sub[2,]$V2)
    children_L = unique(c(children_L$sp1, children_L$sp2))
      
    #Fill both children_R & children_L but it should be the same answer
    dS100_filter = mutate(dS100_filter,
      inv = ifelse(mrca != i, inv, ifelse(sp1 %chin% children_L, "FALSE", "TRUE")), 
      inv2 = ifelse(mrca != i, inv2, ifelse(sp2 %chin% children_R, "FALSE", "TRUE"))
      )
}

#check that all inv == inv2
filter(dS100_filter, inv!=inv2) #should be empty

#Now, inverse those who need to be:
toInvert = filter(dS100_filter, inv==TRUE)
colnames(toInvert) = c("subject", "query", "sStart", "sEnd", "qStart", "qEnd", "dS", "dN", "alnLength", "assembly.2", "assembly.1", "bedBusco.2", "bedBusco.1", "busco.2", "busco.1", "pairR", "pairF", "sp2", "sp1", "divTime", "mrca", "out", "aa2", "aa1", "dbBusco.2", "dbBusco.1", "hit", "inv2", "inv")
toInvert = select(toInvert, colnames(dS100_filter))
toNotInvert  = filter(dS100_filter, inv==FALSE)
# Put back together those we want to invert and the others
dS100_filter = rbind(toInvert, toNotInvert)

write.table(select(dS100_filter, -c( "hit", "inv", "inv2")), "dS/dNdS/dS100AA_filterBusco_organized.txt", quote=F, row.names=F, na="NA", sep='\t')

# We can now calculate the median dS for each busco for each mrca
# We choose to not take into account values of 9.999, which means that it was saturated
# NOTE: for next studies, one should also remove dS<0 (it means it couln't be calculated)
dS100_filter_median = filter(ds100_filter, dS<9) %>% 
  group_by(mrca, busco.1) %>% 
  summarize(
    dSmedian = median(dS), 
    dNmedian = median(dN), 
    n = n(), #number of values used for this calculation
    alnLengthMedian = median(alnLength)
  )

# Now that we have one median dS per busco gene, we obtain a distribution of median dS per mrca

#Add info about divTime
dS100_filter_median = select(pairs, c(mrca,divTime)) %>% 
  distinct() %>% 
  left_join(., dS100_filter_median)

# Save this very important table that we will use to decipher whether a hit of TE is from HT or HV
write.table(dS100_filter_median, "dS/dNdS/dS100AA_filterBusco_median.txt", quote=F, row.names=F, na="NA", sep='\t')



# making a figure to illustrate the rate of synonymous molecular evolution 

# Delete pairs for which we don't know DivTime
# In my case, the one mrca for which we didn't know divTime was set <80Myrs:
mrcaPif = c(unique(filter(dS100_filter_median, divTime<80)$mrca)) 
dS100_filter_median_confidentDiv = filter(dS100_filter_median, !mrca %in% mrcaPif)

# We generate divergence time classes within which we compute average dS (all species pairs considered).

# We use sqrt() to have shorter intervals for low divergence times
timeRanges <- seq(0, sqrt(max(dS100_filter_median_confidentDiv$divTime, na.rm = T) + 1), length.out = 20)^2

# assigns divergence times to the classes
dS100_filter_median_confidentDiv <- mutate(dS100_filter_median_confidentDiv, 
  range = .bincode(divTime, timeRanges))
dSDistrib <- dS100_filter_median_confidentDiv[
    dSmedian < 9 & # ignoring dS values ≥ 9 (oversaturated)
    dSmedian>=0 & # ignoring dS valuees <0 (couln't be calculated)
    !is.na(divTime), #ignoring when na for divTime (we are not sure of the divergence time for that clade)
    .(time = mean(divTime), 
      dS = weighted.mean(dSmedian, w = alnLengthMedian)), 
    by = range 
]


pdf("dSdistribution_of_BUSCO_overTime.pdf")
p <- dSDistrib[, plot(
    x = time,
    y = dS,
    xlab = "Divergence time (My)",
    ylab = "dS",
    pch = 16,
    col = "darkgrey"
)]

# fits a curve assuming initial linear correlation between dS and divergence time (coeff a) and saturation at dS = max
fit <- nls(dS ~ max * a * time / (max + a * time),
    data = dSDistrib,
    start = list(a = 0.01, max = 3)
)

# to add a smooth curve corresponding to the model, we generate 150 divergence time values (x axis)
x <- dSDistrib[, seq(0, max(time), length.out = 150)]
lines(x, predict(fit, list(time = x)))
dev.off()
