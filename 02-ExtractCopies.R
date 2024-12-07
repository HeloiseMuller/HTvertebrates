################################
###### Step 0: annotate TE ######
################################

# On each assembly, run RepeatModeler this way :
  # BuildDatabase -name $assembly -engine rmblast genomes/$i
  # RepeatModeler -database $assembly -LTRStruct -pa 1 -ninja_dir /opt/NINJA-0.95-cluster_only/NINJA/
  # To do in the directory RepeatModeler
  # Before running RepeatModeler on all genomes, check that the 1st output is right: should have IDs begining by >ltr, >rnd-1, ..., >rnd-6

# Several assemblies could end up with the same TE name
  # Give a unique name to each TE :
  # echo RepeatModeler/${assembly}-families.fa.classified | cut -f1 -d'-'`; sed 's/>/>'$assembly':/' $i > RepeatModeler/${assembly}-families.fa.classified.renamed

# Then, concatanate all this TE :
  # cat RepeatModeler/*families.fa.classified.renamed > DB_clustered/TE_cat.fasta
  # Optinnaly, one can add other TE consensus > DB_clustered/TE_cat_all.fasta

# Keep only classified TE
  # Firstly look at the different super families that exist in the DB:
    # grep ">" DB_clustered/TE_cat_all.fasta | cut -f2 -d'#' | awk '{print $1}' | sort | uniq
    # Write as many grep -v as there are superfamilies we want to remove :
      # grep  ">" DB_clustered/TE_cat_all.fasta | grep -v 'LTR/Unknown' | grep -v '?' | grep 'DNA/\|LINE/\|Helitron\|DNA_virus\|LTR/' | sed 's/^>//' > DB_clustered/TE_cat_all.filtered.lst
    #seqtk subseq DB_clustered/TE_cat_all.fasta DB_clustered/TE_cat_all.filtered.lst > DB_clustered/TE_cat_all.fasta_filtered.fasta

# Keep only TE >= 300bp
  # bioawk -c fastx '{ if(length($seq) > 300) { print ">"$name; print $seq }}' DB_clustered/TE_cat_all.fasta_filtered.fasta >  DB_clustered/TE_cat_all.fasta_filtered.long.fasta 

# Clusterize this common DB with MMseq2 
  # mmseqs easy-cluster DB_clustered/TE_cat_all.fasta_filtered.long.fasta DB_clustered/TE_cat_all.filtered.long.clustered  tmp --min-seq-id 0.8 -c 0.8 --cov-mode 1 -v 2 --threads 10 --cluster-reassign

# Run RepeatMasker on each the assemblies with this common database 
  # RepeatMasker -nolow -no_is -norna -engine rmblast -pa 1 -gff -lib TDB_clustered/TE_cat_all.filtered.long.clustered_rep_seq.fasta -dir RepeatMasker/
  # v4.0 has a bug when option dir --> use at least v4.1.0
  # Also create a directory RepeatMasker/copies


################################
################################


library(data.table)
library(dplyr)
library(stringr)
library(stringi)
library(tidyr)
library(parallel)

path = "~/Project/"
setwd(paste0(path, "RepeatMasker"))

source("../HTvFunctions.R")

################################
###### Step 1: list files ######
################################

RMout <- list.files(path=".", pattern = ".out$", full.names = T, recursive = F)

#Get name assemblies
RMouts <- data.table(RMout, assembly = stri_extract(RMout, regex = "[A-Z]+_\\d+.\\d"))

# we add a column for the names of future output files 
# (compressed fastas of TE copies for each species)
RMouts[, out := stri_c("copies/", assembly, ".TEs.fasta")]

# we avoid redoing some work in case the script needs to be relaunched
RMouts <- RMouts[!file.exists(out)]

genome <- list.files("../genomes", pattern = "_genomic.fna$", full.names = T, recursive = T)

genomes <- data.table(genome, assembly = stri_extract(genome, regex = "[A-Z]+_\\d+.\\d"))

# as our function to extract copies will use both these files
# we put these file names in a table, where each row corresponds to a species 
m <- merge(RMouts, genomes, by = "assembly", all.x = T)


# we will process bigger genomes first, not make better use of the CPUs
# (else, one single job may be running at the end for the largest genome)
sizes <- file.size(m$genome)
m <- m[order(sizes, decreasing = T)]

# we split the filename table by row (species) to send to the function below in parallel
# we don't call mcMap() on the table because it apparently has issues returning tables
fileList <- split(m, 1:nrow(m))

#ATTENTION check that genome exist

################################
###### Step 2: Extract TE ######
################################

bedtools <- function(fas,
                  bed,
                  out,
                  ex = "bedtools getfasta",
                  formated = F) {
  # calls bedtools to return sequence from a fasta fas, based on bedfile bed.
  # Write to fasta file out (if specified) or return to R (if not). 
  # if formatted is TRUE, return sequences as a DNAStringset (if not, returns a character vector 
  # corresponding to the fasta) "ex" specifies how bedtools getfasta subseq should be executed

    #e.g. bedtools getfasta -fi ../genomes/GCA_902825295.1_Laufergenomev3_genomic.fna -bed GCA_902825295.1_Laufergenomev3_genomic.fna.out.bed -name  > GCA_902825295.1_Laufergenomev3_genomic.fna.out.fa 
  
  if (missing(out)) {
    seqs <- system(paste(ex, "-fi", fas, "-bed", bed, "-name"), intern = T)
    
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
    system(paste(ex, "-fi", fas, "-bed", bed, "-name >", out))
  }
}


extractCopies <- function(assembly, RMout, out, genome, minLen = 300L) {
    # extracts TE copies on length ≥ minLen of a given species (sp) from its genome (genome) to gzipped fasta file (out), 
    #based on a gzipped repeat masker gff (gff), write also a report of the percentage of consensus sequences in "cons" 
    # (the repeat modeler fasta file) that have masked copies as a filename ending by ".percentMasked.txt".
    # this function returns metrics on the TE compositions as a data.table (see end of function)
       
    # We read all RM outputs, execpt for the lines that contanin "*"
    # Copies with "*" are copies which are included in a higher-scoring match
    RMout <-  fread(cmd=paste("grep -v *$ ", RMout), skip=3, header=F, fill=T)

    #Rename columns
    colnames(RMout) = c("SW_score", "percDiv", "percDel", "percIns", "querySequence", "query_posBegin", "query_posEnd", "query_left", "strand", "TEconsensus", "superF", "TE_posBegin", "TE_posEnd", "TE_left", "ID")

    #Remove rRNA
    RMout = filter(RMout, !str_detect(TEconsensus, "rRNA"))
 
    # selects TE copies that are long enough 
    selectedCopies <- filter(RMout, (query_posEnd-query_posBegin+1L) >= minLen)

    #Make bed to make fasta
    bed = select(selectedCopies, c("querySequence", "query_posBegin", "query_posEnd", "TEconsensus", "strand")) %>%  
        mutate(., query_posBegin = query_posBegin-1) #correct pos
    # Give a unque copy name to each copy: assembly-queryAssembly-startAssembly-endAssembly-strandAssembly-TEconsensus
    bed = mutate(bed, name = paste(assembly, querySequence, query_posBegin, query_posEnd, strand, TEconsensus, sep='-')) 
    bed_out = stri_c("copies/", assembly, ".TEcopies.bed")
    
    write.table(select(bed, -c(strand,TEconsensus)), bed_out, col.names=F, row.names=F, sep='\t', quote=F) 

    # imports sequences directly from bedtools (as a fasta string)
    bedtools(genome, bed_out, out)   #also possible to not give output to not save in on disk, but in variable 
    
    #Proportion of consensus found in this assembly
    prop <- mean(cons[,"TEconsensus"] %chin% selectedCopies$TEconsensus)
    # writes a report about the percent of consensuses that masked TEs
    write.table(data.table(assembly = assembly, prop, nCopies = nrow(selectedCopies)),
        stri_c("TEcomposition/", assembly, ".percentMasked.txt"), row.names = F, quote = F)

    # we return metrics on the TE composition of the genome, per superfamily, for selected copies
    res <- selectedCopies[, .(
        nCopies = .N,
        nCons = length(unique(TEconsensus)), #ERROR unique() s'applique seulement aux vecteurs
        bp = sum(query_posEnd - query_posBegin + 1L)
    ), by = superF]

    # adds assembly name in a new column
    res[, assembly := assembly]
    res
}

#Read IDs of the fasta of TE consensus (ie the file we used to run RepeatMasker)
cons <- system(paste("grep '>'", "../DB_clustered/TE_cat_all.long.clustered_rep_seq.fasta"), intern = T)
cons <- gsub(">", "", splitToColumns(cons, "#"))
colnames(cons) <- c("TEconsensus", "superF")
#N.B. Since some TE consensus that we added to our database are not well formated, 
#they don't have a superfamily name wheras they do correcponf to classifed TE
#This is why we will classify all our TE consensus of interest ourselves latter in the pipeline (at script 6)

TEcomposition <- mclapply(
  X = fileList, #will use the arguments of fileList as arguments of the function extractCopies
  FUN = function(files) do.call(extractCopies, files),
  mc.cores = 20,
  mc.preschedule = F
)


# we stack the data for all genomes in a single table
TEcomposition = rbindlist(TEcomposition)
write.table(TEcomposition, file = "TEcomposition/all.TEcomposition.txt", row.names=F, quote=F, sep='\t')

# we generate the supplementary dataset associated with the paper
# in which we show the composition over all genomes on a per-superfamily basis

TEcompo <- TEcomposition[, .(
    number_of_copies_300bp = sum(nCopies),
    total_nucleotides = sum(bp),
    number_of_consensus = sum(nCons)
), by = .(assembly = assembly, RepeatModeler_superfamily = superF)]


TEcompo <- TEcompo[number_of_copies_300bp > 0L]

fwrite(TEcompo, "TEcomposition/supplementary-data-TEcomposition_per_species.txt", sep='\t') 
