## %######################################################%##
#                                                          #
####            this script processes "raw"             ####
####              blast output between TEs              ####
####              of different species to               ####
####          select hits with sufficient pID           ####
#                                                          #
## %######################################################%##

# This script is run at stage 06-filterTEhits.R


library(data.table)
library(parallel)
library(stringi)
library(stringr)
library(dplyr)

path = "~/Project/"
setwd(path) 


# the only argument is the number of CPUs to use
args <- commandArgs(trailingOnly = TRUE)
nCPUs <- as.integer(args[1])

# where the filtered blast output files will go
dir.create("filtered_HitsTE/")
pairs <- fread("pairsToFilterQ5.txt") 

#List files with outputs 
files = list.files(path="mapCopies", pattern = "tab.filtered", full.names = T, recursive = T)
files1 <- str_match(files, "mmseq2_\\s*(.*?)\\s*_vs")[,2]
files2 <- str_match(files, "vs_\\s*(.*?)\\s*_evalue")[,2]
files = data.frame(assembly1 = files1, assembly2 = files2, out = files)
files$out = as.character(files$out) #data.frame() turned files into interger. Put it back to character

#in case i ran a pair twice, keep only one output
files = files[!duplicated(files$out),]

#Add name output in table pairs
pairs = inner_join(pairs, files, by=c("assembly1", "assembly2"))

#Add name output of the reverse mapping
colnames(files) = c("assembly2", "assembly1", "rev")
pairs = inner_join(pairs, files, by=c("assembly1", "assembly2"))

# paths of filtered output files to be generated
pairs = mutate(pairs, filtered = stri_c("filtered_HitsTE/", stri_extract(out, regex="mmseq.*")))

# this function filter the hits for two reciprocal searches:
filterHits <- function(out, rev, filtered, q) {
    
    # we create an empty data table to avoid a bug in fread() if a search output file is empty (no hit)
    search <- data.table()

    # to report the number of hits
    nr <- 0
    if (file.size(out) > 0) {
        # imports "forward" search between two species
        # we do not name columns as they wont be retained in the output
        search <- fread(
            out,
            sep = "\t",
            header = F,
            drop = c(5, 6, 11)
        )

    print("search read")
        nr <- nr + nrow(search)
        
        if(q!="FALSE"){
            # Keep hits with sufficient pID
            # ie hits whose pID > q where q is 1 - 0.5% quantile of dS of Busco
            search <- search[V3*100 > q, ] #I multiply by 100 to get a percetage (with mmseq 0<pID<1)
        }
    }

    # we do the same for the "reverse" search
    reverse <- data.table()
    if (file.size(rev) > 0) {
        reverse <- fread(
            rev,
            sep = "\t",
            header = F,
            drop = c(5, 6, 11)
        )
    print("reverse read")
        nr <- nr + nrow(reverse)
        if(q!="FALSE"){
            reverse <- reverse[V3*100 > q, ]
        }

        # reversing query and subject fields so we can concatenate forward and reverse hits, and remove reciprocal hits (amongst other things)
        reverse[, c("V1", "V2", "V7", "V8", "V9", "V10") := .(V2, V1, V9, V10, V7, V8)]
    }

    # if there is no hit, we exit
    if (nrow(search) == 0 & nrow(reverse) == 0) {
        return(NULL)
    }
    
    # we stack the forward and reverse hits
    search <- rbind(search, reverse)
    rm(reverse)

    # we put the best hits on top (column 12 is the score)
    setorder(search, -V12)

    # we retain the best hsp per pair of copies, which also removes reciprocal hits
    search <- search[!duplicated(data.table(V1, V2))]

    #Get assemby names
    mat1 <- stri_split(search$V1, fixed = "-", simplify = T)
    mat2 <- stri_split(search$V2, fixed = "-", simplify = T)
 
    search <- mutate(search, assembly1 = unique(mat1[,1]), assembly2= unique(mat2[,1]))
    
    colnames(search) = c("copie1", "copie2", "pID", "length", "qStart", "qEnd", "sStart", "sEnd", "score", "assembly1", "assembly2")


    # we write file of filtered hits (not returned by the function, as we are limited in the memory we can use here)
    write.table(
        x = search,
        file = filtered,
        col.names = F,
        row.names = F,
        quote = F,
        sep = "\t"

    )

    # instead we return the initial number of hits and number of retained hits
    c(basename(filtered), nr, nrow(search))
}

# applies the function with the requested number of CPUs
hitNumbers <- pairs[!file.exists(filtered), mcMap(
    f = filterHits,
    out,
    rev,
    filtered,
    q = q5,
    mc.cores = nCPUs,
    mc.preschedule = F
)] 

# we stack the reports into a single table
hitNumbers <- data.table(do.call(rbind, hitNumbers))

# and write them to file
write.table(
    x = hitNumbers,
    file = "filtered_HitsTE/filteredStats.txt", 
    col.names = F,
    row.names = F,
    quote = F,
    sep = "\t"
)

# we concatenate the output files of filtered hits. I can't do it with system() because too many files
system("cat filtered_HitsTE/*.out > filtered_HitsTE/all.quantile005score200.mmseq2.out")

# NOTE:It might not work if too many species in the dataset. 
# In this case, do the following, adapting the values of hean -n and tail -n in such a way that all files processed (in my case I had 30313 files)
    #cd filtered_HitsTE/
    #ls -lh  | grep "tab.filtered" | head -n15156 > cat1-15156
    #ls -lh  | grep "tab.filtered" | tail -n15157 > cat15157-30313

    #On each file:
    #awk '{print $9}' cat1-15156 > a 
    #mv a cat1-15156
    #cat cat1-15156 |  tr '\n' ' '  > cat1-15156.sh #put file names on the same line
    #sed 's/^/cat /'  cat1-15156.sh -i #add cat at the begining of the line
    #sed 's/$/ >  cat1-15156.quantile005score200.mmseq2.out/'  cat1-15156.sh -i #add name of the concatenated file at the end of the line
    #bash  cat1-15156.sh
        
    #once done all files:  cat cat*out > all.quantile05score200.mmseq2.out


