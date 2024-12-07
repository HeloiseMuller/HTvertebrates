##########################################################
####      This stage prepares the clustering of       ####
###   We will use 2 independent clustering methods :   ###
##         - per pair of species 
##         - per clade
##########################################################

library(data.table)
library(Biostrings)
library(dplyr)
library(stringi)

path = "~/Project/"
setwd(path) 

source("HTvFunctions.R")

#IMPORTANT
#Choose a method
method <- "perClade" 
#method <- "perPairs"

##################################################
##### Set parameters depending on the method chosen
##################################################

if(method=="perClade"){
    out <- "TEs/clustering/selfBlastn_perClade"
    key <- "mrca"
} else {
    out <- "TEs/clustering/selfBlastn_perPairs"
    key <- "assembly"
}

##################################################
######### Read intput = hits to clsuster #########
##################################################
httHits <- fread("TEs/clustering/dcMegablast/occ200dS05_dcMegablast.txt")

##################################################
# STEP 1 : Get copies to self blast
##################################################

copies <- select(httHits, c(ID.1, assembly.1, rep_superF.1, mrca)) %>% 
    rename(., ID.2 = "ID.1", assembly.2 = "assembly.1") %>%
    rbind(., select(httHits, c(ID.2, assembly.2, rep_superF.1, mrca))) %>% 
    rename(., ID = "ID.2", assembly = "assembly.2") %>%
    distinct()

stats_copies <- copies %>%
    #Group by assembly or mrca (in key) depending on the method chosen
    group_by(across(all_of(key)), rep_superF.1) %>%
    #Count the number of copies in each group
    summarize(nbCopy = n())

#The group with the most copies:
maxCopies <- max(stats_copies$nbCopy)
print(paste("Set max_target_seqs >", maxCopies)) #IMPORTANT

    
#We split per assembly or clade (in key) depending on the method chosen
copies_split <- split(copies, copies[,..key])

#Create output directory
dir.create(paste0(out, "/copies/"), recursive = TRUE)

#Function that take a tb in input, and extrat its IDs to a fasta
extractID_f <- function(tb, dir, superF){
#Input is a table woth 4 columns: ID, assembly, rep_superF.1, mrca
#Assembly MUST be th same for all lines of the table
#Only the 1st 2 columns are necessary to use this function
    #Get assembly name
    assembly = unique(tb$assembly)
    #Write in a file the list of ID we want to extract from fasta
    write.table(select(tb, ID), paste0(dir, superF, "_", assembly, ".bed"), quote=F, row.names=F, col.names=F)
    #Extract these ID from fasta
    system(paste0("seqtk subseq TEs/clustering/dcMegablast/copies/", assembly, "_", superF, ".IDs.fasta ", dir, superF, "_", assembly, ".bed > ", dir, superF, "_", assembly, ".IDs.fasta"))
    #Remove the list of ID (not usefull anymore)
    system(paste0("rm ", dir, superF,  "_", assembly, ".bed"))
}


#Put in a same file copies we want to blast
lapply(copies_split, function(x){

    #output for this one :
    dir <- paste0(out, "/copies/", key, "_", unique(x[,..key]), "/")

    #If there is already an output, does do it again
    if(!dir.exists(dir)){ 

        #Create output
        dir.create(dir)

        #Split per super family
        per_superF = split(x, x$rep_superF)
    
        lapply(per_superF, function(y){
            # Get superfamily name
            # But don't keep DNA/, LTR/n etc because "/" would  be interprated like a new directory
            superF <- sub(".*/", "", unique(y$rep_superF))

            if(method=="perClade"){
                #When method = clade, we have several assembly for a same clade
                #So we need to look for copies which are in different files
                per_assembly = split(y, y$assembly)
                #apply the function on each assembly
                lapply(per_assembly, extractID_f, dir = dir, superF = superF)
                #concantenate all the ourputs of this mrca-superF in a singe file
                system(paste0("cat ", dir, superF, "_*", ".IDs.fasta > ", dir, superF, ".fasta"))
                #we can now delete the fasta we concatenated
                system(paste0("rm ", dir, superF, "_GCA*.IDs.fasta"))
            } else {
                #When method = per pair of species, we can directly apply the function
                extractID_f(y, dir = dir, superF = superF)
            }   
        })
    }
})

##################################################
# STEP 2 : Make database 
##################################################


dir.create(paste0(out, "/out/"))
dir.create(paste0(out, "/tmp/"))

#Funciton to make a database on a fasta file
makeDB <- function(files){
  system(paste0("makeblastdb -in ",  dirCopies, files , ".fasta -dbtype nucl"))
  return(NULL)
}

##################################################
# STEP 3 : Generate command lines to launch blast searches
##################################################

#Make an emtpy table telling path for input and output
toDo_all <- data.table(file = as.character(), dirCopies = as.character(), dirOut = as.character())

for(i in unique(copies[, ..key])){
  dirCopies <- paste0(out, "/copies/", key, "_", i, "/")
  dirOut <- paste0(out, "/out/", key, "_", i, "/")
  dir.create(dirOut)

  # we import the table we need to determine the blast searches to launch
  searches <- list.files(dirCopies, pattern=".fasta$")
  searches <- sub(".fasta", "", searches)

  # in case jobs need to be relaunched, we list output files to avoid redoing already completed searches
  done <- list.files(dirOut, pattern="_selfBlastn.out")
  done <- sub("_selfBlastn.out", "", done)
  toDo <- searches[!searches %in% done]
  
  #We fill out table toDo_all
  toDo_all <- rbind(toDo_all,  data.table(file = toDo, dirCopies = rep(dirCopies, length(toDo)), dirOut = rep(dirOut, length(toDo))))
  
  #Run the fucntion makeDB on each of these files, using dirCopies for output for all of them
  mcMap(f = makeDB, files = as.list(toDo), mc.cores = 10, mc.preschedule = T)
}

#Write with one line per blast
sink(paste0("commands_selfblast_", method, ".sh"))
paste0("blastn ",
    "-task dc-megablast ",
    "-query ", toDo_all$dirCopies, toDo_all$file , ".fasta ",
    "-db ", toDo_all$dirCopies, toDo_all$file , ".fasta ",
    "-max_target_seqs ", (maxCopies+1),
    " -max_hsps 1 ",
    "-outfmt '6 qseqid sseqid pident length qstart qend sstart send bitscore'",
    "-num_threads 10 | awk '{if ($4>=100 && $1<$2) print $0}' ",
    "> ", toDo_all$dirOut, toDo_all$file, "_selfBlastn.out & PID=$!"
    )
sink()

#IMPORTANT
#ATTENTION add this lines thanks to bash:
# sed -i 's/.*"blastn/"blastn/' commands_selfblast_perClade.sh
# sed -i 's/"//g' commands_selfblast_perClade.sh
# sed -e 's/$/\nwait $PID/' -i commands_selfblast_perClade.sh


#Then i can split the file to run it in parallel (here 7014 lines, counting the waits)
sed -n '1,1753p' commands_selfblast_perClade.sh > commands_selfblast_perClade.sh1
sed -n '1754,3507p' commands_selfblast_perClade.sh > commands_selfblast_perClade.sh2 
sed -n '3508,5260p' commands_selfblast_perClade.sh > commands_selfblast_perClade.sh3 
sed -n '5261,7014p' commands_selfblast_perClade.sh > commands_selfblast_perClade.sh4 


#NOTE some files have copy of huge size that stuck the jobs
#One might one to remove these copies

##################################################
# STEP 4: run blast searches
##################################################

#This step consist in running the bash scripts generated above
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!! ATTENTION Make sure there is a WAIT between each line!!! VERY IMPORTANT !!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Otherwise could kill the server !


##################################################
# STEP 5: Remove copies too long
##################################################

#Sometimes a copy is anomaly long (eg 167kb) which lead to blast that do not work
#In such a case, remove these copies from TEs/clustering/dcMegablast/occ200ds05_dcMegablast.txt and save it
