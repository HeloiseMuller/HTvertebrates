library(data.table)
library(dplyr)
library(stringr)
library(ape)
library(stringi)

path = "~/Project/"
setwd(path) 

###############################
### STEP ONE: We make pairs ###
###############################

#Read the tree of the dataset to get divergence time
    # Importantly, the tree should not have any polytomy.
    # The absolute time of  divergence is not that important, but the topology should be right
tree = read.tree("datasetTree.nwk") 
distMat <- cophenetic(tree)
distTbl <- setNames(data.table(as.table(distMat)), c("sp1_name","sp2_name","divTime"))

#NOTE: cophenetic calclulates additive distances (so twice more the divergence time since their common ancestor)
#Here I want ultrametric distance, so I divide by 2
distTbl$divTime = distTbl$divTime/2

#Read the table containing all the metadata (eg species name, habitat)
meta = fread("metadata.tbl")
meta$species = gsub(' ', '_', meta$species) #write species name like in the tree, ie Species_genus

#Empty dataframe that will cointain two colomns: species1 & species2  
selectedSpeciesPairs = data.frame(assembly1=character(0), sp1_name=character(0), group1=character(0), assembly2=character(0), sp1_name=character(0), group2=character(0), divTime = numeric(0))

##Go through each genome to fill selectedSpeciesPairs
for (i in 1:nrow(meta)) {
    #Filter 1 is because this genome might already be paired with previous genomes done previously in the loop
    #Filter 2 is because we do not want pair someone with itself
    paired = filter(meta, ! assembly %in% selectedSpeciesPairs$assembly1 & assembly!=meta[i,]$assembly) 
    selectedSpecies_i =  data.frame(assembly1 = rep(meta[i,]$assembly, nrow(paired)),
        sp1_name =  rep(meta[i,]$species, nrow(paired)),
        group1 = rep(meta[i,]$group, nrow(paired)),
        assembly2 = paired$assembly,
        sp2_name=paired$species,
        group2 = paired$group)

    #Add info about divergence time
    selectedSpecies_i = inner_join(selectedSpecies_i, distTbl)
    
    #Add these pairs to the ones of the previous genomes
    selectedSpeciesPairs = rbind(selectedSpeciesPairs, selectedSpecies_i)

}

# In order to reduce the workload, one can choose to not look at all for HTT among young clades
# Here we chose 40Myrs, based on Peccoud et al. 2017 (study on insects)
# But this time may vary depending on the taxa
selectedSpeciesPairs_filtered = filter(selectedSpeciesPairs,  divTime>=40)

# If the time of divergence between 2 species in uncertain, one should do the search in case both species are acutally not that related
# In my case, I wrote 30 Myrs for a pair of species for which I couldn't find any clue of its divergence. But I had to write something to have it in the tree  
# I saved a list of pairs for which I am not sure of divTime. The one woth keepPair==T could be older than 40Myrs
pairsDivTimeUnknown <- fread("pairsDivTimeUnknown", sep = " ", header = T)
keepPairs <- filter(pairsDivTimeUnknown, keepPair=="T") %>%
    select(., c("sp1_name", "sp2_name"))
keepPairsRev = keepPairs
colnames(keepPairsRev) = c("sp2_name", "sp1_name")
keepPairs = rbind(keepPairsRev, keepPairs)
rm(keepPairsRev)
selectedSpeciesPairs_filtered = rbind(selectedSpeciesPairs_filtered, inner_join(keepPairs, selectedSpeciesPairs))

# we will launch the longer blast searches first (those involving the bigger fasta files), to better use the CPUs NB Not totally true, run time also depends a lot on divergence time
selectedSpeciesPairs_filtered$size = file.size(stri_c("RepeatMasker/copies/", selectedSpeciesPairs_filtered$assembly1, ".TEs.fasta")) + file.size(stri_c("RepeatMasker/copies/", selectedSpeciesPairs_filtered$assembly2, ".TEs.fasta")) 

setorder(selectedSpeciesPairs_filtered, -size)

# Save the table contening all the pair of species between which we want to look for HTT
write.table(selectedSpeciesPairs_filtered, "pairsToSearch.txt", row.names=F, quote=F, sep='\t') 

###############################
### STEP TWO: We make db ###
###############################

#Create file that will allow us to run all those searches on the server genotoul
#In reality, we ran some searches on this server, and others on other servers (eg not shown)
sink("commandCreateDB_genotoul")
for(i in meta$assembly){
    if(i %in% selectedSpeciesPairs_filtered$asssembly1 | i %in% selectedSpeciesPairs_filtered$assembly2){
        print(paste("module load system/Miniconda3-4.7.10 && module load bioinfo/mmseqs2-v13.45111 && mmseqs createdb copies/", i, ".TEs.fasta copies/", i, ".TEs.fasta.dbm --dbtype 2", sep=""))
        print(paste("module load system/Miniconda3-4.7.10 && module load bioinfo/mmseqs2-v13.45111 && mmseqs createindex copies/", i, ".TEs.fasta.dbm tmp --search-type 3", sep=""))
    }
}
sink()
#then do: sed -i 's/.*"module/"module/' ../../commandCreateDB_genotoul 
#sed -i 's/"//g' ../../commandCreateDB_genotoul

# And run the job with :
# sarray -J makeDB --mail-type=FAIL commandCreateDB_genotoul 

###############################
### STEP THREE: We run mmseq ###f
###############################

# Here again I show only for genotoul.
# It might differ with other server depending on how it works

sink("commandMMseqs_genotoul")
paste("sbatch 04bis-command_mmseq2_model.sbatch ", selectedSpeciesPairs_filtered$assembly1, " ", selectedSpeciesPairs_filtered$assembly2, " ", sep="") #$1=querry, $2=target
paste("sbatch 04bis-command_mmseq2_model.sbatch ", selectedSpeciesPairs_filtered$assembly2, " ", selectedSpeciesPairs_filtered$assembly1, " ", sep="")
sink()

#then do: sed -i 's/.*"sbatch/"sbatch/' commandMMseqs_genotoul 
#sed -i 's/"//g' commandMMseqs_genotoul

# And run the job with :
sarray -J mmseq --mail-type=ARRAY_TASKS,FAIL  --%=20 commandMMseqs_genotoul  
#ATTENTION I cannot give such a big file (2500 works, 3000 does not)
#--> split it 

# Genotoul use slurm and have a special script array that allow to manage a list of jobs
# For other servers using slurm, one can use array
# For a server not using slurm, one can something similar to one I did for Genotoul, I add a wait $PID between each line
    # This way the next line will start only once the previous one is done
    # The file should be split in several files, so several jobs can run at the same time
