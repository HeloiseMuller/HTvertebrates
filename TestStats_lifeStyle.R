## %#######################################%##
#                                            #
####        Statistical test to           ####
####     test the effect of the lifestyle ####
####         on the number of HTT         ####
#                                            #
## %#######################################%##

# This can be run at any point after stage 11


library(data.table)
library(dplyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(stringr)

path = "~/Project/"
setwd(path) 

source("HTvFunctions.R")

##############################################
# We sampled a total of 19 taxonomic groups, each of which contain at least 1 independent transition of lifestyle

# Choose the number of taxa to use, comment the two others
    # All the taxonomic groups:
nbTaxa = 19
    # Exclude the 2 taxonomic groups composed of aquatic species only (fish & palaeoptera):
nbTaxa = 17
    # Exclude the same 2 taxonomic groups, but also the ones without fully aquatic species:
nbTaxa = 10
    # Keep fish
nbTaxa = 11


##############################################
# Inputs:

httHits <- fread("occ200HitGroup_perPair.txt")
evaluation <- fread("TEs/clustering/hitGroupStats_perPair.txt")

httHits_filtered <- filter(httHits, hitGroup %in% evaluation[retained==TRUE,]$hitGroup)

#Save it for publication
httHits_filtered_toSave <- httHits_filtered
setnames(
    x = httHits_filtered_toSave,
    old = c("copie1", "copie2", "ucomm", "rep_superF.1"),
    new = c("copy1", "copy2", "community", "superfamily")
)
httHits_filtered_toSave <- httHits_filtered_toSave[, c(
    "copy1", "copy2", "ID.1", "ID.2", "TEconsensus.1", "TEconsensus.2", "species.1", "species.2",
    "superfamily", "mrca", "divTime", "pID", "length", "qStart", "qEnd", "sStart", "sEnd", 
    "dN", "dS", "length.aa", "community", "hitGroup")]
fwrite(httHits_filtered_toSave, "HTThitsAssessed_perPairs_retained.txt", sep='\t')

meta <- fread("metadata.tbl")
meta$species <- gsub(" ", "_", meta$species)

##############################################

#Count number of HTT for each pair
dt <- select(httHits_filtered, c(species.1, species.2, hitGroup)) %>%
    #Make sure always the same species on the right or left
    mutate(., species.1b=ifelse(species.1<species.2, species.1, species.2), species.2=ifelse(species.1<species.2, species.2, species.1), species.1=species.1b) %>%
    group_by(species.1, species.2, hitGroup) %>% summarize(n=n()) %>% summarize(n=n())

#Save number of of HTT between pairs (does not contain pairs for which 0 transfers)
fwrite(dt, "nbHTTevents_perPair.txt", sep='\t') 

#Also make table that also contain pairs of 0 HTT

#For that, get all the pairs we investigated
pairs <- fread("Busco_new/pairsALL.txt")
mrca_FilterB <- readRDS("mrca_FilterB.Rds")
pairsInvestigated <- filter(pairs, !(divTime <=80 | mrca %in% mrca_FilterB)) 

#Make sure always the same species on the right or left
pairsInvestigated <- mutate(pairsInvestigated, species.1=ifelse(sp1<sp2, sp1, sp2), species.2=ifelse(sp1<sp2, sp2, sp1)) %>%
    select(., -c(sp1, sp2, assembly.1, assembly.2))

#add in dt all pairs investigated with 0 event
dt <- left_join(pairsInvestigated, dt, by=c("species.1", "species.2")) 

#replace NA by 0
dt[is.na(n),]$n = 0

#That one also  contain pairs for which 0 transfers
fwrite(dt, "nbHTTevents_perPairAll.txt", sep='\t')


############################################################
############# TEST whether more HTT in aquatic  ############
############################################################

if(nbTaxa==19){
    metaIncluded <- meta
} else if(nbTaxa==17){
    metaIncluded <- filter(meta, ! group %in%  c("Actinopterygii", "Coelacanthi", "Palaeoptera"))
} else if(nbTaxa==10){
    metaIncluded <- filter(meta,
        group %in% unique(filter(meta, stage_habitat=="full" & habitat=="aquatic")$group) & stage_habitat=="full" &
        ! group %in%  c("Actinopterygii", "Coelacanthi", "Palaeoptera")) 
} else if(nbTaxa==11){
  metaIncluded <- filter(meta,
                         group %in% unique(filter(meta, stage_habitat=="full" & habitat=="aquatic")$group) & stage_habitat=="full" &
                            group != "Palaeoptera")
}
#Only keep columns I need
metaIncluded <- metaIncluded[, c("species", "phylum", "group", "habitat")]

#Here, we tried several methods, each has its function

#Function to sample our dataset and see whether more HTT in aquatic or terrestrial
samplingTest <- function(simulationN=1, nbTaxa, theoretical=F, metric="mean"){

  print(paste("Starting simulation", simulationN))

  #We want 1000 samplings
  for(i in 1:1000){

      #Randomly sample 1 species per group/habitat
      meta_sub <- metaIncluded %>% group_by(group, habitat) %>% sample_n(1) %>% data.table()

      #If fish (only when nbTaxa = 19), we actualy just want one, not 1 actino + 1 coelacanthi
      if(nbTaxa==19){
        # Filrtly delete the 2 sampled fish above
        meta_sub <- meta_sub[! group %in% c("Actinopterygii", "Coelacanthi"),]
        # Then pick just one (can be actino or coelacanthi) 
        fish <- metaIncluded[group %in%  c("Actinopterygii", "Coelacanthi"),]  %>% sample_n(1) %>%
            #and rename its group
            mutate(., group="Euteleostomi") 
        #Put in the same table, the 18 * 2 species sampled above with the 1 fish
        meta_sub <- rbind(data.table(meta_sub), fish)
        #Get the 2 species name that were sampled among the nested groups (fish & palaeoptera)
        speciesNested <- meta_sub[group %in% c("Euteleostomi", "Palaeoptera"),]$species
      }

      #Subset from table of HTT per pairs only those that involve species we sampled
      dt_sub_i <-  dt[species.1 %in% meta_sub$species & species.2 %in% meta_sub$species,]

      #For each sampled species, count its number of HTT events when taking into account only those sampled
      dt_sub_i <- rbind(
            rename(dt_sub_i[, c("species.1", "n")], species="species.1"),
            rename(dt_sub_i[, c("species.2", "n")], species="species.2")
            ) %>%  group_by(species) %>% summarize(nTot = sum(n))%>%
            right_join(., meta_sub, by = "species") %>% mutate(.,nTot = ifelse(is.na(nTot), 0, nTot)) %>% data.table()

      #For each group, write how many HTT in the aquatic species and in the terrestrial one
      if(nbTaxa!=19){ #If nested taxon are not included, easy
       
        dt_sub_i <- left_join(
            rename(dt_sub_i[habitat=="aquatic", c("group", "nTot")], nA="nTot"),
            rename(dt_sub_i[habitat=="terrestrial", c("group", "nTot")], nT = "nTot"),
            by = "group"
            )
   
       } else { #For nested, we need a value for terrestrial, so more complex
       
        #For euteleostomi: A = fish, T = terrestrial chordata. So for T, a total of 6 species were sampled since 6 groups
        #--> 7 lines for euteleostomi
        dt_sub_i_euteleostomi <- dt_sub_i[phylum=="Chordata" & habitat=="terrestrial",] %>%  
            rbind(., dt_sub_i[group=="Euteleostomi",]) %>% mutate(., group="Euteleostomi")

        #For palaeoptera: A = palaeoptera, T = all other terretrial insects. So for T, a total of 7 species were sampled since 7 groups
        #--> 8 lines for palaeoptera
        dt_sub_i_insecta <- dt_sub_i[phylum=="Arthropoda" & habitat=="terrestrial" & (!group %in% c("Chelicerata", "Crustacea")),] %>% 
            rbind(., dt_sub_i[group=="Palaeoptera",]) %>% mutate(., group="Insecta")

        #Cat 3 tables togehter :
        dt_sub_i <- rbind(
            rbind(dt_sub_i[!species %in% speciesNested,], #Number of HTT in each species of the 17 taxa (17*2 species)
                dt_sub_i_euteleostomi
                ),
            dt_sub_i_insecta #Number of HTT in the sampled palaeoptera + the 7 terrestrial insects (also duplication with above)
            )  
        #For aquatic species, it is easy since one per group. Simply, select columns and lines of interest
        Ai <- rename(dt_sub_i[habitat=="aquatic",  c("group", "nTot")], nA = "nTot")
        #For terrrestial, as easy for the 17 easy group. But need to cacluate the average for the nested groups
        Ti <- dt_sub_i[habitat=="terrestrial", c("group", "nTot")] %>% group_by(group) %>% summarize(nT = mean(nTot)) #For terrestrial, some group have several species so I do the average

        #Write in same table HTT in A vs T 
        dt_sub_i <- left_join(Ai, Ti, by = "group")

    }

    #Calculate the difference between A-T
    dt_sub_i[, delta_t := nA-nT]

    #Permute A / T
    if(theoretical==T){
      #1/2 chances to shuffle tag A|T
      dt_sub_i[, shuffle := sample(c(-1,1), nbTaxa, replace=T)]      
      
      #if 1, we keep the same tag for that group (so delta is the same)
      #otherwise (if 2), we shuffle the tag, so delta change sign
      dt_sub_i[, delta_t_th := delta_t*shuffle]

      delta_t <- dt_sub_i$delta_t_th
      
    } else {
      delta_t <- dt_sub_i$delta_t
    }


    if(metric=="mean"){
      delta = mean(delta_t)
    } else if(metric=="median"){
      delta = median(delta_t)
    }
    
    #save all these values of delta 
    if(i==1){
      out = delta
    } else {
      out = rbind(out, delta)
    }
    
  }
  
  if(metric=="mean"){
    return(mean(out))
  } else if(metric=="median"){
    return(median(out))
  }
}

#Try change sign mean delta instead each delta_t
samplingTest2 <- function(simulationN=1, nbTaxa, theoretical=F, metric="mean"){

  print(paste("Starting simulation", simulationN))

  #We want 1000 samplings
  for(i in 1:1000){

      #Randomly sample 1 species per group/habitat
      meta_sub <- metaIncluded %>% group_by(group, habitat) %>% sample_n(1) %>% data.table()

      #If fish (only when nbTaxa = 19), we actualy just want one, not 1 actino + 1 coelacanthi
      if(nbTaxa==19){
        # Filrtly delete the 2 sampled fish above
        meta_sub <- meta_sub[! group %in% c("Actinopterygii", "Coelacanthi"),]
        # Then pick just one (can be actino or coelacanthi) 
        fish <- metaIncluded[group %in%  c("Actinopterygii", "Coelacanthi"),]  %>% sample_n(1) %>%
            #and rename its group
            mutate(., group="Euteleostomi") 
        #Put in the same table, the 18 * 2 species sampled above with the 1 fish
        meta_sub <- rbind(data.table(meta_sub), fish)
        #Get the 2 species name that were sampled among the nested groups (fish & palaeoptera)
        speciesNested <- meta_sub[group %in% c("Euteleostomi", "Palaeoptera"),]$species
      }

      #Subset from table of HTT per pairs only those that involve species we sampled
      dt_sub_i <-  dt[species.1 %in% meta_sub$species & species.2 %in% meta_sub$species,]

      #For each sampled species, count its number of HTT events when taking into account only those sampled
      dt_sub_i <- rbind(
            rename(dt_sub_i[, c("species.1", "n")], species="species.1"),
            rename(dt_sub_i[, c("species.2", "n")], species="species.2")
            ) %>%  group_by(species) %>% summarize(nTot = sum(n))%>%
            right_join(., meta_sub, by = "species") %>% mutate(.,nTot = ifelse(is.na(nTot), 0, nTot)) %>% data.table()

      #For each group, write how many HTT in the aquatic species and in the terrestrial one
      if(nbTaxa!=19){ #If nested taxon are not included, easy
       
        dt_sub_i <- left_join(
            rename(dt_sub_i[habitat=="aquatic", c("group", "nTot")], nA="nTot"),
            rename(dt_sub_i[habitat=="terrestrial", c("group", "nTot")], nT = "nTot"),
            by = "group"
            )
   
       } else { #For nested, we need a value for terrestrial, so more complex
       
        #For euteleostomi: A = fish, T = terrestrial chordata. So for T, a total of 6 species were sampled since 6 groups
        #--> 7 lines for euteleostomi
        dt_sub_i_euteleostomi <- dt_sub_i[phylum=="Chordata" & habitat=="terrestrial",] %>%  
            rbind(., dt_sub_i[group=="Euteleostomi",]) %>% mutate(., group="Euteleostomi")

        #For palaeoptera: A = palaeoptera, T = all other terretrial insects. So for T, a total of 7 species were sampled since 7 groups
        #--> 8 lines for palaeoptera
        dt_sub_i_insecta <- dt_sub_i[phylum=="Arthropoda" & habitat=="terrestrial" & (!group %in% c("Chelicerata", "Crustacea")),] %>% 
            rbind(., dt_sub_i[group=="Palaeoptera",]) %>% mutate(., group="Insecta")

        #Cat 3 tables togehter :
        dt_sub_i <- rbind(
            rbind(dt_sub_i[!species %in% speciesNested,], #Number of HTT in each species of the 17 taxa (17*2 species)
                dt_sub_i_euteleostomi
                ),
            dt_sub_i_insecta #Number of HTT in the sampled palaeoptera + the 7 terrestrial insects (also duplication with above)
            )  
        #For aquatic species, it is easy since one per group. Simply, select columns and lines of interest
        Ai <- rename(dt_sub_i[habitat=="aquatic",  c("group", "nTot")], nA = "nTot")
        #For terrrestial, as easy for the 17 easy group. But need to cacluate the average for the nested groups
        Ti <- dt_sub_i[habitat=="terrestrial", c("group", "nTot")] %>% group_by(group) %>% summarize(nT = mean(nTot)) #For terrestrial, some group have several species so I do the average

        #Write in same table HTT in A vs T 
        dt_sub_i <- left_join(Ai, Ti, by = "group")

    }

    #Calculate the difference between A-T
    dt_sub_i[, delta_t := nA-nT]


    if(metric=="mean"){
      delta = mean(dt_sub_i$delta_t)
    } else if(metric=="median"){
      delta = median(dt_sub_i$delta_t)
    }

    #Permute A / T
    if(theoretical==T){
      #1/2 chances to shuffle tag A|T
      shuffle := sample(c(-1,1), 1)    
      
      delta <- delta * suffle
    }
    
    #save all these values of delta 
    if(i==1){
      out = delta
    } else {
      out = rbind(out, delta)
    }
    
  }
  
  #Calcualte mean (or average) of the 1000 delta
  if(metric=="mean"){
    return(mean(out))
  } else if(metric=="median"){
    return(median(out))
  }
}

#Rosina's suggestion:
samplingTestAlt <- function(nbTaxa, sampling_n, sim_n){
  # sampling_n is the number of samplings
  # sim_n is the number of simulations per sampling
  
  #We want 1000 samplings
  for(i in 1:sampling_n){
    # print out every 100th number
    if(i%%100 == 0) {
      print(paste("Sampling", i, "out of", sampling_n))
    }
    
    #Randomly sample 1 species per group/habitat
    meta_sub <- metaIncluded %>% group_by(group, habitat) %>% sample_n(1) %>% data.table()
    
    #Subset from table of HTT per pairs only those that involve species we sampled
    dt_sub_i <-  dt[species.1 %in% meta_sub$species & species.2 %in% meta_sub$species,]
    
    #For each sampled species, count its number of HTT events when taking into account only those sampled
    dt_sub_i <- rbind(
      rename(dt_sub_i[, c("species.1", "n")], species="species.1"),
      rename(dt_sub_i[, c("species.2", "n")], species="species.2")
    ) %>%  group_by(species) %>% summarize(nTot = sum(n))%>%
      right_join(., meta_sub, by = "species") %>% mutate(.,nTot = ifelse(is.na(nTot), 0, nTot)) %>% data.table()
    
    #For each group, write how many HTT in the aquatic species and in the terrestrial one      
    dt_sub_i <- left_join(
      rename(dt_sub_i[habitat=="aquatic", c("group", "nTot")], nA="nTot"),
      rename(dt_sub_i[habitat=="terrestrial", c("group", "nTot")], nT = "nTot"),
      by = "group"
    )
    #Calculate the difference between A-T
    dt_sub_i[, delta_t := nA-nT]
    
    # median of the 17 deltas for the current sampling
    obs_delta <- median(dt_sub_i$delta_t)
    # calculate sim_n simulated deltas
    sim_deltas <- c()
    for (sim in 1:sim_n) {
      shuffle <- sample(c(-1,1), nbTaxa, replace=T)
      #if 1, we keep the same tag for that group (so delta is the same)
      #otherwise (if -1), we shuffle the tag, so delta change sign
      shuffled_deltas <- dt_sub_i[, delta_t] * shuffle
      sim_deltas <- c(sim_deltas, median(shuffled_deltas))
    }
    # calculate P
    P <- mean(sim_deltas <= obs_delta) * 2
    expected <- median(sim_deltas)
    if (P>1) {
      P <- 1
    }
    if(i==1){
      out_obs <- obs_delta
      Ps <- P
      out_exp <- expected
    } else {
      out_obs <- rbind(out_obs, obs_delta)
      Ps <- rbind(Ps, P)
      out_exp <- rbind(out_exp, expected)
    }
  }
  
  print("Number significant:")
  print(sum(alt_results$Pvalues < 0.05))
  print("Median effect across samplings:")
  print(median(alt_results$observed))
  
  return(list("observed" = out_obs, "Pvalues" = Ps, "expected" = out_exp))
}

#Rosina's suggestion but compatible with nbTaxa=19
samplingTestAlt_HM <- function(nbTaxa, sampling_n, sim_n){
  # sampling_n is the number of samplings
  # sim_n is the number of simulations per sampling
  
  #We want 1000 samplings
  for(i in 1:sampling_n){
    # print out every 100th number
    if(i%%100 == 0) {
      print(paste("Sampling", i, "out of", sampling_n))
    }
    
    #Randomly sample 1 species per group/habitat
    meta_sub <- metaIncluded %>% group_by(group, habitat) %>% sample_n(1) %>% data.table()
    
    if(nbTaxa==19){
      # Filrtly delete the 2 sampled fish above
      meta_sub <- meta_sub[! group %in% c("Actinopterygii", "Coelacanthi"),]
      # Then pick just one (can be actino or coelacanthi) 
      fish <- metaIncluded[group %in%  c("Actinopterygii", "Coelacanthi"),]  %>% sample_n(1) %>%
        #and rename its group
        mutate(., group="Euteleostomi") 
      #Put in the same table, the 18 * 2 species sampled above with the 1 fish
      meta_sub <- rbind(data.table(meta_sub), fish)
      #Get the 2 species name that were sampled among the nested groups (fish & palaeoptera)
      speciesNested <- meta_sub[group %in% c("Euteleostomi", "Palaeoptera"),]$species
    }    

    #Subset from table of HTT per pairs only those that involve species we sampled
    dt_sub_i <-  dt[species.1 %in% meta_sub$species & species.2 %in% meta_sub$species,]
    
    #For each sampled species, count its number of HTT events when taking into account only those sampled
    dt_sub_i <- rbind(
      rename(dt_sub_i[, c("species.1", "n")], species="species.1"),
      rename(dt_sub_i[, c("species.2", "n")], species="species.2")
    ) %>%  group_by(species) %>% summarize(nTot = sum(n))%>%
      right_join(., meta_sub, by = "species") %>% mutate(.,nTot = ifelse(is.na(nTot), 0, nTot)) %>% data.table()
    
    #For each group, write how many HTT in the aquatic species and in the terrestrial one      

    if(nbTaxa!=19){ #If nested taxon are not included, easy
      
      dt_sub_i <- left_join(
        rename(dt_sub_i[habitat=="aquatic", c("group", "nTot")], nA="nTot"),
        rename(dt_sub_i[habitat=="terrestrial", c("group", "nTot")], nT = "nTot"),
        by = "group"
      )
      
    } else { #For nested, we need a value for terrestrial, so more complex
      
      #For euteleostomi: A = fish, T = terrestrial chordata. So for T, a total of 6 species were sampled since 6 groups
      #--> 7 lines for euteleostomi
      dt_sub_i_euteleostomi <- dt_sub_i[phylum=="Chordata" & habitat=="terrestrial",] %>%  
        rbind(., dt_sub_i[group=="Euteleostomi",]) %>% mutate(., group="Euteleostomi")
      
      #For palaeoptera: A = palaeoptera, T = all other terretrial insects. So for T, a total of 7 species were sampled since 7 groups
      #--> 8 lines for palaeoptera
      dt_sub_i_insecta <- dt_sub_i[phylum=="Arthropoda" & habitat=="terrestrial" & (!group %in% c("Chelicerata", "Crustacea")),] %>% 
        rbind(., dt_sub_i[group=="Palaeoptera",]) %>% mutate(., group="Insecta")
      
      #Cat 3 tables togehter :
      dt_sub_i <- rbind(
        rbind(dt_sub_i[!species %in% speciesNested,], #Number of HTT in each species of the 17 taxa (17*2 species)
              dt_sub_i_euteleostomi
        ),
        dt_sub_i_insecta #Number of HTT in the sampled palaeoptera + the 7 terrestrial insects (also duplication with above)
      )  
      #For aquatic species, it is easy since one per group. Simply, select columns and lines of interest
      Ai <- rename(dt_sub_i[habitat=="aquatic",  c("group", "nTot")], nA = "nTot")
      #For terrrestial, as easy for the 17 easy group. But need to cacluate the average for the nested groups
      Ti <- dt_sub_i[habitat=="terrestrial", c("group", "nTot")] %>% group_by(group) %>% summarize(nT = mean(nTot)) #For terrestrial, some group have several species so I do the average
      
      #Write in same table HTT in A vs T 
      dt_sub_i <- left_join(Ai, Ti, by = "group")
      
    }
    
    
    #Calculate the difference between A-T
    dt_sub_i[, delta_t := nA-nT]
    
    # median of the 17 deltas for the current sampling
    obs_delta <- median(dt_sub_i$delta_t)
    # calculate sim_n simulated deltas
    sim_deltas <- c()
    for (sim in 1:sim_n) {
      shuffle <- sample(c(-1,1), nbTaxa, replace=T)
      #if 1, we keep the same tag for that group (so delta is the same)
      #otherwise (if -1), we shuffle the tag, so delta change sign
      shuffled_deltas <- dt_sub_i[, delta_t] * shuffle
      sim_deltas <- c(sim_deltas, median(shuffled_deltas))
    }
    # calculate P
    P <- mean(sim_deltas >= obs_delta) #One tailed, and check only whether excess aquatic
    expected <- median(sim_deltas)
    if (P>1) {
      P <- 1
    }
    if(i==1){
      out_obs <- obs_delta
      Ps <- P
      out_exp <- expected
    } else {
      out_obs <- rbind(out_obs, obs_delta)
      Ps <- rbind(Ps, P)
      out_exp <- rbind(out_exp, expected)
    }
  }
  
  print("Number significant:")
  print(sum(alt_results$Pvalues < 0.05))
  print("Median effect across samplings:")
  print(median(alt_results$observed))
  
  return(list("observed" = out_obs, "Pvalues" = Ps, "expected" = out_exp))
}

#ALso allow nbTaxa=11
samplingTestAlt_HM2 <- function(nbTaxa, sampling_n, sim_n){
  # sampling_n is the number of samplings
  # sim_n is the number of simulations per sampling
  
  #We want 1000 samplings
  for(i in 1:sampling_n){
    # print out every 100th number
    if(i%%100 == 0) {
      print(paste("Sampling", i, "out of", sampling_n))
    }
    
    #Randomly sample 1 species per group/habitat
    meta_sub <- metaIncluded %>% group_by(group, habitat) %>% sample_n(1) %>% data.table()
    
    if("Actinopterygii" %chin% meta_sub$group){
      # Firstly delete the 2 sampled fish above
      meta_sub <- meta_sub[! group %in% c("Actinopterygii", "Coelacanthi"),]
      # Then pick just one (can be actino or coelacanthi) 
      fish <- metaIncluded[group %in%  c("Actinopterygii", "Coelacanthi"),]  %>% sample_n(1) %>%
        #and rename its group
        mutate(., group="Euteleostomi") 
      #Put in the same table, the 18 * 2 species sampled above with the 1 fish
      meta_sub <- rbind(data.table(meta_sub), fish)
      #Get the 2 species name that were sampled among the nested groups (fish & palaeoptera)
      speciesNested <- meta_sub[group %in% c("Euteleostomi", "Palaeoptera"),]$species
    }    
    
    #Subset from table of HTT per pairs only those that involve species we sampled
    dt_sub_i <-  dt[species.1 %in% meta_sub$species & species.2 %in% meta_sub$species,]
    
    #For each sampled species, count its number of HTT events when taking into account only those sampled
    dt_sub_i <- rbind(
      rename(dt_sub_i[, c("species.1", "n")], species="species.1"),
      rename(dt_sub_i[, c("species.2", "n")], species="species.2")
    ) %>%  group_by(species) %>% summarize(nTot = sum(n))%>%
      right_join(., meta_sub, by = "species") %>% mutate(.,nTot = ifelse(is.na(nTot), 0, nTot)) %>% data.table()
    
    #For each group, write how many HTT in the aquatic species and in the terrestrial one      
    
    if(nbTaxa!=19 & nbTaxa!=11){ #If nested taxon are not included, easy
      
      dt_sub_i <- left_join(
        rename(dt_sub_i[habitat=="aquatic", c("group", "nTot")], nA="nTot"),
        rename(dt_sub_i[habitat=="terrestrial", c("group", "nTot")], nT = "nTot"),
        by = "group"
      )
      
    } else { #For nested, we need a value for terrestrial, so more complex
      
      #For euteleostomi: A = fish, T = terrestrial chordata. So for T, a total of 6 species were sampled since 6 groups
      #--> 7 lines for euteleostomi
      dt_sub_i_euteleostomi <- dt_sub_i[phylum=="Chordata" & habitat=="terrestrial",] %>%  
        rbind(., dt_sub_i[group=="Euteleostomi",]) %>% mutate(., group="Euteleostomi")
      
      if(nbTaxa==19){
      
        #For palaeoptera: A = palaeoptera, T = all other terretrial insects. So for T, a total of 7 species were sampled since 7 groups
        #--> 8 lines for palaeoptera
        dt_sub_i_insecta <- dt_sub_i[phylum=="Arthropoda" & habitat=="terrestrial" & (!group %in% c("Chelicerata", "Crustacea")),] %>% 
        rbind(., dt_sub_i[group=="Palaeoptera",]) %>% mutate(., group="Insecta")
      
        #Cat 3 tables togehter :
        dt_sub_i <- rbind(
          rbind(dt_sub_i[!species %in% speciesNested,], #Number of HTT in each species of the 17 taxa (17*2 species)
                dt_sub_i_euteleostomi
          ),
          dt_sub_i_insecta #Number of HTT in the sampled palaeoptera + the 7 terrestrial insects (also duplication with above)
        )
      } else {
        #Cat 2 tables togehter :
        dt_sub_i <- rbind(
          dt_sub_i[!species %in% speciesNested,], #Number of HTT in each species of the 17 taxa (17*2 species)
          dt_sub_i_euteleostomi
          )
      }
      
      #For aquatic species, it is easy since one per group. Simply, select columns and lines of interest
      Ai <- rename(dt_sub_i[habitat=="aquatic",  c("group", "nTot")], nA = "nTot")
      #For terrrestial, as easy for the 17 easy group. But need to cacluate the average for the nested groups
      Ti <- dt_sub_i[habitat=="terrestrial", c("group", "nTot")] %>% group_by(group) %>% summarize(nT = median(nTot)) #For terrestrial, some group have several species so I do the average
      
      #Write in same table HTT in A vs T 
      dt_sub_i <- left_join(Ai, Ti, by = "group")
      
    }
    
    
    #Calculate the difference between A-T
    dt_sub_i[, delta_t := nA-nT]
    
    # median of the 17 deltas for the current sampling
    obs_delta <- median(dt_sub_i$delta_t)
    # calculate sim_n simulated deltas
    sim_deltas <- c()
    for (sim in 1:sim_n) {
      shuffle <- sample(c(-1,1), nbTaxa, replace=T)
      #if 1, we keep the same tag for that group (so delta is the same)
      #otherwise (if -1), we shuffle the tag, so delta change sign
      shuffled_deltas <- dt_sub_i[, delta_t] * shuffle
      sim_deltas <- c(sim_deltas, median(shuffled_deltas))
    }
    # calculate P
    P <- mean(sim_deltas >= obs_delta) #One tailed, and check only whether excess aquatic
    expected <- median(sim_deltas)
    if (P>1) {
      P <- 1
    }
    if(i==1){
      out_obs <- obs_delta
      Ps <- P
      out_exp <- expected
    } else {
      out_obs <- rbind(out_obs, obs_delta)
      Ps <- rbind(Ps, P)
      out_exp <- rbind(out_exp, expected)
    }
  }
  
  print("Number significant:")
  print(sum(alt_results$Pvalues < 0.05))
  print("Median effect across samplings:")
  print(median(alt_results$observed))
  
  return(list("observed" = out_obs, "Pvalues" = Ps, "expected" = out_exp))
}

##########################################

#Run my function:
obs_mean_17 <- samplingTest(1, 17)
th_mean_17 <- rbind(mclapply(1:1000, samplingTest, nbTaxa=17, theoretical=T, mc.cores = 6))
th_mean_17 <- unlist(th_mean_17) 

th <- th_mean_17
obs <- obs_mean_17

limMinus <- min(c(th, obs))
limMax <- max(c(th, obs))
lim <- round(max(limMinus, limMax)+1)

hist(th, col = "gray",border = F, breaks = 50, xlab = "mean delta", xlim = c(-lim,lim))
abline(v=obs, col = "red")


##########################################

#Run Rosina function
alt_results <- samplingTestAlt(17, 1000, 1000)
save(alt_results, file="simulationsAT/samplingTestAlt_17taxa.RData")
save(alt_results, file="simulationsAT/samplingTestAlt_17taxa_1tailed.RData")
#Read with load()

alt_results <- samplingTestAlt(10, 1000, 1000)
save(alt_results, file="simulationsAT/samplingTestAlt_10taxa.RData")
save(alt_results, file="simulationsAT/samplingTestAlt_10taxa_1tailed.RData")


alt_results <- samplingTestAlt_HM(19, 1000, 1000)
save(alt_results, file="simulationsAT/samplingTestAlt_19taxa.RData")
save(alt_results, file="simulationsAT/samplingTestAlt_19taxa_1tailed.RData")

alt_results <- samplingTestAlt_HM2(11, 1000, 1000)
save(alt_results, file="simulationsAT/samplingTestAlt_11taxa.RData")
save(alt_results, file="simulationsAT/samplingTestAlt_11taxa_1tailed.RData")


#I prefer to have the plot outside the function for modification afterward
excess_over_expected <- alt_results$observed - alt_results$expected
par(mfrow = c(1,2))
hist(excess_over_expected, xlab = "Observed - expected delta",
     main = "", col = "red3", border = NA, breaks = 20)
perc_sig <- mean(alt_results$Pvalues < 0.05) * 100
print("Number significant:")
print(sum(alt_results$Pvalues < 0.05))
print("Median effect across samplings:")
print(median(alt_results$observed))
hist(alt_results$Pvalues,
     xlim = c(0,1),
     xlab = "P-value",
     col = "steelblue3", border = NA, breaks = 20,
     main = paste("% significant:", signif(perc_sig, 2)))
abline(v = 0.05, col = "black", lty = 3)

#OR just plot observed and no stats
hist(alt_results$observed, breaks = 20, main = "")
#Save as samplingTestAlt_19taxa_justObserved

#Fr sup:
par(mfrow = c(1,2))
hist(alt_results$observed, xlab = "∆",
     main = "", col = "red3", border = NA, breaks = 20)
perc_sig <- mean(alt_results$Pvalues < 0.05) * 100
print("Number significant:")
print(sum(alt_results$Pvalues < 0.05))
print("Median effect across samplings:")
print(median(alt_results$observed))
hist(alt_results$Pvalues,
     xlim = c(0,1),
     xlab = "P-value",
     col = "steelblue3", border = NA, breaks = 20,
     main = paste("% significant:", signif(perc_sig, 2)))
abline(v = 0.05, col = "black", lty = 3)
