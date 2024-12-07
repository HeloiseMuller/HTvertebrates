library(data.table)
library(ggplot2)
library(dplyr)

##Metadata
meta <- fread("metadata.tbl")
meta$species <- gsub(" ", "_", meta$species)

##All pairs among which we looked for HTT and their number of HT
transfers <- fread("nbHTTevents_perPairAll")

##Table with all possible pairs of species of the dataset, even those for which we did not look for HTT:
	#(247*247-247)/2 + 247 rows
pairs <- fread("Busco_new/pairsALL.txt")

##Put in a same taxon species among which we didn't look for HTT
#the clade too related  :
mrca_FilterB <- readRDS("mrca_FilterB.Rds") #dS Busco<0.85
mrca_young <- unique(pairs[divTime<80 & divTime!=0,]$mrca) #clades diverged <40 Myrs ago

#The species to group in a same taxa are those ones:
sp_toGrp <- pairs[mrca %in% c(mrca_FilterB, mrca_young), c("sp1", "mrca", "divTime")] %>% rename(., sp2="sp1") %>% rbind(., pairs[mrca %in% c(mrca_FilterB, mrca_young), c("sp2", "mrca", "divTime")])  %>% rename(., sp="sp2")
#for each species, select the oldest clade among which we didn't look for HTT
setorder(sp_toGrp, -divTime)
sp_toGrp <- sp_toGrp[!duplicated(sp_toGrp$sp), -"divTime"] %>% rename(., spClade="mrca")

#If need to make a tree with one species per clade, I have to save somewhere a representative sp for each clade (here representative just mean that it is its name we will keep in the tree)
rep <- data.table(spRep = meta[!species %in% sp_toGrp[duplicated(spClade),]$sp,]$species)
rep <- left_join(rep, sp_toGrp, by = c("spRep"="sp"))
rep <- mutate(rep, spClade = ifelse(is.na(spClade), spRep, spClade))
fwrite(rep, "spRep_clade.txt", sep='\t')

##For each species, give it a clade: itself if no related species, or the name of the common node with related species otherwise
sp_clades <- data.frame(sp = meta[!species %in% sp_toGrp$sp,]$species, spClade = meta[!species %in% sp_toGrp$sp,]$species) %>% rbind(., sp_toGrp)

##Go back to data about HTT, and add info about clade
transfers <- left_join(transfers, sp_clades, by=c("species.1"= "sp")) %>% left_join(., sp_clades, by=c("species.2"="sp"), suffix=c(".1", ".2"))

transfers[, pairClade := ifelse(spClade.1<spClade.2, paste0(spClade.1, "-", spClade.2),  paste0(spClade.2, "-", spClade.1))]

fwrite(transfers, "nbHTTevents_perPairAll_withClades.txt", sep='\t')

#Keep only info we need
transfers2 <- transfers[, c("divTime", "pairClade", "n")]

##Sample a unique pair of species per pair of clade
sampling1 <- transfers2 %>% group_by(pairClade) %>% sample_n(1) %>% data.table

linear_model <- lm(sampling1$n ~ sampling1$divTime + sampling1$same_habitat)


#####
#Calculate median HTT for pairs of species of similar divTime
#For this, make window of divTime

#I make window of 100myrs, overlapping of 50yrs
sampling1[, bin:=cut_width(divTime, width=100, boundaries=0)]
sampling1[, bin2:=cut_width(divTime, width=100, boundaries=0, center=min(divTime)/2)]

#Because most pair of species have 0 HTT, most median would be 0
#I cannot use mean either because of some extrem values
# --> Separate info : 
  #1) a plot that contains portion of pair of species with 0 HTT (p1)
  #2) a plot that shows the median number if HTT when keeping only pairs of species involed in >0 HTT (p2)

#Start with p1
zero <- sampling1[n==0,]
zero_a <- zero %>% group_by(bin) %>% summarize(nbPairs_at0=n(), medianTime = median(divTime)) %>% data.table
zero_b <- zero %>% group_by(bin2) %>% summarize(nbPairs_at0=n(), medianTime = median(divTime)) %>% data.table %>% rename(., bin="bin2")

#how many pairs per bin?
nbPairs_a <-  sampling1 %>% group_by(bin) %>% summarize(nbPairs = n()) %>% data.table
nbPairs_b <-  sampling1 %>% group_by(bin2) %>% summarize(nbPairs = n()) %>% data.table %>% rename(., bin="bin2")

zero_ratio <- left_join(rbind(zero_a,zero_b), rbind(nbPairs_a, nbPairs_b)) %>% mutate(., ratioPairs_at0 = nbPairs_at0/nbPairs)

p1 <- ggplot(zero_ratio, aes(x=medianTime, y=ratioPairs_at0, size=nbPairs)) + geom_point(color="red") + theme_bw() + ylim(0,1) +
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(1, 5300), 
                        breaks =c(1,seq(10,100, 10), 500, seq(1500,5500, 1000)))

#Then plot p2
sampling1_noExtr <- sampling1[n>0,] 
setorder(sampling1_noExtr, divTime)

a <- sampling1_noExtr %>% group_by(bin) %>% summarize(meanHTT=mean(n), medianHTT=median(n), meanTime=mean(divTime), medianTime = median(divTime), nbPairs = n()) %>% data.table
b <- sampling1_noExtr %>% group_by(bin2) %>% summarize(meanHTT=mean(n), medianHTT=median(n), meanTime=mean(divTime), medianTime = median(divTime), nbPairs = n()) %>% data.table %>% rename(., bin="bin2")

p2 <- ggplot(rbind(a,b), aes(x=medianTime, y=medianHTT, size=nbPairs)) + geom_point(color="blue") + theme_bw() +
    scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(5, 800), #give min and max values nbPairs
                        breaks = c(5, 25,50, 100, 200, 300, 400, 500, 600, 700, 800))
                        
  
#########                      
#We can also plot both plots together

max(rbind(a,b)$medianHTT) #11.5 so scale 2nd axis at 12
scaleAxe = 12

pdf("Figure3a.pdf", width=9, height=7)
ggplot() +
    geom_point(data=rbind(a,b), aes(x=medianTime, y=medianHTT, size=nbPairs), col="darkolivegreen3")  +
    geom_point(data=zero_ratio, aes(x=medianTime, y=ratioPairs_at0*scaleAxe, size=nbPairs), col="magenta") +
    scale_y_continuous(name = "Median HTT above 0", sec.axis = sec_axis(~./scaleAxe, name="Proportion of pairs of species with 0 HTT") ) +
    theme_bw() +
    xlab("Median time of divergence in the bin") +
    theme(
    axis.title.y = element_text(color = "darkolivegreen3", size=13),
    axis.title.y.right = element_text(color = "magenta", size=13)
  ) +
    scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(1, 5300), 
                        breaks =c(1,10, 50, 100, 500, seq(1500,5500, 1000))) +
    #I want the legend in black but this does not work..
     guides(size = guide_legend(fill="black"))
dev.off()

#Calculate correlation with weight since number of pair of species very different in eac data point 
weighted_corr <- cov.wt(zero_ratio[, c("ratioPairs_at0","medianTime")], wt = log(zero_ratio$nbPairs), cor = TRUE) #log because max nbPairs is very high