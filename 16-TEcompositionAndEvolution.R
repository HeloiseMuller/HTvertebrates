## %######################################################%##
#                                                          #
####          This stage makes figures showing          ####
####       TE composition in species, hit groups,       ####
####           and HTT events, and           ####
####         analyses the dN/dS ratios of TEs.          ####
#                                                          #
## %######################################################%##


# This can be run at any point after stage 13

require(RColorBrewer)
require(data.table)
library(dplyr)
library(stringi)
library(stringr)

path = "~/Project/"
setwd(path) 

source("HTvFunctions.R")

# data file provided with the paper, which is a table of hits representing HTT
retainedHits <- fread("HTThitsAssessed_perClade_retained.txt")

# STEP ONE, we collect stats on TE composition in species, to build the barplots (figure 4) ---------------------------------------------------

# we import reports on numbers  of copies of at least 300 bp for the different genomes and 
# super families, generated in stage 2
compo <- fread("RepeatMasker/TEcomposition/supplementary-data-TEcomposition_per_species.txt")

#Group TE of same superfamilies under same name
compo = mutate(compo, superF = ifelse(str_detect(RepeatModeler_superfamily, "TcMar") | RepeatModeler_superfamily=="DNA/Mariner/Tc1", "DNA/Mariner",
    ifelse(str_detect(RepeatModeler_superfamily, "hAT") | str_detect(RepeatModeler_superfamily, "DNA/haT"), "DNA/hAT",
    ifelse(str_detect(RepeatModeler_superfamily, "Helitron"), "DNA/Helitron",
    ifelse(str_detect(RepeatModeler_superfamily, "Sola"), "DNA/Sola",
    ifelse(RepeatModeler_superfamily== "DNA/piggyBac" | RepeatModeler_superfamily=="DNA/PiggyBac" , "DNA/PiggyBac",
    ifelse(str_detect(RepeatModeler_superfamily, "DNA/MuLE") | str_detect(RepeatModeler_superfamily, "DNA/MULE"), "DNA/MuLE",
    ifelse(str_detect(RepeatModeler_superfamily, "DNA/Kolobok"), "DNA/Kolobok",
    ifelse(str_detect(RepeatModeler_superfamily, "DNA/Academ"), "DNA/Academ",
    ifelse(str_detect(RepeatModeler_superfamily, "CMC") | RepeatModeler_superfamily %in% c("DNA/EnSpm/CACTA", "DNA/Transib", "DNA/Chapaev"), "DNA/CMC",
    ifelse(str_detect(RepeatModeler_superfamily, "Crypton"), "DNA/Crypton",
    ifelse(str_detect(RepeatModeler_superfamily, "PIF") | RepeatModeler_superfamily %in% c("DNA/Harbinger", "DNA/ISL2EU"), "DNA/PIF-Harbinger",
    ifelse(str_detect(RepeatModeler_superfamily, "Ginger"), "DNA/Ginger",
    ifelse(str_detect(RepeatModeler_superfamily, "LTR/BEL") | str_detect(RepeatModeler_superfamily, "LTR/Pao"), "LTR/BEL-Pao",
    ifelse(str_detect(RepeatModeler_superfamily, "LTR/ERV"), "LTR/ERV",
    ifelse(str_detect(RepeatModeler_superfamily, "LTR/Gypsy"), "LTR/Gypsy",
    ifelse(str_detect(RepeatModeler_superfamily, "LINE/RTE") | str_detect(RepeatModeler_superfamily, "Proto2"), "LINE/RTE",
    ifelse(str_detect(RepeatModeler_superfamily, "LINE/L1") | RepeatModeler_superfamily=="LINE/Tx1", "LINE/L1",
    ifelse(str_detect(RepeatModeler_superfamily, "LINE/R2") | str_detect(RepeatModeler_superfamily, "LINE/CRE") | str_detect(RepeatModeler_superfamily, "R4") | RepeatModeler_superfamily %in% c("LINE/NeSL", "LINE/Hero"), "LINE/R2",
    ifelse(str_detect(RepeatModeler_superfamily, "LINE/L2") | str_detect(RepeatModeler_superfamily, "LINE/CR1") | RepeatModeler_superfamily %in% c("LINE/Daphne", "LINE/Crack", "LINE/Rex1"), "LINE/Jockey",
    ifelse(str_detect(RepeatModeler_superfamily, "LINE/R1") | RepeatModeler_superfamily %in% c("LINE/I-Nimb", "LINE/Loa" , "LINE/LOA", "LINE/Nimb", "LINE/Tad1", "LINE/Ingi"), "LINE/I",
    ifelse(RepeatModeler_superfamily %in% c("LTR/Caulimovirus", "LTR/Lentivirus"), "LTR/Retrovirus",
    #NA #replace by NA at the begining to check what left
    RepeatModeler_superfamily
))))))))))))))))))))))

#Check what has not been replaced
filter(compo, is.na(superF)) %>% select(., RepeatModeler_superfamily) %>% distinct()

#Make table of correspondance 
corres = select(compo, RepeatModeler_superfamily, superF) %>% distinct

#Use it for the table of hits
retainedHits[, "superF" := corres[chmatch(retainedHits$superfamily,RepeatModeler_superfamily), superF]]

compo[, class := ifelse(str_detect(superF, "DNA/"), "DNA", "RNA")]


# ------------------------------------------


# to make the baplots, we sum counts per super family
perSuperF <- compo[, .(
    nCopies = sum(number_of_copies_300bp),
    nCons = sum(number_of_consensus), 
    bp = sum(as.numeric(total_nucleotides))
), by = .(superF,  class)]

# we count hit groups and independent htt event per superfamily
HTTperSuperF <- retainedHits[, .(
  nTr = length(unique(hitGroup[independent])), 
  nHitGroup = length(unique(hitGroup))
  ), by = superF]

# the merge below removes super families not involved in HTT since we don't set all = T
perSuperF <- merge(perSuperF, HTTperSuperF, by = "superF", all = F) 

# we pool super families that are involved in only few HTT, for each TE class
# we define the number of independent transfer below which we pool super families within TE classes
lower <- 50L #Modify this number so less superfamily to plot than possible number of colors with brewer.pal()

# we sum numbers from these pooled superfamilies in a new table
smallFam <- perSuperF[nTr < lower, .(
    nCopies = sum(nCopies, na.rm = T),
    nCons = sum(nCons),
    bp = sum(bp),
    nTr = sum(nTr),
    nHitGroup = sum(nHitGroup)
), by = class]

# we create a table that has the small super families pooled per class
forBarplots <- rbind(
    data.table(
        superF = paste("Other", smallFam$class),
        smallFam
    ),
    perSuperF[nTr >= lower]
) 

# we sort by class then number of independent transfers, for plotting
forBarplots <- forBarplots[order(-class, !grepl("Other", superF), nTr)]

# we set colours. Barplot sectors of DNA transposons have darker colors
palDNA <- c(grey(0.5), brewer.pal(nrow(forBarplots[class=="DNA",])-1, "Paired"))
palRNA <- c(grey(0.7), brewer.pal(nrow(forBarplots[class=="RNA",])-1, "Pastel1"))

# we add a new column that tells which color a TE superfamily should have on the plot
forBarplots[, col := c(palRNA, palDNA)]
# we keep this table around as we do not draw the plots yet

#Also add normalization
#This allow to correct for the fact that some TE are better annotated than other
forBarplots[, nTr_perMillionCopies := nTr/(nCopies/10^6)]
forBarplots[, nTr_perThousandConsensus := nTr/(nCons/1000)]

# STEP TWO, we process results on Ka/Ks of TEs, to show on the same figure as the barplots -------------------------------------------------

# we import results from ka ks value between TEs within communities generated by stage 13-TEKaKsWinthinGenomes.R. 
alldNdS <- fread("TEs/TEevolution/TEdNdSndDistance.txt")

# we remove values that relate to hit communities in hitgroups that we did not retain
alldNdS <- alldNdS[com %in% retainedHits$community]

# we similar super families together
alldNdS[, superfamily:= .(corres[chmatch(alldNdS$superF, RepeatModeler_superfamily), superF])]

# we add a column denoting the TE class 
alldNdS[, class := ifelse(str_detect(superfamily, "DNA/"), "DNA", "RNA")]

# we again pool super families involved in less than 10 (or whatever number chose) independent transfers
alldNdS[superfamily %chin% perSuperF[nTr < lower, superF], 
        superfamily := stri_c("Other ", class)]


# we compute density curves of dN/dS ratios -----------------------------------

# the function below returns xy coordinates to drawn density curves 
# that all reach the same height in the Y axis.
# The difference in Y coordinates between super families on the plot will be 1 unit, 
# the x and y coordinates of the density curves (polygons) are placed in a single
# large data table (one row per vertex). For these curves, we use individual
# ratio rates (not averages per hit group)
# Instead of dN/dS, use dN-dS/dN+dS, which is borned by -1 and +1
densities <- alldNdS[dS > 0L & dN>0L & dN<9L & dS<9L, #Remove aberant values
  dens((dN-dS) / (dN+dS), length), 
  by = .(superfamily, htt, class)]

# we order the curves to match the future barplots
densities <- densities[order(chmatch(superfamily, forBarplots$superF), htt)]

# we compute the Y position of each super family on the density plots (+1 unit per super family)
densities[, Yoffset := toInteger(superfamily)]

# we add a column denoting the color of ka/ks curves, to match that of barplots. 
# Colors get replicated for each vertex, but we can afford that
densities[, col := forBarplots[chmatch(densities$superfamily, superF), col]]

# the curves showing rate of TEs involved in htt will be drawn negatively
densities[htt == T, y := -y]


# figure 4 of the paper
pdf("Figure_EvolTE.pdf", width = 7, height = 4.5)

# there will be six plots on the figure
par(
    lwd = 0.25,
    bty = "n",
    mfrow = c(1, 6), 
    las = 1,
    mai = c(0.6, 0.4, 0.2, 0.2)
)

layout(matrix(rep(1:6, each = 2), nrow = 2), widths = c(1.3, 1.3, 1.7, 1.3, 1.3, 4.1))

# border colours will be grey
border <- "grey"

#we draw the first barplot for TE composition in genomes
forBarplots[, barplot(
    cbind(nCopies / 10^6),
    names.arg = "Copies \n(millions)",
    xlim = 0:1,
    col = col,
    border = border,
    cex.axis = 0.7,
    cex.names = 0.7
)]

# the second barplot for number of consensus
 forBarplots[, barplot(
    cbind(nCons / 1000),
    names.arg = "Consensus \n(thousand)",
    xlim = 0:1,
    col = col,
    border = border,
    cex.axis = 0.7,
    cex.names = 0.7
)]

# the third barplot for htt counts
forBarplots[, barplot(
    cbind(nTr),
    names.arg = "Transfers",
    xlim = 0:1,
    width = 0.6,
    col = col,
    border = border,
    cex.axis = 0.7,
    cex.names = 0.7
)]

#Add 2 plots for normalization
forBarplots[, barplot(
    cbind(nTr_perMillionCopies),
    names.arg = "Transfers per \n million copies",
    xlim = 0:1,
    col = col,
    border = border,
    cex.axis = 0.7,
    cex.names = 0.7
)]

forBarplots[, barplot(
    cbind(nTr_perThousandConsensus),
    names.arg = "Transfers per \n thousand consensus",
    xlim = 0:1,
    col = col,
    border = border,
    cex.axis = 0.7,
    cex.names = 0.7
)]

# we plot Ka/Ks density polygons
# we prepare the empty plot on which polygons will be drawn
plot(
    x = c(-1.1, 1.1), #this one if use the other ratio
    y = c(0.5, max(densities$Yoffset + 0.5)),
    type = "n",
    yaxt = "n",
    ylab = "",
    xlab = "dN-dS/dN+dS", 
    bty = "n"
)

# because we draw the vertical line at ka/ks = 0 behind the polygons
abline(v = 0, lwd = 0.25, col = "grey") 

# we now draw the density polygons
# we split density curves per super family
# calling polygon() from the whole data table does not draw the correct
# density curve for certain super families, sometimes (bug with data.table?)
splitDens = split(densities, f = densities$superfamily)
p = lapply(splitDens, FUN = function(vertices) vertices[,polygon(
    x = x,
    y = y + Yoffset,
    col = col[1],
    border = border,
    lwd = 0.25
    )]
)

# we add the Y axis of superfamily names
p <- densities[!duplicated(superfamily), axis(2, at = Yoffset, labels = superfamily, las = 1, cex.axis=0.8)]

dev.off()