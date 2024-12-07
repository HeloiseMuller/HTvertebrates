## %######################################################%##
#                                                          #
####           This stages makes the convave            ####
####               tree showing HTT events.             ####
#                                                          #
## %######################################################%##

# This can be run at any point after stage 12

require(RColorBrewer)
require(ape)
require(dplyr)
require(data.table)
library(stringr)
library(stringi)
library(bit64)

path = "~/Project/"
setwd(path) 

source("HTvFunctions.R")
source("circularPlots.R")

#choose weither HTT events should be colors by TE types or habitat or just black
colEvents = "TEtype"
#colEvents = "habitat"
#colEvents = "black"

#Function from: http://blog.phytools.org/2012/01/function-to-get-descendant-node-numbers.html
getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}

########################""

# this script uses:
# - the data file provided with the paper: a table of hits representing HTT
retainedHits <- fread("HTThitsAssessed_perClade_retained.txt") #file save at the end script 12

# - the timetree
tree <-  read.tree("datasetTree.nwk")
meta = fread("metadata.tbl")


meta = mutate(meta, assembly = sub("[.].*", "", assembly))
meta = mutate(meta, species = gsub(" ","_", species))


# we get the data ready for the plot-----------------------------

# we create a table of the best hit per transfer, the one we will show on the tree
# we thus place hits with highest pID on top
setorder(retainedHits, -pID)

# we extract these hits, and the columns we need
connections <- retainedHits[!duplicated(hitGroup), .(species.1, species.2, superfamily, hitGroup, independent)]

# we will show different colors for the TE class, we determine them based on superfamily name
connections[, class := ifelse(test = str_detect(superfamily, "DNA/") | str_detect(superfamily, "RC/") , yes = "DNA",  no = "RNA")]

connections <- left_join(connections, meta[,data.frame(species, habitat)] , by=c("species.1"="species")) %>% left_join(., meta[,data.frame(species, habitat)] , by=c("species.2"="species"), suffix=c(".1", ".2")) 

connections <- mutate(connections, compEnv = ifelse(habitat.1==habitat.2, habitat.1, "inter"))

#Either I color the connections according to the habitat

if(colEvents == "habitat"){

    connections <- connections[
        sample(.N), # We shuffle rows to ensure that arcs will be drawn in random order,
        # to avoid visual bias between TE classes
        data.table(independent, # We retain the "independent" column to leave us the
                            # option to show independent HTTs or all hit groups.
                            # Species names are converted to integer that
                            # correspond to tip numbers on the tree.
            tip1 = match(species.1, tree$tip.label),
            tip2 = match(species.2, tree$tip.label),
            col = ifelse(compEnv == "terrestrial",
                rgb(0.4, 0.9, 0.4, 0.4), # green 
                ifelse(compEnv == "aquatic",
                    rgb(0.4, 0.4, 0.9, 0.4),   # blue
                    rgb(0.9, 0.4, 0.4, 0.4))) #ref
            
        )
    ] 
} else if(colEvents == "TEtype"){ #Or I color the connections according to the class of the TE 
    connections <- connections[
        sample(.N), # We shuffle rows to ensure that arcs will be drawn in random order,
        # to avoid visual bias between TE classes
        data.table(independent, # We retain the "independent" column to leave us the
                            # option to show independent HTTs or all hit groups.
                            # Species names are converted to integer that
                            # correspond to tip numbers on the tree.
            tip1 = match(species.1, tree$tip.label),
            tip2 = match(species.2, tree$tip.label),
            col = ifelse(test = class == "DNA",
                yes = rgb(0.9, 0.4, 0.4, 0.4), # red
                no = rgb(0.4, 0.4, 0.9, 0.4)   # blue   
            )
        )
    ] 
}




# this will be the position of species along the circle
xTipPos <- node.height(tree)

# we obtain the max depth of the tree (age of mrca)
treeDepth <- max(node.depth.edgelength(tree))

# see function. This uses ape to return lines (vertex coordinates) to draw curved branches
brancheLines <- linesFromTree(tree)

# Y position of the tree tips = distance from the center of the plot area
tipPos <- treeDepth * 1.7 + 10

# we locate clades younger than 40 My within which HTT was not inferred
# (excluding clades composed of just one species, as we don't need to outline those)
outlinedClades <- cladesOfAge(
    tree = tree,
    age = 40, 
    withTips = T,
    names = F,
    singletons = F
)

#There is an execption though: mrca 434
#Even though I write its divTime is 40Myrs, I don't know what it really is, so I looked for HTT there
outlinedClades = filter(outlinedClades, node !=434)

#I also removed all HTT among clades that did not pass filter B,
    #ie clade whose dS busco < 0.85
mrca_fitlerB = readRDS("mrca_FilterB.Rds")
edges <- data.table(tree$edge)
subClades <- edges[V1 %in% mrca_fitlerB, ]
colnames(subClades) <- c("node", "tip")
outlinedClades <- rbind(outlinedClades, subClades)

#chek that all children of a deleted clade were also deleted:
res = c()
for(i in mrca_fitlerB){
        res = append(res, getDescendants(tree,i)) 
}
res = unique(res)
filter(retainedHits, mrca %in% res & !mrca %in% mrca_fitlerB) #should be empty


pairs <- fread("Busco/pairsALL.txt")

#Add info
pairs = left_join(pairs, select(meta, -c(species,dbBusco)), by=c("assembly.1"="assembly")) %>% 
    left_join(., select(meta, -c(species,dbBusco)), by=c("assembly.2"="assembly"), suffix=c(".1", ".2"))

#MFind the deepest node for each group
taxa = filter(pairs, group.1==group.2)[,max(divTime), by=.(group.1)] %>% 
    inner_join(., select(pairs, c(group.1, divTime, mrca)), by=c("V1"="divTime", "group.1"="group.1")) %>% distinct()


# we compute the age of the mrca of the taxa, which is use to
# find the right position to write taxa names on the tree
taxa$age = taxa$V1/2

taxa$onTree = TRUE

#Give colors to each taxon. They are in the following file:
colTaxa <- fread("colorsTaxa.txt")
taxa[, col := colTaxa[match(group.1, group.1), col]]


#add capital letter
taxa$group.1 <- str_to_title(taxa$group.1) 


pdf("FigureTREE.pdf", 7, 7) 

# see circularPlots.R for the function used
initializePlot(tipPos + treeDepth, -1, max(xTipPos) + 2L, 180)

#I nned to give xOffset & yOffset: give 0 to everyone to begin with , see what plot lookks like, and correct
taxa$xOffset = 0
taxa$yOffset = 0

#ATTENTION following to set on PDF, because totally deffirent in R window!
#If all diptera in one groupe :
taxa$xOffset = c(-10,-4,-8,-3,0.5,0,1,-6,0,3,0,0,-5,-1.5,-1,2,-3,0,0,4) 
taxa$yOffset = c(50,30,0,10,0,40,0,450,0,0,55,-30,40,0,50,2,15,0,7,0) 

# Draw name clades with colors
for (i in which(taxa$onTree)) {
    with(taxa[i, ], outlineTaxon(
        tree, mrca, group.1,
        x = xTipPos[mrca] + xOffset,
        y = tipPos + age + 20 + yOffset,
        xMargin = 0.5, yMargin = 5,
        col = col,
        tipPos = tipPos - 3,
        cex = 0.5, border = NA
    ))
}


# we draw the timescale ----------------------------------------------
# the ages in My:
ages <- seq(100, 700, 100) 

# we draw the tick-marks of the ages
circSegments(
    x0 = -2.3,          # the position is just before the first tip (hence negative)
    y0 = tipPos + ages, # the "depth" of the tick is defined by the ages
    x1 = -1.7,          # to create segments of appropriate length
    curved = F,         # these segments need not be curved according to the tree
    lend = 1, lwd = 0.5
)

# we draw the ages themselves
circText(
    x = rep(-0.8, length(ages)),
    y = tipPos + ages,
    labels = ages, #Jean had / 10
    cex = 0.4,
    correct = F
)


# we draw the concave tree -------------------------------------------------
l <- lapply(
    X = brancheLines,
    FUN = function(v) {
        circLines(
            x = v$y, # x takes column y because fundamentally, this tree is drawn horizontally
            y = tipPos + treeDepth - v$x,
            lwd = 1, lend = 2, col = grey(0.3)
        )
    }
)


# we add colored points at the node of outlined clades ------------------------
taxa[
    onTree == T & age > 1, # age > 1 to avoid drawing points for single-tip taxa
    circPoints(
        x = xTipPos[mrca],
        y = tipPos + age,
        pch = 21,
        col = grey(0.3),
        bg = col,
        cex = 0.5
    )
]

# we add curved grey segments above tips of clades too young to look for HTT
outlinedClades[, circSegments(
    x0 = min(xTipPos[tip]) - 0.2,
    y0 = tipPos - 40,
    max(xTipPos[tip]) + 0.2,
    col = grey(0.5),
    lwd = 3,
    lend = 1
),
by = node
]

# we  add color for tips depending on T or A
pairs[
    sp1==sp2, # age > 1 to avoid drawing points for single-tip taxa
    circPoints(
        x = xTipPos[mrca],
        y = tipPos - 10,
        col = ifelse(habitat.1=="aquatic", "blue", "green"),
        bg = ifelse(habitat.1=="aquatic", "blue", "green"),
        cex = 0.5
    )
]

# we show HTT event with arcs. ---------------------------------------------------
# Note that we only show "independent" transfers but we could draw
# all hit groups by removing the filter on the rows if we wanted


if(colEvents = "black"){
    co <- connections[independent == T, connectionArcs(
        x0 = xTipPos[tip1],
        y0 = tipPos - 60,
        x1 = xTipPos[tip2],
        col = rgb(0,0,0,0.1), 
        lwd = 1
    )]  
} else {
    co <- connections[independent == T, connectionArcs(
        x0 = xTipPos[tip1],
        y0 = tipPos - 60,
        x1 = xTipPos[tip2],
        col = col
        lwd = 1
    )]      
}

#Add legend for color lines (when colors represent class of TE)
if(colEvents = "TEtype"){
    legend( x="bottomright", 
        legend=c("HT of class I","HT of class II"),
        col=c("#6666E666","#E6666666"), lwd=2, bty='n')
}
        
dev.off()