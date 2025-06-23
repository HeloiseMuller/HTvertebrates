library(ggplot2)

path = "~/Project/"
setwd(path) 

source("HTvFunctions.R")

tree = read.tree("datasetTree.nwk")
htt = fread("nbHTTevents_perPairAll.txt")
mrcas = mrca(tree)
htt[,mrca := mrcas[cbind(species.1, species.2)]]


# comparison intra vs inter-habitat  ----------------
species = fread("metadata.txt")
pairs = htt[,.(same = sum(habitat.1 == habitat.2), different = sum(habitat.1 != habitat.2), n = sum(n)), by = mrca]
mrcas = pairs[same >=1 & different >=1 & n > 0, mrca]

htt[,same := habitat.1 == habitat.2]
means = htt[mrca %in% mrcas,.(intra = mean(n[same]), inter = mean(n[!same])), by = mrca]
means[,wilcox.test(intra, inter, paired = T)]

intraColor = "deepskyblue3"
interColor = "tomato"

divTime = max(node.depth.edgelength(tree))*2 - node.depth.edgelength(tree)*2
yPos = as.factor(rank(divTime[means$mrca]))
breaks = as.character(seq(1, length(yPos), 2))

ggplot(means, aes(y = yPos, yend = yPos, col = intra > inter)) + scale_color_manual(values = c(interColor, intraColor)) +
  geom_segment(aes(x = intra, xend = inter), lwd = 1, show.legend = F) + ylab("MRCA") + scale_y_discrete(breaks = breaks) +
  xlab("Mean number of transfers per species pair") +  theme_minimal() + theme(axis.text=element_text(size=8))

ggsave('figure3b.pdf', width = 4, height = 4)

means[,wilcox.test(intra, inter, paired = T, alternative = "greater")]
means[,binom.test(table(intra < inter), alternative = "greater")]

# showing MRCAs on the tree
pdf("figureS7.pdf", width = 7, height = 15)
par(mai = c(0.1, 0.1, 0.1, 0.1))
plot(tree, show.tip.label = T, cex = 0.4, label.offset = 5)
X =  max(node.depth.edgelength(tree))
Y = node.height(tree)[1:length(tree$tip.label)]
cols = species[match(tree$tip.label, species), ifelse(habitat == "aquatic", "lightskyblue", "darkolivegreen3")]
points(rep(X, length(Y)), Y, pch = 19, col = cols, cex = 0.5)

Y = node.height(tree)[means$mrca]
X = node.depth.edgelength(tree)[means$mrca] 
cols = means[,ifelse(intra > inter, intraColor , ifelse(intra < inter,interColor, "grey"))]
points(X, Y, pch = 19, col = cols, cex = 2)
text(X, Y, labels = yPos, cex= 0.8)

dev.off()

