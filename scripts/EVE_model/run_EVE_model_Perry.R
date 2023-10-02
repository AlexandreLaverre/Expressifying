library(evemodel)
library(ape)

path <- "/Users/alaverre/Documents/Expressifying/data/Perry - Supplemental.Database.primateRNAseq/Gene.expression.data/"

################################################################################
# Gene expression
gene.express <- read.table(paste0(path, "Normalized.expression.data.txt"), h=T, row.names = 1)
gene.express <- gene.express[,33:ncol(gene.express)] # remove Length and GC columns

exprMat <- as.matrix(gene.express)
colSpecies <- sub("_.*$","",colnames(exprMat))
message("Original data: ", nrow(exprMat), " genes, ", ncol(exprMat), " samples and ", length(unique(colSpecies)), " species.")

# Filters
# Reduce missing values
exprMat <- exprMat[rowSums(is.na(exprMat)) < 4, ] # keep gene with less than 4 samples with NA 
#exprMat <- exprMat[1:1000,] # Keep the first 1000 genes

# Vector of species 
colSpecies <- sub("_.*$","",colnames(exprMat))
message("Filtered data: ", nrow(exprMat), " genes,", ncol(exprMat), " samples and ", length(unique(colSpecies)), " species.")
table(colSpecies)

################################################################################
# Species Tree from Zoonomia
speciesTree <- read.tree("/Users/alaverre/Documents/Detecting_positive_selection/data/species_trees/241-mammals.nk")

missing <- c("Monodelphis_domestica", "Galago moholi", " Varecia variegata", "Eulemur mongos", "Eulemur coronatus")
# remove missing : Opossum, Galago, BWLemur
exprMat <- subset(exprMat, select = -c(43:46, 19:22, 7:10))
colSpecies <- sub("_.*$","",colnames(exprMat))

sp <- c("Monodelphis_domestica", "Mus_musculus", "Dasypus_novemcinctus", "Tupaia_chinensis",  "Nycticebus_coucang",
        "Daubentonia_madagascariensis", "Propithecus_coquereli", "Eulemur_fulvus", "Eulemur_flavifrons", "Callithrix_jacchus",
        "Chlorocebus_sabaeus", "Macaca_mulatta", "Pan_troglodytes", "Homo_sapiens")

# keep only species with gene expression data
common.species <- intersect(sp, speciesTree$tip.label)
speciesTree <- keep.tip(speciesTree, common.species)

used.species.names <- c("SlowLoris", "AyeAye", "Sifaka", "MongooseLemur", "CrownedLemur",
                        "Marmoset", "Macaque", "Vervet", "Chimp", "Human", "TreeShrew", "Mouse", "Armadillo")

speciesTree$tip.label <- used.species.names
plot(speciesTree)

################################################################################
# Running beta shared test
res <- betaSharedTest(tree = speciesTree, gene.data = exprMat, colSpecies = colSpecies)
saveRDS(res, file=paste0(path, "betaSharedTest.RDS"))

# plot likelihood ratio test statistic histogram
hist(res$LRT,freq = F)

# Plot the chi-squared distribution with one degree of freedom
x = seq(0.5,15,length.out = 100)
y = dchisq(x,df = 1)
lines(x,y,col="red")

pval= pchisq(res$LRT,df = 1,lower.tail = F)
hist(p.adjust(pval))


