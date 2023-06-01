library(BgeeDB)
library(XML)

path <- "/Users/alaverre/Documents/Expressifying/"
species <- read.csv(paste0(path, "/data/Bgee_libraries/species_list.csv"), header=T, fill=T)
mammals <- species[which(species$Large_Taxa == "Mammal" | species$Large_Taxa == "Marsupial" ), "complete_name"]
#mammals[mammals != "Monodelphis_domestica"]

################################################################################
# Get numbers of experiments per IDs
IDs <- list()
for (sp in mammals){
  bgee <- Bgee$new(species=sp, dataType = "rna_seq", pathToData=paste0(path, "/data/Bgee_libraries/"))
  
  annotation_bgee <- getAnnotation(bgee)$sample.annotation
  annotation_bgee$ID <- paste(annotation_bgee$Anatomical.entity.name, annotation_bgee$Anatomical.entity.ID)
  
  IDs[[sp]] <- as.data.frame(table(annotation_bgee$ID))
  colnames(IDs[[sp]]) <- c("ID", sp)
}

# Merging tables
IDs_count <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), IDs)
rownames(IDs_count) <- IDs_count$ID
IDs_count <- IDs_count[, -1]

# Counting the number of species with at least 2 experiments
IDs_count$Nb_sp_var <- apply(IDs_count, 1, function(x) sum(x > 1, na.rm = T))
IDs_count <- IDs_count[order(IDs_count$Nb_sp_var, decreasing=T),]
head(IDs_count)

# Get species list for tissues with the highest variance
liver_variance <- head(colnames(IDs_count[which(IDs_count[1,] > 1)]), -1) # UBERON:0002107 
adult_kidney_variance <- head(colnames(IDs_count[which(IDs_count[2,] > 1)]), -1) # UBERON:0000082  
cerebellum_variance <- head(colnames(IDs_count[which(IDs_count[3,] > 1)]), -1) #UBERON:0002037
testis_variance  <- head(colnames(IDs_count[which(IDs_count[4,] > 1)]), -1) # UBERON:0000473

species_list <- list(liver_variance, adult_kidney_variance, cerebellum_variance, testis_variance)
names(species_list) <- c("liver", "kidney", "cerebellum", "testis")
UBERON <- c("UBERON:0002107", "UBERON:0000082", "UBERON:0002037", "UBERON:0000473")
names(UBERON) <- names(species_list) 

################################################################################
# Retrieve gene expression data for each tissue
for (tissu in names(species_list)){
  print(tissu)
  gene.expression <- list()
  for (sp in species_list[[tissu]]){
    print(sp)
    bgee <- Bgee$new(species=sp, dataType = "rna_seq", pathToData=paste0(path, "/data/Bgee_libraries/"))
    data <- getData(bgee, anatEntityId = UBERON[[tissu]])
    simplified.data <- formatData(bgee, data, callType = "present", stats = "fpkm")
    all_gene_expression <- simplified.data@assayData[["exprs"]]
    
    # Remove genes with less than 2 measures
    gene.expression[[sp]] <- all_gene_expression[rowSums(!is.na(all_gene_expression)) >= 2,]
    
  }
  saveRDS(gene.expression, file=paste0(path, "/data/gene_expression/mammals_", tissu, "_gene_expression.Rds"))
}

################################################################################