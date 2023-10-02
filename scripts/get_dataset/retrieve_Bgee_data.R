library(BgeeDB)

################################################################################

path <- "/Users/alaverre/Documents/Expressifying/"

# BgeeDB functions to properly retrieve gene expression rank
source(paste0(path, "/scripts/get_dataset/modified_bgee_functions.R"))

# Species list from Tarcisio and manually annotated for "Taxa"
species <- read.csv(paste0(path, "/data/Bgee_libraries/species_list.csv"), header=T, fill=T)
mammals <- species[which(species$Large_Taxa == "Mammal" | species$Large_Taxa == "Marsupial" ), "complete_name"]
#mammals[mammals != "Monodelphis_domestica"]

################################################################################
## Prepare samples to download 
# Get numbers of experiments per UBERON IDs
IDs <- list()
for (sp in mammals){
  bgee <- Bgee$new(species=sp, dataType = "rna_seq", pathToData=paste0(path, "/data/Bgee_libraries/"))
  annotation_samp <- getAnnotation(bgee)$sample.annotation
  annotation_samp$ID <- paste(annotation_samp$Anatomical.entity.name, annotation_samp$Anatomical.entity.ID, annotation_samp$Sex)
  
  IDs[[sp]] <- as.data.frame(table(annotation_samp$ID))
  colnames(IDs[[sp]]) <- c("ID", sp)
}

# Merging tables
IDs_count <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), IDs)
rownames(IDs_count) <- IDs_count$ID
IDs_count <- IDs_count[, -1]

# Counting the number of species with at least 2 libraries
IDs_count$Nb_sp_var <- apply(IDs_count, 1, function(x) sum(x >= 2, na.rm = T))

# Top tissue with maximum number of species
IDs_count <- IDs_count[order(IDs_count$Nb_sp_var, decreasing=T),]
head(IDs_count[,c(1,28)])

# Get species list for selected tissues (from head(IDs_count))
liver_variance <- head(colnames(IDs_count[which(IDs_count[1,] > 1)]), -1) # UBERON:0002107 
adult_kidney_variance <- head(colnames(IDs_count[which(IDs_count[2,] > 1)]), -1) # UBERON:0000082  
cerebellum_variance <- head(colnames(IDs_count[which(IDs_count[3,] > 1)]), -1) #UBERON:0002037
testis_variance  <- head(colnames(IDs_count[which(IDs_count[4,] > 1)]), -1) # UBERON:0000473

species_list <- list(liver_variance, adult_kidney_variance, cerebellum_variance, testis_variance)
names(species_list) <- c("liver", "kidney", "cerebellum", "testis")
UBERON <- c("UBERON:0002107", "UBERON:0000082", "UBERON:0002037", "UBERON:0000473")
names(UBERON) <- c("liver", "kidney", "cerebellum", "testis") 

################################################################################
# Retrieve gene expression data for all species in selected tissue 
for (tissu in names(species_list)){
  print(tissu)
  output.expression <- paste0(path, "/data/gene_expression/mammals_", tissu, "_gene_expression_log2TPM.Rds")
  output.score <- paste0(path, "/data/gene_expression/mammals_", tissu, "_gene_expression_score.Rds")
  
  if (file.exists(output.expression)){
    print("Already done!")
  }else{
  
    gene.expression <- list()
    gene.score <- list()
    for (sp in species_list[[tissu]]){
      print(sp)
      bgee <- Bgee$new(species=sp, dataType = "rna_seq", pathToData=paste0(path, "/data/Bgee_libraries/"))
      data <- getData(bgee, anatEntityId = UBERON[[tissu]])
      
      data$Rank <- as.numeric(data$Rank)
      simplified.data <- formatData(bgee, data, callType = "present", stats = "tpm")
      all_gene_expression <- log2(0.01+simplified.data@assayData[["exprs"]])
      
      simplified.data <- formatData(bgee, data, callType = "present", stats = "rank")
      all_gene_rank <- simplified.data@assayData[["exprs"]]
      
      # New score based on rank normalized per library
      all_gene_rank_normalised <- apply(all_gene_rank, 2, function(x) 100-(x*100)/max(x, na.rm=T))
      
      # Remove genes with less than 2 measures
      gene.expression[[sp]] <- all_gene_expression[rowSums(!is.na(all_gene_expression)) >= 2,]
      gene.score[[sp]] <- all_gene_rank_normalised[rowSums(!is.na(all_gene_rank_normalised)) >= 2,]
      
    }
    saveRDS(gene.expression, file=output.expression)
    saveRDS(gene.score, file=output.score)
  }
}

################################################################################
