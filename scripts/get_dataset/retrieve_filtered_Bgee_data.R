library(BgeeDB)
library(dplyr)
library(ape)

path <- "/Users/alaverre/Documents/Expressifying/"

# BgeeDB functions to properly retrieve gene expression rank
source(paste0(path, "/scripts/get_dataset/modified_bgee_functions.R"))

################################################################################
# Bgee annotations
bgee.species <- listBgeeSpecies(release="15.0")
bgee.species$Species <- paste(bgee.species$GENUS, bgee.species$SPECIES_NAME, sep="_")

annotation <- lapply(as.character(bgee.species$ID), function(id) 
                     getAnnotation(Bgee$new(species=id, dataType="rna_seq",
                                            pathToData=paste0(path, "/data/Bgee_libraries/"),
                                            release="15.0")))
names(annotation) <- bgee.species$Species

# Get experiments and samples annotations
annotation.exp <- lapply(annotation, function(i) i$experiment.annotation) %>% bind_rows(.id="Species")
annotation.sample <- lapply(annotation, function(i) i$sample.annotation) %>% bind_rows(.id="Species")

################################################################################
## Filters
# Remove experiments with confounding factors
confounding.factors <- "circadian|diurnal|time.*point|before.*after|feeding|diet|sectioning|pooled"

exp.remove <- annotation.exp$Experiment.ID[
              grepl(confounding.factors, ignore.case=TRUE, annotation.exp$Experiment.name) |
              grepl(confounding.factors, ignore.case=TRUE, annotation.exp$Experiment.description)]

# Keep only post-juvenile (UBERON:0000113 and descendants) or "life cycle" (UBERON:0000104) samples
post.juvenile.terms <- read.csv(paste0(path, "/data/Bgee_libraries/desc_UBERON_0000113.csv"), header=TRUE)
colnames(post.juvenile.terms) <- c("Stage.ID", "Stage.name")

post.juvenile.terms$Stage.ID <- gsub("http://purl.obolibrary.org/obo/","", post.juvenile.terms$Stage.ID) %>% gsub("_",":",.)

# Remove also experiment without Sex information
samp.remove <-  annotation.sample[
                !(annotation.sample$Stage.ID %in% c("UBERON:0000104", unlist(post.juvenile.terms$Stage.ID))) 
                | annotation.sample$Experiment.ID %in% exp.remove
                | is.na(annotation.sample$Sex), "Library.ID"]

annotation.sample.adult <- annotation.sample[!(annotation.sample$Library.ID %in% samp.remove),]

################################################################################
## Number of libraries per species - organ - sex
annotation.sample.adult$ID <- paste(annotation.sample.adult$Anatomical.entity.name,
                                    annotation.sample.adult$Anatomical.entity.ID, 
                                    annotation.sample.adult$Sex, sep="_")

list.count.by.ID <- tapply(annotation.sample.adult$ID, as.factor(annotation.sample.adult$Species),
               function(x) as.data.frame(table(x, dnn="ID")))

# Reduce in one data frame 
count.by.ID <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), list.count.by.ID)
colnames(count.by.ID) <- c("ID", names(list.count.by.ID))
rownames(count.by.ID) <- count.by.ID$ID
count.by.ID <- count.by.ID[, -1]

# Remove species not in Zoonomia
Zoonomia.tree <- read.tree("/Users/alaverre/Documents/Detecting_positive_selection/data/species_trees/241-mammals.nk")
common.species <- c(intersect(colnames(count.by.ID), Zoonomia.tree$tip.label), "Canis_lupus familiaris")

count.by.ID <- count.by.ID[,common.species]

# Get species list for organ-sex with at least 2 libraries per species
species.lists <- apply(count.by.ID, 1, function(x) colnames(count.by.ID[which(x > 1)]))

################################################################################
# Retrieve gene expression data for each organ-sex containing at least 10 species 
for(i in 1:length(species.lists)){
  if (length(species.lists[[i]]) > 10){
    ID = names(species.lists)[i]
    species=species.lists[[i]]
    print(paste(ID, length(species), "species."))
    
    tissue=gsub(" ", "_", unlist(strsplit(ID, "_"))[1])
    UBERON=unlist(strsplit(ID, "_"))[2]
    sex=unlist(strsplit(ID, "_"))[3]
    
    output.expression <- paste0(path, "/data/gene_expression/filtered/mammals_", tissue, "_", sex, "_gene_expression_log2TPM.Rds")
    output.score <- paste0(path, "/data/gene_expression/filtered/mammals_", tissue, "_", sex, "_gene_expression_score.Rds")
    
    if (file.exists(output.expression)){
      print("Already done!")
    }else{
      
      gene.expression <- list()
      gene.score <- list()
      for (sp in species){
        print(sp)
        bgee <- Bgee$new(species=sp, dataType = "rna_seq", pathToData=paste0(path, "/data/Bgee_libraries/"))
        data <- getData(bgee, anatEntityId=UBERON, sex=sex)
        
        # Remove filtered samples
        data <- data[which(!data$Library.ID %in% samp.remove),]
        
        # Retrieve log2(TPM)
        simplified.data <- formatData(bgee, data, callType = "present", stats = "tpm")
        all_gene_expression <- log2(0.01+simplified.data@assayData[["exprs"]])
        
        # Calculate a score based on rank normalized per library
        data$Rank <- as.numeric(data$Rank)
        simplified.data <- formatData(bgee, data, callType = "present", stats = "rank")
        all_gene_rank <- simplified.data@assayData[["exprs"]]
        all_gene_rank_normalised <- apply(all_gene_rank, 2, function(x) 100-(x*100)/max(x, na.rm=T))
        
        # Remove genes with less than 2 measures
        gene.expression[[sp]] <- all_gene_expression[rowSums(!is.na(all_gene_expression)) >= 2,]
        gene.score[[sp]] <- all_gene_rank_normalised[rowSums(!is.na(all_gene_rank_normalised)) >= 2,]
        
      }
      warsaveRDS(gene.expression, file=output.expression)
      saveRDS(gene.score, file=output.score)
    }
  }
}

################################################################################
    
  
