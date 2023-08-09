library(progress)
path = "/Users/alaverre/Documents/Expressifying/"

all_orthogroups = read.csv(paste0(path, "data/gene_orthologies/one2one_orthogroups.csv"), row.names = 1)
tissues <- c("testis", "kidney", "cerebellum")
min_species=15

for(tissue in tissues){
  print(tissue)
  
  gene.expression = readRDS(paste0(path, "data/gene_expression/mammals_", tissue, "_gene_expression_logTPM.Rds"))
  gene.score = readRDS(paste0(path, "data/gene_expression/mammals_", tissue, "_gene_expression_score.Rds"))
  
  # Correct dog species name 
  dog_index <- which(names(gene.expression) == "Canis_lupus familiaris")
  names(gene.expression)[dog_index] <- "Canis_lupus_familiaris"
  names(gene.score)[dog_index] <- "Canis_lupus_familiaris"
  
  # Remove the naked-mole rat because OMA annotations (male) are not the same as Bgee (female) 
  gene.expression[["Heterocephalus_glaber"]] <- NULL
  gene.score[["Heterocephalus_glaber"]] <- NULL
  
  species = names(gene.expression)
  orthogroups = all_orthogroups[,species]
  
  # Select orthogroups with genes for at least 15 species or all species
  orthogroups_all_sp = orthogroups[rowSums(!is.na(orthogroups)) >= min_species,]
  #orthogroups_all_sp =  orthogroups[complete.cases(orthogroups),]
  
  # Create empty data.frame for gene expression in each orthogroup
  samples <- unlist(sapply(gene.expression, function(x) colnames(x)))
  orthogroup_IDs <- row.names(orthogroups_all_sp)
  expression_ortho <- data.frame(matrix(ncol=length(orthogroup_IDs)+1, nrow=length(samples)))
  
  colnames(expression_ortho) <- c("species", orthogroup_IDs)
  row.names(expression_ortho) <- samples
  score_ortho <- expression_ortho
  
  ### Complete the data.frame
  # Create a progress bar
  pb <- progress_bar$new(total = length(orthogroup_IDs), format = "[:bar] :percent :elapsed")
  
  for (orthogroup in orthogroup_IDs){
    for (sp in colnames(orthogroups_all_sp)){
      gene = as.character(orthogroups_all_sp[orthogroup,sp])
      samples = colnames(gene.expression[[sp]])
      
      # Get the expression of this gene if present in data
      if (gene %in% rownames(gene.expression[[sp]])){
        expression = gene.expression[[sp]][gene,]
        score = gene.score[[sp]][gene,]
      }else{
        expression = NA
        score = NA
      }
      
      expression_ortho[samples, orthogroup] = expression
      expression_ortho[samples, "species"] = sp
      score_ortho[samples, orthogroup] = score
      score_ortho[samples, "species"] = sp
    }
    pb$tick()
  }
  
  ## Filters to get gene expression for each species in each orthogroup
  # Count the number of non-NA values per species
  species_without_na <- aggregate(!is.na(expression_ortho[, -1]), by = list(expression_ortho$species), FUN = sum)
  
  # Get orthogroups with at least 2 samples for each species
  nb_sp_ortho <- apply(species_without_na[, -1], 2, function(x) sum(x >= 2, na.rm = T))
  ortho_all_sp <- names(nb_sp_ortho[which(nb_sp_ortho >= min_species)])
  
  express_ortho_all_sp <- expression_ortho[,c("species", ortho_all_sp)]
  score_ortho_all_sp <- score_ortho[,c("species", ortho_all_sp)] 
  
  # Print summary and save final data.frame per tissue
  print(paste(length(species), "species;", 
              nrow(express_ortho_all_sp), "samples;", 
              ncol(express_ortho_all_sp)-1, "genes."))
  
  # add sample as column
  express_ortho_all_sp <- cbind(sample=row.names(express_ortho_all_sp), express_ortho_all_sp)
  score_ortho_all_sp <- cbind(sample=row.names(score_ortho_all_sp), score_ortho_all_sp)
  
  write.csv(express_ortho_all_sp, file=paste0(path, "results/gene_expression/mammals_", tissue, "_gene_expression_logTPM_orthogroups.csv"), row.names = F)
  write.csv(score_ortho_all_sp, file=paste0(path, "results/gene_expression/mammals_", tissue, "_gene_expression_score_orthogroups.csv"), row.names = F)
}



### All species
# liver = 24 species, 585 samples, 481 genes.
# cerebellum = 21 species, 112 samples, 577 genes.
# kidney = 21 species, 153 samples, 600 genes. 
# testis = 19 species, 160 samples, 750 genes.

### At least 15 species 
# liver: 24 species; 585 samples; 8936 genes.
# cerebellum: 21 species; 112 samples; 7676 genes.
# kidney: 21 species; 153 samples; 7958 genes.
# testis: 19 species; 160 samples; 7027 genes.

sp = "Homo_sapiens"
gene = "ENSG00000109606"
orthogroup = row.names(orthogroups[which(orthogroups[[sp]] == gene),])

# Directly from orthogroup
expression = expression_ortho[[orthogroup]]

par(mar = c(7, 4, 2, 2) + 0.1)
a <- boxplot(expression~expression_ortho$species, outline=F, xaxt = "n",
             ylab="Gene expression (FPKM)", xlab="", main=paste("DHX15 in", tissue), cex.lab=1)


axis(1, at = 1:length(a$names), labels = NA, cex.axis = 1, srt=45, tck = -0.02)
text(x = 1:length(a$names), y = par("usr")[3] - 3, labels = a$names, srt = 45, adj = 1, xpd = TRUE, cex=0.8)
