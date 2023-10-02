library(progress)
path = "/Users/alaverre/Documents/Expressifying/"

all_orthogroups = read.csv(paste0(path, "data/gene_orthologies/one2one_orthogroups_2022_mammals.csv"), row.names = 1)
tissues <- c("adult_mammalian_kidney_male", "liver_female", "liver_male", "testis_male")
min_species=10

for(tissue in tissues){
  print(tissue)
  
  gene.expression = readRDS(paste0(path, "data/gene_expression/filtered/mammals_", tissue, "_gene_expression_log2TPM.Rds"))
  gene.score = readRDS(paste0(path, "data/gene_expression/filtered/mammals_", tissue, "_gene_expression_score.Rds"))
  
  # Correct dog species name 
  dog_index <- which(names(gene.expression) == "Canis_lupus familiaris")
  names(gene.expression)[dog_index] <- "Canis_lupus_familiaris"
  names(gene.score)[dog_index] <- "Canis_lupus_familiaris"
  
  # Remove the naked-mole rat because OMA annotations (male) are not the same as Bgee (female) 
  gene.expression[["Heterocephalus_glaber"]] <- NULL
  gene.score[["Heterocephalus_glaber"]] <- NULL
  
  species = names(gene.expression)
  orthogroups = all_orthogroups[,species]
  
  # Select orthogroups with genes for at least 10 species or all species
  orthogroups_all_sp = orthogroups[rowSums(!is.na(orthogroups)) >= min_species,]
  #orthogroups_all_sp =  orthogroups[complete.cases(orthogroups),]
  
  # Background genes
  orthogroups_human <- orthogroups_all_sp$Homo_sapiens
  orthogroups_human <- orthogroups_human[!is.na(orthogroups_human)]
  write.table(orthogroups_human, file=paste0(path, "results/gene_list/mammals_", tissue, "_1-1_orthogroups_filtered.txt"), row.names = F, quote=F, col.names=F)
  
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
  nb_ortho = ncol(express_ortho_all_sp)
  print(paste(length(species), "species;", 
              nrow(express_ortho_all_sp), "samples;", 
              nb_ortho-1, "genes."))
  
  # Background expressed genes
  orthogroups_human <- orthogroups_all_sp[colnames(express_ortho_all_sp)[2:nb_ortho],"Homo_sapiens"]
  orthogroups_human <- orthogroups_human[!is.na(orthogroups_human)]
  write.table(orthogroups_human, file=paste0(path, "results/gene_list/mammals_", tissue, "_expressed_filtered.txt"), row.names = F, quote=F, col.names=F)
  
  # add sample as column
  express_ortho_all_sp <- cbind(sample=row.names(express_ortho_all_sp), express_ortho_all_sp)
  score_ortho_all_sp <- cbind(sample=row.names(score_ortho_all_sp), score_ortho_all_sp)
  
  write.csv(express_ortho_all_sp, file=paste0(path, "results/gene_expression/filtered/mammals_", tissue, "_gene_expression_log2TPM_orthogroups_2022_only_mammals.csv"), row.names = F)
  write.csv(score_ortho_all_sp, file=paste0(path, "results/gene_expression/filtered/mammals_", tissue, "_gene_expression_score_orthogroups_2022_only_mammals.csv"), row.names = F)
}
