library(ggplot2)
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(dplyr)
library(RColorBrewer)
library(heatmaply)

########################################################################
path = "/Users/alaverre/Documents/Expressifying/results/"
tissues <- c("cerebellum", "testis", "liver", "kidney")

measures=c("logTPM", "score") # logTPM or score

for (measure in measures){
  
# Get all data
scores <- list()
median_scores <- list()
Species <- c()
for (tissu in tissues){
  score = read.csv(paste0(path, "gene_expression/mammals_", tissu, "_gene_expression_", measure, "_orthogroups.csv"), row.names = 1)
  species <- score$species
  row.names(score) = paste(species, row.names(score))
  score = score[,-1]
  scores[[tissu]] <- score
  median = as.data.frame(apply(score, 2, function(x) tapply(x, species, function(y) median(y, na.rm=T))))
  row.names(median) <- paste(row.names(median), tissu)
  median_scores[[tissu]] <- median
  
  Species <- c(Species, species)
}

########################################################################
### ACP all samples
all_scores <- bind_rows(scores, .id = "column_label")
Tissues <- all_scores$column_label
res.pca <- PCA(all_scores[,-1], graph = F)

pdf(paste0(path, "figures/ACP_mammals_", measure, ".pdf"))
p <- fviz_pca_ind(res.pca, habillage = as.factor(Tissues), label="none", addEllipses=T, ellipse.level=0.95) + theme_minimal()
print(p)

# ACP with colors and shapes
basic_plot <- fviz_pca_ind(res.pca)
ggplot(cbind(basic_plot$data, Species, Tissues), aes(x=x,y=y,col=Species,shape=Tissues)) + stat_ellipse(aes(x=x,y=y,col=Tissues)) + geom_point() + theme_minimal() 

dev.off()

### ACP median per sp and tissu 
all_scores_median <- bind_rows(median_scores) #, .id = "column_label")

res.pca.med <- PCA(all_scores_median, graph = F)
basic_plot <- fviz_pca_ind(res.pca.med)
Species_med = sapply(strsplit(row.names(basic_plot$data), " "), `[[`, 1)
Tissues_med = sapply(strsplit(row.names(basic_plot$data), " "), `[[`, 2)

pdf(paste0(path, "figures/ACP_mammals_", measure, "_median.pdf"))
p <- fviz_pca_ind(res.pca.med, habillage = as.factor(Tissues_med), label="none", addEllipses=T,
                  ellipse.level=0.95) + theme_minimal()
print(p)
# ACP with colors and shapes
basic_plot <- fviz_pca_ind(res.pca.med)
ggplot(cbind(basic_plot$data, Species_med, Tissues_med), aes(x=x,y=y,col=Species_med,shape=Tissues_med)) + stat_ellipse(aes(x=x,y=y,col=Tissues_med)) + geom_point() + theme_minimal() 
dev.off()

########################################################################
# Heatmap 
# Models correlations clustering
get_corr <- function(data){
  rcorrDat <- Hmisc::rcorr(t(data), type="pearson")
  cormat <- rcorrDat$r
  
  # Reorder the matrix: use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  return(cormat)
}

mycolors <- colorRampPalette(c("dodgerblue4", "white", "firebrick"))
# Plot Heatmap
heatmaply(get_corr(all_scores[,-1]), k_row = 4, k_col = 4, fontsize_row=4, fontsize_col = 4, 
          colors = mycolors, colorbar_len =0.2, 
          file=paste0(path, "figures/heatmap_mammals_", measure, "_dendrogram.pdf"))

heatmaply(get_corr(all_scores_median), k_row = 4, k_col = 4, fontsize_row=4, fontsize_col = 4, 
          colors = mycolors, colorbar_len =0.2, 
          file=paste0(path, "figures/heatmap_mammals_", measure, "_median_dendrogram.html"))
########################################################################
}

