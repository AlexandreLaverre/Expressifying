
path <- "/Users/alaverre/Documents/Expressifying/"
EVE_liver <- read.table(paste0(path, "results/EVE_model/liver_log2TPM_betaSharedTest.txt"))

var_log <- read.table("/Users/alaverre/Library/CloudStorage/Dropbox/Shared/Expressifying/results/table_var_inter_intra/liver_logTPM.tsv", h=T)
var_log$dataset <- "liver"
tissues = c("testis", "kidney", "cerebellum")

for (tissue in tissues){
  var <- read.table(paste0("/Users/alaverre/Library/CloudStorage/Dropbox/Shared/Expressifying/results/table_var_inter_intra/", tissue, "_logTPM.tsv"), h=T)
  var$dataset <- tissue
  var_log <- rbind(var_log, var)
}

pdf(paste0(path, "results/figures/Boxplots_Variation_tissues.pdf"))
plot((table(var_log[which(var_log$ratio>1),]$dataset)/table(var_log$dataset))*100, ylab="% genes with ratio>1", main="", las=1, cex.lab=1.3, cex.axis=1.3)
boxplot(var_log$var_between~var_log$dataset, outline=F, notch=T, xlab="", ylab="Between species variation", las=1, cex.lab=1.3, cex.axis=1.3)
boxplot(var_log$var_within~var_log$dataset, outline=F, notch=T,xlab="", ylab="Within species variation", las=1, cex.lab=1.3, cex.names=1.3)
boxplot(var_log$ratio~var_log$dataset, outline=F, notch=T,xlab="", ylab="Ratio", las=1, cex.lab=1.3, cex.names=1.3)
dev.off()


Bayes_rank <- read.table("/Users/alaverre/Library/CloudStorage/Dropbox/Shared/Expressifying/results/table_posterior_probs/merge_Bayes_score.tsv", h=T)
Bayes_log <- read.table("/Users/alaverre/Library/CloudStorage/Dropbox/Shared/Expressifying/results/table_posterior_probs/merge_Bayes_logTPM.tsv", h=T)


for (tissue in unique(Bayes_log$dataset)){
  genes.in.tissue <- Bayes_log[which(Bayes_log$dataset==tissue), "trait"]
  write.table(genes.in.tissue, file=paste0(path, "results/gene_list/mammals_", tissue, "_ratio_sup_1.txt"), row.names = F, quote=F, col.names=F)
}

Bayes_log <- Bayes_log[order(Bayes_log$ratio,decreasing=TRUE),]
write.table(Bayes_log$trait, file=paste0(path, "results/gene_list/mammals_ordered_ratio_sup_1.txt"), row.names = F, quote=F, col.names=F)


orthogroups_human <- orthogroups_all_sp[colnames(express_ortho_all_sp)[2:nb_ortho],"Homo_sapiens"]
orthogroups_human <- orthogroups_human[!is.na(orthogroups_human)]
write.table(orthogroups_human, file=paste0(path, "results/gene_list/mammals_", tissue, "_expressed_1-1_orthogroups.txt"), row.names = F, quote=F, col.names=F)


par(mfrow=c(2,2))
par(mai=c(0.6,0.6,0.3,0.3), mgp=c(2.2,0.8,0))
plot(table(Bayes_log$dataset), ylab="Nb genes with p>1", main="log(TPM)", las=2)
boxplot(Bayes_log$var_between~Bayes_log$dataset, outline=F, notch=T, xlab="", ylab="Between var", las=2)
boxplot(Bayes_log$var_within~Bayes_log$dataset, outline=F, notch=T,xlab="", ylab="Within var", las=2)
boxplot(Bayes_log$ratio~Bayes_log$dataset, outline=F, notch=T,xlab="", ylab="Ratio", las=2)

par(mfrow=c(2,2))
plot(table(Bayes_rank$dataset), ylab="Nb genes with p>1", main="Rank", las=2)
boxplot(Bayes_rank$var_between~Bayes_rank$dataset, outline=F, notch=F, xlab="", ylab="Between var", las=3)
boxplot(Bayes_rank$var_within~Bayes_rank$dataset, outline=F, notch=F, xlab="", ylab="Within var", las=3)
boxplot(Bayes_rank$ratio~Bayes_rank$dataset, outline=F, notch=F, xlab="", ylab="Ratio", las=2)

# plot likelihood ratio test statistic histogram
hist(EVE_results$LRT,freq = F)

# Plot the chi-squared distribution with one degree of freedom
x = seq(0.5,15,length.out = 100)
y = dchisq(x,df = 1)
lines(x,y,col="red")

pval = pchisq(res$LRT,df = 1,lower.tail = F)
hist(p.adjust(pval), breaks=50, xlab="adjusted p-values", main="Liver 2,702 genes")

####
EVE_liver_signif <- EVE_liver[which(EVE_liver$p.adjusted < 0.05),]

common <- Bayes[which(EVE_liver_signif$Human_ID %in% Bayes$trait),]

Bayes_liver_signif <- Bayes[which(Bayes$dataset == "liver" & Bayes$pp_ratio_greater_1 > 0.95),]
