path = "/Users/alaverre/Documents/Expressifying/"

tissue="liver"

# Expression data
express <- read.csv(paste0(path, "results/gene_expression/mammals_", tissue, "_gene_expression_log2TPM_orthogroups.csv"), h=T)
score <- read.csv(paste0(path, "results/gene_expression/mammals_", tissue, "_gene_expression_score_orthogroups.csv"), h=T)

# Gene orthology
ortho <- read.csv(paste0(path, "data/gene_orthologies/one2one_orthogroups.csv"), row.names = 1)

################################################################################
#### Expression of a given orthogroup in all species 
orthogroup = sample(colnames(express), 1)
human.gene = ortho[orthogroup, "Homo_sapiens"]

par(mfrow=c(1,1))
par(mai=c(1.5,0.8,0.4,0.3), mgp=c(2.2,0.7,0))
a <- boxplot(express[[orthogroup]]~express$species, outline=F, xaxt = "n", las=1,
             ylab="log(TPM)", xlab="", main=paste(orthogroup, "in", tissue), cex.lab=1)

axis(1, at = 1:length(a$names), labels = NA, cex.axis = 1, srt=45, tck = -0.02)
text(x = 1:length(a$names), y = par("usr")[3] - 0.3, labels = a$names, srt = 45, adj = 1, xpd = TRUE, cex=1)

################################################################################
#### Expression of a given species 
sp="Homo_sapiens"

nb_ortho = length(colnames(express))
expression.sp = t(express[which(express$species == sp),3:nb_ortho])
score.sp = t(score[which(score$species == sp),3:nb_ortho])
colnames(expression.sp) <- express[which(express$species == sp),1]
colnames(score.sp) <- express[which(express$species == sp),1]

median_expression <- apply(expression.sp, 1, function(x) median(x, na.rm=T))
median_score <- apply(score.sp, 1, function(x) median(x, na.rm=T))

sd_expression <- apply(expression.sp, 1, function(x) sd(x, na.rm=T))
sd_score <- apply(score.sp, 1, function(x) sd(x, na.rm=T))

# all samples
#plot(expression.sp~score.sp,cex=0.05, xlab="Expression score", ylab="log(TPM)", main=paste(sp, tissue))

# one sample
par(mfrow=c(1,1))
par(mai=c(1.5,0.8,0.4,0.3), mgp=c(2.2,0.7,0))
sample = sample(colnames(expression.sp),1)
plot(expression.sp[,sample]~score.sp[,sample],cex=0.05, xlab="Normalised rank score", ylab="log(TPM)", las=1, cex.lab=1.2,
     main=paste(sp, tissue, sample))

test <- cor.test(expression.sp[,25],score.sp[,25])
mtext(paste("R =",signif(test$estimate,3)), line=-2, cex=1.2)

# 10 random genes
orthogroups = c(orthogroup, sample(rownames(expression.sp), 9))
par(mai=c(1.3,0.8,0.5,0.3), mgp=c(2.2,1,0))

par(mfrow=c(2,1))
par(mai=c(0.4,0.8,0.1,0.3), mgp=c(2.2,0.8,0))
a <- boxplot(t(expression.sp[orthogroups,]), outline=F, notch=T, las=2, xaxt = "n", ylab="log(TPM)")

par(mai=c(1,0.8,0,0.3), mgp=c(2.2,0.8,0))
boxplot(t(score.sp[orthogroups,]), outline=F, notch=T, ylab="Normalised rank score", las=1, main="", xaxt = "n")
axis(1, at = 1:length(a$names), labels = NA, cex.axis = 1, srt=45, tck = -0.02)
text(x = 1:length(a$names), y = par("usr")[3] - 5, labels = a$names, srt = 45, adj = 1, xpd = TRUE, cex=0.8)

# Median all samples

plot(median_expression~median_score, cex=0.1, las=1, cex.lab=1.2,
     xlab="median Normalised rank score", ylab="median log(TPM)", main=paste("All", tissue, "samples of", sp))
test <- cor.test(median_expression, median_score)
mtext(paste("R =",signif(test$estimate,3)), line=-2, cex=1.2)

# SD all samples
plot(sd_expression~sd_score, cex=0.1, las=1, cex.lab=1.2,
     xlab="Standard Deviation Score", ylab="Standard Deviation log(TPM)", main=paste("All", tissue, "samples of", sp))

test <- cor.test(sd_expression, sd_score)
mtext(paste("R =",signif(test$estimate,3)), line=-2, cex=1.2)


plot(sd_expression~median_expression, cex=0.1, las=1, cex.lab=1.2,
     xlab="log(TPM)", ylab="Standard Deviation log(TPM)", main=paste("All", tissue, "samples of", sp))

plot(sd_score~median_score, cex=0.1, las=1, cex.lab=1.2,
     xlab="Score", ylab="Standard Deviation Score", main=paste("All", tissue, "samples of", sp))

################################################################################

sp2 = "Pan_paniscus"
expression.sp2 = t(express[which(express$species == sp2),3:nb_ortho])
score.sp2 = t(score[which(score$species == sp2),3:nb_ortho])
colnames(expression.sp2) <- express[which(express$species == sp2),1]
colnames(score.sp2) <- express[which(express$species == sp2),1]

median_expression2 <- apply(expression.sp2, 1, function(x) median(x, na.rm=T))
median_score2 <- apply(score.sp2, 1, function(x) median(x, na.rm=T))

common.ortho <- intersect(names(median_expression), names(median_expression2))

test <- cor.test(median_expression[common.ortho], median_expression2[common.ortho])
plot(median_expression[common.ortho], median_expression2[common.ortho], cex=0.1, las=1, cex.lab=1.2, 
     xlab=paste(sp, "log(TPM)"), ylab=paste(sp2, "log(TPM)"), main=paste("R =",signif(test$estimate,3)))

test <- cor.test(median_score[common.ortho], median_score2[common.ortho])
plot(median_score[common.ortho], median_score2[common.ortho], cex=0.1, las=1, cex.lab=1.2,
     xlab=paste(sp, "Score"), ylab=paste(sp2, "Score"), main=paste("R =",signif(test$estimate,3)))


med_low <- median_expression[common.ortho]
med_low <- med_low[med_low<8]
ortho_low <- names(med_low)

# minus high expressed
test <- cor.test(median_expression[ortho_low], median_expression2[ortho_low])
plot(median_expression[ortho_low], median_expression2[ortho_low], cex=0.1, las=1, cex.lab=1.2, 
     xlab=paste(sp, "log(TPM)"), ylab=paste(sp2, "log(TPM)"), main=paste("R =",signif(test$estimate,3)))

test <- cor.test(median_score[ortho_low], median_score2[ortho_low])
plot(median_score[ortho_low], median_score2[ortho_low], cex=0.1, las=1, cex.lab=1.2,
     xlab=paste(sp, "Score"), ylab=paste(sp2, "Score"), main=paste("R =",signif(test$estimate,3)))


