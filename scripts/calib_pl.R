require(ape)
args <- commandArgs(trailingOnly = T)
in_phy <- as.character(args[1])
out_chro <- as.character(args[2])

phy <- read.tree(in_phy)
### correlated rate model with penalized likelihood smoothing param = 1:
chr <- chronos(phy)

write.tree(chr, out_chro)