library("pophelper")

q_filename <- snakemake@input[[1]]
fam_filename <- snakemake@input[[2]]
meta_filename <- snakemake@input[[3]]
out_filename <- snakemake@output[[1]]
K <- snakemake@wildcards[["K"]]

#colors <- brewer.pal(K, "Set2")

qvals <- read.table(q_filename)
head(qvals)
fam <- read.table(fam_filename, sep=" ")
head(fam)
meta <- read.table(meta_filename, sep="\t",header=TRUE)
head(meta)

rownames(qvals) <- fam$V1
head(qvals)

pdf(out_filename)
#[order(qvals$V1),]
barplot(t(as.matrix(qvals)), col=rainbow(K),xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()
