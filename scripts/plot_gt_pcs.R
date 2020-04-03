# Writes the first 10 genotype principal components to file for later inclusion 
# into the linear model. Computes various visualizations, e.g. PCs, heatmaps and
# scree plot.

library("latticeExtra")
library("scatterplot3d")

# Get the name of the input file
filename <- snakemake@input[[1]]
annotation_file <- snakemake@input[[2]]
path <- snakemake@params[[1]]
anno_param <- as.integer(snakemake@params[[2]])
fname_pca_1vs2 <- snakemake@output[[1]]
fname_pca_1vs3 <- snakemake@output[[2]]
fname_pca_1vs4 <- snakemake@output[[3]]
fname_pca_2vs3 <- snakemake@output[[4]]
fname_pca_2vs4 <- snakemake@output[[5]]
fname_pca_3vs4 <- snakemake@output[[6]]
fname_scree <- snakemake@output[[7]]
fname_pca_3d <- snakemake@output[[8]]

########## Get Eigenstrat PCs ##########

# Lese die Eigenvalues von der ersten Zeile (Kommentarzeile)
eigenvalues <- as.vector(read.table(filename,nrows=1,comment.char = "", row.names=1))
print(dim(eigenvalues))
pcs <- read.table(filename,row.names=1,col.names=c("rowname","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","pheno"))


########## Sample annotation ##########

# Obtain the sample annotation
sample_annotation <- read.table(annotation_file, header=FALSE, sep="\t")

# Bring samples in genotype and in sample annotation in the same order
order <- sample_annotation[,1]
pcs <- pcs[as.character(order),]
# Make sure that rownames and sample annotation match
all(as.character(sample_annotation[,1]) == row.names(pcs))

########## Genotype principal components ##########

# Plot the PCs
eig_pc1 <- eigenvalues[1]
eig_pc2 <- eigenvalues[2]
eig_pc3 <- eigenvalues[3]
eig_pc4 <- eigenvalues[4]

# Rewrite empty info to "Unknown"
group <- as.character(sample_annotation[,anno_param])
print(group)
group[group==""] <- "Unknown"
group <- as.factor(group)
pc_col <- rainbow(length(levels(group)))

# Get the PC values of the Egyptians
#egyptian_ids <- sample_annotation[sample_annotation[,6]=="Egypt",]$V1
pc_egypt <- pcs[sample_annotation[,6]=="Egypt",]
head(pc_egypt)

# PC1/PC2
pdf(fname_pca_1vs2)
# pch=16
fig <- xyplot(PC2~PC1,groups=group, data=as.data.frame(pcs),pch=16,cex=0.3,col=pc_col,main=draw.key(key=list(rect=list(col=pc_col),text=list(levels(group)))),xlab="PC1",ylab="PC2")
## Overlay Egyptian samples
fig <- fig + layer(panel.points(pc_egypt$PC1, pc_egypt$PC2, pch=1, cex=1, col="black"))
fig
dev.off()

# PC1/PC3
pdf(fname_pca_1vs3)
fig <- xyplot(PC3~PC1,groups=group, data=as.data.frame(pcs),pch=16,cex=0.3,col=pc_col,main=draw.key(key=list(rect=list(col=pc_col),text=list(levels(group)))),xlab="PC1",ylab="PC3")
## Overlay Egyptian samples
fig <- fig + layer(panel.points(pc_egypt$PC1, pc_egypt$PC3, pch=21, cex=1, col="black"))
fig
dev.off()

# PC1/PC4
pdf(fname_pca_1vs4)
fig <- xyplot(PC4~PC1,groups=group, data=as.data.frame(pcs),pch=16,cex=0.3,col=pc_col,main=draw.key(key=list(rect=list(col=pc_col),text=list(levels(group)))),xlab="PC1",ylab="PC4")
## Overlay Egyptian samples
fig <- fig + layer(panel.points(pc_egypt$PC1, pc_egypt$PC4, pch=21, cex=1, col="black"))
fig
dev.off()

# PC2/PC3
pdf(fname_pca_2vs3)
fig <- xyplot(PC3~PC2,groups=group, data=as.data.frame(pcs),pch=16,cex=0.3,col=pc_col,main=draw.key(key=list(rect=list(col=pc_col),text=list(levels(group)))),xlab="PC2",ylab="PC3")
## Overlay Egyptian samples
fig <- fig + layer(panel.points(pc_egypt$PC2, pc_egypt$PC3, pch=21, cex=1, col="black"))
fig
dev.off()

# PC2/PC4
pdf(fname_pca_2vs4)
fig <- xyplot(PC4~PC2,groups=group, data=as.data.frame(pcs),pch=16,cex=0.3,col=pc_col,main=draw.key(key=list(rect=list(col=pc_col),text=list(levels(group)))),xlab="PC2",ylab="PC4")
## Overlay Egyptian samples
fig <- fig + layer(panel.points(pc_egypt$PC2, pc_egypt$PC4, pch=21, cex=1, col="black"))
fig
dev.off()

# PC3/PC4
pdf(fname_pca_3vs4)
fig <- xyplot(PC4~PC3,groups=group, data=as.data.frame(pcs),pch=16,cex=0.3,col=pc_col,main=draw.key(key=list(rect=list(col=pc_col),text=list(levels(group)))),xlab="PC3",ylab="PC4")
## Overlay Egyptian samples
fig <- fig + layer(panel.points(pc_egypt$PC3, pc_egypt$PC4, pch=21, cex=1, col="black"))
fig
dev.off()

# Scree plot
pdf(fname_scree)
plot(1:20,eigenvalues,type="b",main="Scree plot genotype PCA")
dev.off()

# 3D Plot

# Make a vector of colors for the dots (samples) to be plotted
num_samples <- length(group)
sample_names <- sample_annotation[,1]
cols <- rep("",num_samples)
for (i in 1:num_samples){
  sample_class <- group[i]
  cols[i] <- pc_col[match(sample_class,as.vector(levels(group)))]
}

pdf(fname_pca_3d,useDingbats=FALSE)
s3d <- scatterplot3d(pcs$PC1, pcs$PC2, pcs$PC3,
              color='black', pch=21, box=FALSE, tick.marks=FALSE, bg=cols,
              xlab="PC1",
              ylab="PC2",
              zlab="PC3")
legend("topleft",inset=.05,bty="n",title="",levels(group),fill=pc_col,cex=0.8)
dev.off()
