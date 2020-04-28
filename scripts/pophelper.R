library("pophelper")


meta_filename <- snakemake@input[[1]]
#print(meta_filename)
fam_filename <- snakemake@input[[2]]
#print(fam_filename)
cv_filename <- snakemake@input[[3]]
#print(cv_filename)
title_filename <- snakemake@input[[4]]
#print(title_filename)
q_filenames <- rep("",25)
for (n in 1:25){
	q_filenames[n] <- snakemake@input[[n+4]]
}
#print(q_filenames)
out_path <- snakemake@params[[1]]
#print(out_path)
out_filename <- snakemake@params[[2]]
#print(out_filename)
out_meta <- snakemake@output[[2]]
#print(out_meta)

# Read in CV numbers
cv_values <- read.table(cv_filename,sep=" ")[10:25,4]
cv_value_labs <- paste(paste("K",10:25,sep=""),cv_values)
#print(cv_value_labs)

# Read in the admixture output files
slist <- readQ(q_filenames,filetype="auto")

fam <- read.table(fam_filename, sep=" ",stringsAsFactors=F)
#print(fam)
sample_anno <- read.table(meta_filename, sep="\t",header=TRUE,stringsAsFactors=F)
#print(sample_anno)

# Select from Metafile the annotation in the correct order
meta <- sample_anno[sample_anno[,1] %in% fam[,1],]
meta_ordered <- meta[match(fam[,1],meta$SAMPLE),]
# Check whether everything is in correct order
all(meta_ordered$SAMPLE == fam[,1])
groups <- meta_ordered[,2:6]
# Maker sure these are characters
sapply(groups, is.character)
#print(groups)
write.table(meta_ordered,file=out_meta,sep="\t",quote=FALSE,row.names=FALSE)

# Read file with number of variants, make title
numvariants_info <- readLines(title_filename)
title_numbers <- paste(length(meta$SAMPLE)," individuals, ",numvariants_info[6]," variants (",numvariants_info[4]," filtered, ",numvariants_info[2]," intersected)",sep="")
#print(title_numbers)

# Colors "shiny" and "kelly_22 (1:5) plus black for having 25 colors overall
#shiny_colors <- c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E","#000000","#F2F3F4","#222222","#F3C300","#875692","#F38400")
# Now, we use the paper colors, as specified by Michael (plus black for 25th component)
#shiny_colors <- c("#9467BD","#F7B6D2","#7F7F7F","#9EDAE5","#98DF8A","#C49C94","#2CA02C","#C7C7C7","#FF7F0E","#C5B0D5","#BCBD22","#F8DC0B","#D62728","#AEC7E8","#2CA02C","#17BECF","#E377C2","#000000","#FF9896","#FF7F0E","#8C564B","#1F77B4","#1F77B4","#DBDB8D","#FFBB78","#000000")
shiny_colors <- c("#9467BD","#F7B6D2","#7F7F7F","#9EDAE5","#98DF8A","#C49C94","#32CD32","#C7C7C7","#DAA520","#C5B0D5","#BCBD22","#F8DC0B","#D62728","#AEC7E8","#2CA02C","#17BECF","#E377C2","#000000","#FF9896","#FF7F0E","#8C564B","#FFFFF0","#1F77B4","#DBDB8D","#FFBB78")


plotQ(slist[10:25],imgoutput="join",grplab=groups,ordergrp=T,sharedindlab=FALSE,showlegend=T,showsp=FALSE,sortind="all",imgtype="pdf",exportpath=out_path,outputfilename=out_filename,splab=as.character(cv_value_labs),clustercol=shiny_colors,showtitle=TRUE,titlelab=title_numbers)

