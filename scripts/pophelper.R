library("pophelper")


meta_filename <- snakemake@input[[1]]
#print(meta_filename)
fam_filename <- snakemake@input[[2]]
#print(fam_filename)
cv_filename <- snakemake@input[[3]]
#print(cv_filename)
title_filename <- snakemake@input[[4]]
#print(title_filename)
q_filenames <- rep("",15)
for (n in 1:15){
	q_filenames[n] <- snakemake@input[[n+4]]
}
#print(q_filenames)
out_path <- snakemake@params[[1]]
#print(out_path)
out_filename <- snakemake@params[[2]]
#print(out_filename)

# Read in CV numbers
cv_values <- read.table(cv_filename,sep=" ")[5:15,4]
cv_value_labs <- paste(paste("K",5:15,sep=""),cv_values)
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

# Read file with number of variants, make title
numvariants_info <- readLines(title_filename)
title_numbers <- paste(length(meta$SAMPLE)," individuals, ",numvariants_info[6]," variants (",numvariants_info[4]," filtered, ",numvariants_info[2]," intersected)",sep="")
#print(title_numbers)

# Colors "shiny" plus black for having 20 colors overall
shiny_colors <- c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E")

plotQ(slist[5:15],imgoutput="join",grplab=groups,ordergrp=T,sharedindlab=FALSE,showlegend=T,showsp=FALSE,sortind="all",imgtype="pdf",exportpath=out_path,outputfilename=out_filename,splab=as.character(cv_value_labs),clustercol=shiny_colors,showtitle=TRUE,titlelab=title_numbers)

