roh_filename <- snakemake@input[[1]]
#print(roh_filename)
meta_filename <- snakemake@input[[2]]
#print(meta_filename)
meta_col <- as.integer(snakemake@params[[1]])
#print(meta_col)
out_filename <- snakemake@output[[1]]
#print(out_filename)

# PREPARATION -------------------------------------------------------------
source("scripts/_population_functions.R")

# IMPORT ------------------------------------------------------------------
EGY_ROH <- data.table::fread(input = roh_filename, header=TRUE)
#head(EGY_ROH)
# Extract the first column(sample) and the annoation specified annotation column
pop_info_all <- data.table::fread(input = meta_filename, sep="\t", header=FALSE)
pop_info <- pop_info_all[,c(1,meta_col),with=FALSE]
#head(pop_info) # always C1: Samples; C2: Group/Population

# PLOTTING ----------------------------------------------------------------
# transform input for plotting
pframe <- prepROHPlot(EGY_ROH, pop_info)
# plot to File and also return plot object for further use
p_1 <- plotROH(pframe, ptitle=out_filename, box.plot=T, p.value=FALSE)
