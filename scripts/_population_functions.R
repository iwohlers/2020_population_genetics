# ROH FUNCTIONS -----------------------------------------------------------
## Prepare Plink-ROH-output for plotting
# @param: ROH.list, Plink output - data.table::fread(input="EGY_ROH.hom", header=TRUE)
# @param: pop.info, Sample to Group mapping with two columns
# @return: classified ROHs long format - colnames: "SAMPLE", "class", "MB", "POPULATION"
prepROHPlot <- function(ROH.list, pop.info) {
  require(tibble)
  require(dplyr)
  ## 1. add 'MB' column to the list
  ROH.list$MB <- ROH.list[[grep("KB",colnames(ROH.list), ignore.case=TRUE)]] / 1000
  ## 2. add 'class' column to the list
  ROH.list$class <- factor(sapply(ROH.list$MB, function(value) ifelse(value <= 0.155, "short", ifelse(value <= 1.606, "medium", "long"))), levels=c("short","medium","long","all"))
  ## 3. equalize colnames in lists
  pop.info <- tibble::tibble("FID"=pop.info[[1]], "IID"=pop.info[[2]])
  ## get sum of ROHs per individual
  df <- dplyr::left_join(aggregate(MB ~ FID + class, ROH.list, sum), pop.info, by="FID")
  ## get sum of all classes
  tmp <- dplyr::left_join(data.frame("FID"=unique(df$FID), "class"="all", stringsAsFactors=FALSE), aggregate(MB ~ FID + IID, df, sum), by="FID")
  tmp <- tmp[,c("FID","class","MB","IID")]
  ## join frame for final output
  df <- rbind.data.frame(df, tmp)
  # rename before returning
  colnames(df) <- c("SAMPLE", "class", "MB", "POPULATION")
  
  return(df)
}

## Plot classified ROHs for populations - prints two *.pdf for classiificatiion and frequency
# @param: plot.df, output of prepROHPlot()
# @param: ptitle, part of output-filename
# @param: box.plot, TRUE to construct box-plots instead of violin-plots
# @param: p.value, TRUE to calculate p.values
# @return: list containing plots 'classification' and 'frequency'
plotROH <- function(plot.df, ptitle="untitled", box.plot=FALSE, p.value=FALSE) {
  library(ggpubr)
  # define things to compare - we'll go with groups aka 'IID' for now
  comp.select <- "POPULATION"
  comp.text <- "Populations"
  c.vector <- RColorBrewer::brewer.pal(length(unique(plot.df[[eval(comp.select)]])), "Dark2")
  my_comparisons <- combn(unique(plot.df[[eval(comp.select)]]), m=2, function(x) c(x), simplify = FALSE)
  if(box.plot) {
    p_ROH_all <- ggplot(plot.df, aes(x = plot.df[,eval(comp.select)], y = MB)) + 
      geom_boxplot(notch = TRUE) +
      #scale_y_continuous(trans='log2') + 
      geom_jitter(aes(color=plot.df[,eval(comp.select)]), alpha=0.25, position=position_jitter(0.2)) +
      scale_color_manual(values = c.vector) +
      facet_grid(rows=vars(class), scales="fixed", space="fixed") +
      # Remove x axis title
      theme_bw() +
      theme(axis.text.x=element_text(angle=90, hjust=1)) +
      theme(axis.title.x = element_blank()) + 
      theme(legend.position="bottom") +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      labs(colour = "Population", y = "Total length of ROHs (Mb)")
    #facet_wrap(~ class, ncol=2, nrow=2) +
  } else {
    ## 1. facet_grid plot all the groups
    p_ROH_all <- ggpubr::ggviolin(plot.df, x = eval(comp.select), y = "MB", trim = T,add.params = list(fill = "white"),
                                  fill = eval(comp.select), add = c("boxplot"), palette = c.vector, alpha = 0.65,
                                  ylab = "Total length of ROHs (Mb)", xlab = comp.text,
                                  title = "All ROH comparison") +
      #scale_y_continuous(trans='log2') +
      theme_bw() +
      theme(axis.text.x=element_text(angle=90, hjust=1)) +
      theme(axis.title.x = element_blank()) + 
      theme(legend.position="bottom") +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      labs(fill = comp.text) + 
      facet_grid(rows=vars(class), scales="fixed", space="fixed") 
    #facet_wrap(~ class, ncol=2, nrow=2) +
  }
  
  ## 2. add statistical comparison
  if(p.value) {
    p_ROH_all <- p_ROH_all + ggpubr::stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
      ggpubr::stat_compare_means() 
  }
  
  ### plot frequency
  pops <- unique(plot.df[,eval(comp.select)])
  popfreq.list <- list()
  for( pop.idx in 1:length(pops) ) {
    ## 1. select population to plot
    tmp <- plot.df[which(plot.df[,eval(comp.select)] %in% pops[pop.idx] & !(plot.df[,"class"] %in% "all")),]
    ## 2. get number of samples in population
    sample.scale <- length(unique(tmp$SAMPLE))
    
    popfreq.list[[eval(pops[pop.idx])]] <- ggplot(tmp, aes(MB, color=eval(comp.select))) +
      geom_histogram(aes(y=..count../ sample.scale), alpha=0.75, bins=100) +
      scale_color_manual(values = c.vector[pop.idx]) +
      theme_bw() +
      theme(axis.text.x=element_text(angle=0, hjust=1)) +
      theme(axis.title.x = element_blank()) +
      theme(axis.title.y = element_blank()) + 
      theme(legend.position = "none") + 
      labs(title=pops[pop.idx])
    
    
  }
  
  n <- length(popfreq.list)
  nCol <- floor(sqrt(n))
  pdf(file=paste(ptitle,sep=""), width=6 * nCol, height=length(unique(plot.df$class))*4 / nCol, onefile=FALSE)
  popfreq.plot <- grid.arrange(grobs=popfreq.list, ncol=nCol, left="ROH frequency scaled by number of samples")
  dev.off()
  
  pdf(file=paste("ROH_classification_",ptitle,".pdf",sep=""), width=6, height=length(unique(plot.df$class))*6, onefile=FALSE)
  print(p_ROH_all)
  dev.off()
  
  ## prepare return
  p.list <- list()
  p.list[["classification"]] <- p_ROH_all
  p.list[["frequency"]] <- popfreq.plot
  return(p.list)
}
