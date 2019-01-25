#functions
#from https://github.com/green-striped-gecko/dartR.git
#https://github.com/green-striped-gecko/dartR/blob/master/R/gl.pcoa.plot.r
gl.pcoa.plot <- function(glPca, data, type, title, scale=FALSE, ellipse=FALSE, p=0.95, labels="pop", hadjust=1.5, 
                         vadjust=1, xaxis=1, yaxis=2) {
  
  if(class(glPca)!="glPca" | class(data)!="genlight") {
    cat("Fatal Error: glPca and genlight objects required for glPca and data parameters respectively!\n"); stop()
  }
  
  # Tidy up the parameters
  #  if (labels=="smart") { hadjust <- 0; vadjust <- 0 }
  
  # Create a dataframe to hold the required scores
  m <- cbind(glPca$scores[,xaxis],glPca$scores[,yaxis])
  df <- data.frame(m)
  
  # Convert the eigenvalues to percentages
  s <- sum(glPca$eig)
  e <- round(glPca$eig*100/s,1)
  
  # Labels for the axes
  xlab <- paste("PCoA Axis", xaxis, "(",e[xaxis],"%)")
  ylab <- paste("PCoA Axis", yaxis, "(",e[yaxis],"%)")
  
  # If individual labels
  
  if (labels == "ind") {
    cat("Plotting individuals\n")
    ind <- indNames(data)
    pop <- factor(pop(data))
    df <- cbind(df,ind,pop)
    colnames(df) <- c("PCoAx","PCoAy","ind","pop")
    
    # Plot
    p <- ggplot(df, aes(x=df$PCoAx, y=df$PCoAy, group=pop, colour=ind)) +
      geom_point(size=1,aes(colour=ind)) +
      #geom_dl(aes(label=ind),method="first.points") +
      #ggtitle(paste("PCoA Plot")) +
      theme(axis.title=element_text(face="bold.italic",size="20", color="black"),
            axis.text.x  = element_text(face="bold",angle=0, vjust=0.5, size=10),
            axis.text.y  = element_text(face="bold",angle=0, vjust=0.5, size=10),
            legend.title = element_text(colour="black", size=18, face="bold"),
            legend.text = element_text(colour="black", size = 16, face="bold")
      ) +
      labs(x=xlab, y=ylab) +
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0)
    # Scale the axes in proportion to % explained, if requested
    if(scale==TRUE) { p <- p + coord_fixed(ratio=e[yaxis]/e[xaxis]) }
    # Add ellipses if requested
    if(ellipse==TRUE) {p <- p + stat_ellipse(aes(colour=pop), type="norm", level=0.95)}
  } 
  
  # If population labels
  
  if (labels == "pop") {
    #cat("Plotting populations\n")
    ind <- indNames(data)
    pop <- factor(pop(data))
    df <- cbind(df,ind,pop)
    colnames(df) <- c("PCoAx","PCoAy","ind","pop")
    
    # Plot
    p <- ggplot(df, aes(x=df$PCoAx, y=df$PCoAy)) +
      geom_point(size=1.5,aes(colour=pop, shape=type)) +
      ggtitle(paste("PCoA Plot", title, sep=" ")) +
      theme(axis.title=element_text(face="bold.italic",size="10", color="black"),
            axis.text.x  = element_text(face="bold",angle=0, vjust=0.5, size=10),
            axis.text.y  = element_text(face="bold",angle=0, vjust=0.5, size=10),
            legend.title = element_text(colour="black", size=10, face="bold"),
            legend.text = element_text(colour="black", size=10, face="bold"),
            title = element_text(colour="black", size=10, face="bold")
      ) +
      labs(x=xlab, y=ylab) +
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0) +
      theme(legend.position="right")+
      scale_color_manual(values=c("red","blue"))
    # Scale the axes in proportion to % explained, if requested
    if(scale==TRUE) { p <- p + coord_fixed(ratio=e[yaxis]/e[xaxis]) }
    # Add ellipses if requested
    if(ellipse==TRUE) {p <- p + stat_ellipse(aes(colour="black"), type="norm", level=0.95)}
  }
  
  # If interactive labels
  
  if (labels=="interactive" | labels=="ggplotly") {
    #cat("Preparing a plot for interactive labelling, follow with ggplotly()\n")
    ind <- as.character(indNames(data))
    pop <- as.character(pop(data))
    df <- cbind(df,pop,ind)
    colnames(df) <- c("PCoAx","PCoAy","pop","ind")
    #df$ind <- as.character(df$ind)
    #df$pop <- as.character(df$pop)
    x <- df$PCoAx
    y <- df$PCoAy
    
    # Plot
    p <- ggplot(df, aes(x=x, y=y)) +
      geom_point(size=2,aes(colour=pop, fill=ind)) +
      #geom_dl(aes(label=pop),method="smart.grid") +
      #ggtitle(paste("PCoA Plot")) +
      theme(axis.title=element_text(face="bold.italic",size="20", color="black"),
            axis.text.x  = element_text(face="bold",angle=0, vjust=0.5, size=10),
            axis.text.y  = element_text(face="bold",angle=0, vjust=0.5, size=10),
            legend.title = element_text(colour="black", size=18, face="bold"),
            legend.text = element_text(colour="black", size = 16, face="bold")
      ) +
      labs(x=xlab, y=ylab) +
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0) +
      theme(legend.position="none")
    # Scale the axes in proportion to % explained, if requested
    if(scale==TRUE) { p <- p + coord_fixed(ratio=e[yaxis]/e[xaxis]) }
    # Add ellipses if requested
    if(ellipse==TRUE) {p <- p + stat_ellipse(aes(colour=pop), type="norm", level=0.95)}
    cat("Ignore any warning on the number of shape categories\n")
  }  
  
  # If labels = legend
  
  if (labels == "legend") {
    cat("Plotting populations identified by a legend\n")
    pop <- factor(pop(data))
    df <- cbind(df,pop)
    colnames(df) <- c("PCoAx","PCoAy","pop")
    
    # Plot
    p <- ggplot(df, aes(x=df$PCoAx, y=df$PCoAy,colour=pop)) +
      geom_point(size=2,aes(colour=pop)) +
      #geom_dl(aes(label=ind),method="first.points") +
      #ggtitle(paste("PCoA Plot")) +
      theme(axis.title=element_text(face="bold.italic",size="20", color="black"),
            axis.text.x  = element_text(face="bold",angle=0, vjust=0.5, size=10),
            axis.text.y  = element_text(face="bold",angle=0, vjust=0.5, size=10),
            legend.title = element_text(colour="black", size=18, face="bold"),
            legend.text = element_text(colour="black", size = 16, face="bold")
      ) +
      labs(x=xlab, y=ylab) +
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0)
    # Scale the axes in proportion to % explained, if requested
    if(scale==TRUE) { p <- p + coord_fixed(ratio=e[yaxis]/e[xaxis]) }
    # Add ellipses if requested
    if(ellipse==TRUE) {p <- p + stat_ellipse(aes(colour=pop), type="norm", level=0.95)}
  } 
  
  # If labels = none
  
  if (labels == "none" | labels==FALSE) {
    cat("Plotting points with no labels\n")
    pop <- factor(pop(data))
    df <- cbind(df,pop)
    colnames(df) <- c("PCoAx","PCoAy","pop")
    
    # Plot
    p <- ggplot(df, aes(x=df$PCoAx, y=df$PCoAy,colour=pop)) +
      geom_point(size=2,aes(colour=pop)) +
      #geom_dl(aes(label=ind),method="first.points") +
      #ggtitle(paste("PCoA Plot")) +
      theme(axis.title=element_text(face="bold.italic",size="20", color="black"),
            axis.text.x  = element_text(face="bold",angle=0, vjust=0.5, size=10),
            axis.text.y  = element_text(face="bold",angle=0, vjust=0.5, size=10),
            legend.title = element_text(colour="black", size=18, face="bold"),
            legend.text = element_text(colour="black", size = 16, face="bold")
      ) +
      labs(x=xlab, y=ylab) +
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0)+
      theme(legend.position="none")
    # Scale the axes in proportion to % explained, if requested
    if(scale==TRUE) { p <- p + coord_fixed(ratio=e[yaxis]/e[xaxis]) }
    # Add ellipses if requested
    if(ellipse==TRUE) {p <- p + stat_ellipse(aes(colour=pop), type="norm", level=0.95)}
  }
  
  # If interactive labels
  
  #if (labels=="interactive" | labels=="ggplotly") {
  #  ggplotly(p)
  #} else {
  p
  #}
  
  return (p)
}

gl.filter.monomorphs <- function (x, v=2, pb=FALSE) {
  
  # ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }
  
  if (v < 0 | v > 5){
    cat("    Warning: verbosity must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    v <- 2
  }
  
  if (v > 0) {
    cat("Starting gl.filter.monomorphs: Deleting monomorphic loci\n")
  }
  
  # Create a vector to hold test results
  a <- vector(mode="logical", length=nLoc(x))
  for (i in 1:nLoc(x)) {a[i] <- NA}
  # Set up the progress counter
  if (v > 1 && pb == TRUE) {
    progress <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
    getTxtProgressBar(progress)
  }
  # Identify polymorphic, monomorphic and 'all na' loci
  # Set a <- TRUE if monomorphic, or if all NAs
  xmat <-as.matrix(x)
  for (i in (1:nLoc(x))) {
    a[i] <- all(xmat[,i]==0,na.rm=TRUE) || all(xmat[,i]==2,na.rm=TRUE)
    if (all(is.na(xmat[,i]))) {a[i] <- NA}
    if (v > 1 && pb == TRUE) {setTxtProgressBar(progress, i/nLoc(x))}
  }
  # Count the number of monomorphic loci (TRUE), polymorphic loci (FALSE) and loci with no scores (all.na)
  counts <- plyr::count(a)
  if (v > 2) {
    if (pb) {cat("\n")}
    cat("  Polymorphic loci:", counts[1,2], "\n")
    cat("  Monomorphic loci:", counts[2,2], "\n")
    if (is.na(counts[3,2])) {counts[3,2] <- 0}
    cat("  Loci with no scores (all NA):" , counts[3,2] ,"\n")
  }
  
  #Treat all na loci as monomorphic
  a[is.na(a)] <- TRUE
  # Write the polymorphic loci to a new genlight object
  if (v > 1) {cat("  Deleting monomorphic loci and loci with all NA scores\n")}
  
  x <- x[,(a==FALSE)]
  #x@other$loc.metrics <- x@other$loc.metrics[(a==FALSE),]
  
  if (v > 0) {
    cat("Completed gl.filter.monomorphs\n\n")
  }
  
  return (x)
}

utils.recalc.freqhomsnp <- function(x, v=2) {
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for utils.recalc.freqhomsnp!\n"); stop("Execution terminated\n")
  }
  if (v > 0) {
    cat("Starting utils.recalc.freqhomref: Recalculating frequency of homozygotes, alternate allele\n")
  }
  if (is.null(x@other$loc.metrics$FreqHomSnp)) {
    x@other$loc.metrics$FreqHomSnp <- array(NA,nLoc(x))
    if (v >= 3){
      cat("  Locus metric FreqHomSnp does not exist, creating slot @other$loc.metrics$FreqHomSnp\n")
    }
  }
  
  # Do the deed
  t <- as.matrix(x)
  for (i in 1:nLoc(x)) {
    x@other$loc.metrics$FreqHomSnp[i] <- length(which(t[,i] == 2))/(nInd(x)-length(which(is.na(t[,i]))))
  }
  
  if (v > 0) {
    cat("Completed utils.recalc.freqhomsnp\n\n")
  }
  
  return(x)
}

utils.recalc.freqhomref <- function(x, v=2) {
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for utils.recalc.freqhomref!\n"); stop("Execution terminated\n")
  }
  if (v > 0) {
    cat("Starting utils.recalc.freqhomref: Recalculating frequency of homozygotes, reference allele\n")
  }
  if (is.null(x@other$loc.metrics$FreqHomRef)) {
    x@other$loc.metrics$FreqHomRef <- array(NA,nLoc(x))
    if (v >= 3){
      cat("  Locus metric FreqHomRef does not exist, creating slot @other$loc.metrics$FreqHomRef\n")
    }
  }  
  
  # Do the deed
  t <- as.matrix(x)
  for (i in 1:nLoc(x)) {
    x@other$loc.metrics$FreqHomRef[i] <- length(which(t[,i] == 0))/(nInd(x)-length(which(is.na(t[,i]))))
  }
  
  if (v > 0) {
    cat("Completed utils.recalc.freqhomref\n\n")
  }
  
  return(x)
}

utils.recalc.freqhets <- function(x, v=2) {
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for utils.recalc.freqhets!\n"); stop("Execution terminated\n")
  }
  if (v > 0) {
    cat("Starting utils.recalc.freqhets: Recalculating frequency of heterozygotes\n")
  }
  if (is.null(x@other$loc.metrics$FreqHets)) {
    x@other$loc.metrics$FreqHets <- array(NA,nLoc(x))
    if (v >= 3){
      cat("  Locus metric FreqHets does not exist, creating slot @other$loc.metrics$FreqHets\n")
    }
  }
  
  # Do the deed
  t <- as.matrix(x)
  for (i in 1:nLoc(x)) {
    x@other$loc.metrics$FreqHets[i] <- length(which(t[,i] == 1))/(nInd(x)-length(which(is.na(t[,i]))))
  }
  
  if (v > 0) {
    cat("Completed utils.recalc.freqhets\n\n")
  }
  
  return(x)
}

utils.recalc.maf <- function(x, v=2) {
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.drop.pop.r!\n"); stop("Execution terminated\n")
  }
  
  if (v > 0) {
    cat("Starting gl.report.maf: Minimum Allele Frequency\n")
  }
  
  if (v < 0 | v > 5){
    cat("Warning: Verbosity must take on an integer value between 0 and 5, set to 3\n")
    v <- 3
  }
  
  # Recalculate the relevant loc.metrics
  
  if (v >= 3) {cat("  Removing monomorphic loci and recalculating FreqHoms and FreqHets\n")}
  
  x <- gl.filter.monomorphs(x, v = v)
  x <- utils.recalc.freqhets(x,v=v)
  x <- utils.recalc.freqhomref(x,v=v)
  x <- utils.recalc.freqhomsnp(x,v=v)
  
  # Calculate and plot overall MAF
  
  if (v >= 3) {cat("Calculating MAF\n")}
  if (is.null(x@other$loc.metrics$maf)) {
    if (v >= 3){
      cat("  Locus metric maf does not exist, creating slot @other$loc.metrics$maf\n")
    }
    x@other$loc.metrics$maf <- array(NA,nLoc(x))
  } else {
    if (v >= 3){cat("  Recalculating  minor allele frequency\n")}
  }
  
  homref <- x@other$loc.metrics$FreqHomRef
  homalt <- x@other$loc.metrics$FreqHomSnp
  het <- x@other$loc.metrics$FreqHets
  
  for (i in 1:nLoc(x)){
    x@other$loc.metrics$maf[i] <- min((homref[i]*2 + het[i]), (homalt[i]*2 + het[i]))/2
  }
  
  if (v > 0) {
    cat("Completed gl.filter.maf\n\n")
  }
  
  return(x)
}  

gl.filter.maf <- function(x, threshold=0.01, v=2) {
  
  # ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.drop.pop.r!\n"); stop("Execution terminated\n")
  }
  
  if (v < 0 | v > 5){
    cat("Warning: Verbosity must take on an integer value between 0 and 5, set to 2\n")
    v <- 2
  }
  
  if (v > 0) {
    cat("Starting gl.report.maf: Minimum Allele Frequency\n")
  }
  
  if (threshold > 0.5 | threshold <= 0) {
    cat("Warning: threshold must be in the range (0,0.5], but usually small, set to 0.05\n")
    threshold <- 0.05
  }
  
  # Recalculate the relevant loc.metrics
  
  if (v >= 3) {cat("  Removing monomorphic loci and recalculating FreqHoms and FreqHets\n")}
  print(dim(x@other$loc.metrics$maf))
  x <- utils.recalc.maf(x,v=v)
  print(dim(x@other$loc.metrics$maf))
  
  # Remove loci with NA count <= 1-threshold
  index <- x@other$loc.metrics$maf >= threshold
  x2 <- x[ ,index]
  #x2@other$loc.metrics <- x@other$loc.metrics[index,]
  x2 <- utils.recalc.maf(x2,v=v)
  
  if (v > 2) {
    cat("  Initial number of loci:", nLoc(x), "\n")
    cat("    Number of loci deleted:", nLoc(x) - nLoc(x2), "\n")
    cat("  Final number of loci:", nLoc(x2), "\n")    
  }
  
  if (v > 0) {
    cat("Completed gl.filter.maf\n\n")
  }
  
  return(x2)
}

