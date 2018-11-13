#functions

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
    p <- ggplot(df, aes(x=df$PCoAx, y=df$PCoAy, group=ind, colour=pop)) +
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
  
  # If population labels
  
  if (labels == "pop") {
    cat("Plotting populations\n")
    ind <- indNames(data)
    pop <- factor(pop(data))
    df <- cbind(df,ind,pop)
    colnames(df) <- c("PCoAx","PCoAy","ind","pop")
    
    # Plot
    p <- ggplot(df, aes(x=df$PCoAx, y=df$PCoAy)) +
      geom_point(size=3,aes(colour=type, shape=pop, group=type)) +
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
    cat("Preparing a plot for interactive labelling, follow with ggplotly()\n")
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