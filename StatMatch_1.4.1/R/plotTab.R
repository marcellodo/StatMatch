plotTab <- function(data.A, data.B, xlab.A, xlab.B=NULL, w.A=NULL, w.B=NULL){
    
    x <- Freq <- Sample <- NULL
    if(is.null(xlab.B)) xlab.B <- xlab.A
    if(!is.null(xlab.B)) if(length(xlab.B) != length(xlab.A)) stop("Different number of variables")
    
    p <- length(xlab.A)
    
    fxA <- paste(xlab.A, collapse="+")
    fxB <- paste(xlab.B, collapse="+")
    
    if(is.null(w.A)) ffxA <- paste0("~", fxA)
    else ffxA <- paste(w.A, fxA, sep="~")
    
    if(is.null(w.B)) ffxB <- paste0("~", fxB)
    else ffxB <- paste(w.B, fxB, sep="~")
    
    # table on A
    tA <- prop.table(xtabs(as.formula(ffxA), data=data.A))

    # table on B
    tB <- prop.table(xtabs(as.formula(ffxB), data=data.B))

    tvd <- 100*(1/2)*sum(abs(tA-tB))
    labtvd <- paste("tvd ", 
                    paste(round(tvd, 2), "%", sep=" "), 
                    sep="= ")
    
    
    dfA <- cbind(data.frame(tA), sample="A")
    dfB <- cbind(data.frame(tB), sample="B")
    
    if(all(xlab.B == xlab.A)) xlab <- xlab.A
    else xlab <- paste0("x", 1:p)
    
    colnames(dfA) <- colnames(dfB) <- c(xlab, "Freq", "Sample")
    df <- rbind(dfA, dfB)
    if(p>1) df1 <- apply(df[,1:p], 1, paste, collapse="*")
    else df1 <- df[,1]
    
    newdf <- cbind(x=df1, df[, c("Freq", "Sample")])
    
    btmlab <- paste(paste(xlab, collapse="*"), labtvd, sep=", ")
    
    ggplot2::ggplot(data = newdf, aes(x = x, y = Freq, fill = Sample)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(x = btmlab, y="Rel. freq.")
}
