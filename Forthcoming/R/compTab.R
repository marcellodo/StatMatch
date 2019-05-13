compTab <- function(data.A, data.B, xlab, wA=NULL, wB=wA){
    
    # table on A
    if(is.null(wA)) tA <- prop.table(table(data.A[,xlab]))
    else {
        fxw <- paste(wA, xlab, sep="~")
        tA <- prop.table(xtabs(as.formula(fxw), data=data.A))
    }
    
    # table on B
    if(is.null(wB)) tB <- prop.table(table(data.B[,xlab]))
    else {
        fxw <- paste(wB, xlab, sep="~")
        tB <- prop.table(xtabs(as.formula(fxw), data=data.B))
    }    
    tvd <- 100*(1/2)*sum(abs(tA-tB))
    labtvd <- paste("tvd ", 
                    paste(round(tvd, 2), "%", sep=" "), 
                    sep="= ")
    btmlab <- paste(xlab, labtvd, sep=", ")
    
    dfA <- cbind(data.frame(tA), sample="A")
    dfB <- cbind(data.frame(tB), sample="B")
    colnames(dfA) <- colnames(dfB) <- c("Var1", "Freq", "sample")
    df <- rbind(dfA, dfB)
    ggplot(data = df, 
           aes(x = Var1, y = Freq, fill = sample)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(x = btmlab, y="rel freq")
}
