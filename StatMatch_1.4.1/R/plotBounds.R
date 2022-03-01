plotBounds <- function(outFB){
    p <- length(outFB) 
    if(p>2){
        df <- data.frame(outFB[[1]])
        colnames(df) <- c("Y", "Z", "low.u")
        
        if(!is.null(outFB$low.cx)) df$low.cx <- c(outFB$low.cx)
        df$CIA <- c(outFB$CIA)
        if(is.null(outFB$CIA)) df$CIA <- c(outFB$IA)
        if(!is.null(outFB$up.cx)) df$up.cx <- c(outFB$up.cx)
        df$up.u <- c(outFB$up.u)
    }
    else df <- outFB$bounds
    
    n <- nrow(df)
    k <- ncol(df)
    
    if(k==5) {
        colnames(df) <- c("Y", "Z", "low.u", "CIA", "up.u")
    }
    
    #### plot
    mdf <- data.matrix(df)
    
    # determine position on the X axes
    eps <- tapply(mdf[,"up.u"]-mdf[,"low.u"], mdf[,2], max)
    aa <- round(max(eps),1) + 0.05
    aa <- aa*(1:max(mdf[,2]))
    xs <- rep(aa, each=max(mdf[,1]))
    xxs <- xs + mdf[, -(1:2)]
    
    # determine postion on the Y axes
    yys <- max(mdf[,1])+1-mdf[,1]
    
    # are of the plot
    par(mar=c(2.1, 7.1, 4.1, 1.1))
    
    # plot CIA
    plot(x = xxs[,"CIA"], y=yys, pch=18,
         xlab="",ylab="", xaxt="n", yaxt="n",axes = F,
         xlim = c(min(xxs)-0.1, max(xxs)+0.1),
         ylim = c(min(yys)-1, max(yys)+0.5)
         
    )
    
    # plot bounds
    for(i in 1:n){
        lines(x = c(xxs[i, "low.u"], xxs[i, "up.u"]), y=c(yys[i], yys[i]), lwd=1, lty=3)
        if(k>5) lines(x = c(xxs[i, "low.cx"], xxs[i, "up.cx"]), y=c(yys[i], yys[i]), lwd=2,lty=1)
    }
    
    # labels to put on the axes
    lab1 <- paste(names(df)[1], unique(df[,1]), sep=" = ")
    lab2 <- paste(names(df)[2], unique(df[,2]), sep=" = ")
    
    posCIA <- tapply(xxs[,"CIA"], df[,2], mean)
    axis(side = 3, at = posCIA, labels = lab2)
    axis(side = 2, at = unique(yys), labels = lab1, las=2)
    
    # add the CIA value and the bounds' width
    wdtu <- df$up.u - df$low.u
    if(k>5) {
        wdtcx <- df$up.cx - df$low.cx
        labwdt <- paste(round(wdtcx,3), " (", round(wdtu,3), ")", sep="")
    }
    else labwdt <- round(wdtu,3)
    text(xxs[,"CIA"], y=yys, labels = round(df$CIA,3), pos = 3, cex=0.8) 
    text(xxs[,"CIA"], y=yys, labels = labwdt, pos = 1, cex=0.8) 
    
    ### end
}