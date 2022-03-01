plotCont <- function(data.A, data.B, xlab.A, xlab.B=NULL, w.A=NULL, w.B=NULL,
                     type="density", ref=FALSE){
###################################################
    hist.bks <- function(x, w=NULL, n=NULL){
        
        xx <- x[!is.na(x)]
        
        if(is.null(w)) {
            iqr.x <- IQR(xx)
            deff.uw <- 1
        }
        else{
            ww <- w[!is.na(x)]
            wQ <- Hmisc::wtd.quantile(x = xx, weights = ww, probs = c(0.25, 0.75))
            iqr.x <- c(wQ[2] - wQ[1])
            deff.uw <- length(ww) * sum(ww^2)/(sum(ww)^2)
        }
        if(is.null(n)) nn <- length(xx) / deff.uw
        else nn <- n
        wid.FD <- 2 * iqr.x* 1/(nn^(1/3))
        bins <- (max(xx)-min(xx))/wid.FD
        bins <- floor(bins)+1
        if(bins==1){
            bins <- 2
            wid.FD <- (max(xx)-min(xx))/bins
        }
        seq(from=min(xx), to=max(xx), by=wid.FD)
    }
#################################################################    
# preparation    
    nA <- nrow(data.A)
    nB <- nrow(data.B)
    n <- min(nA, nB)
    if(is.null(xlab.B)) xlab.B <- xlab.A
    
    if(xlab.A==xlab.B) pr.lab <- xlab.A
    else pr.lab <- "x"
    
    xA <- data.A[, xlab.A]
    xB <- data.B[, xlab.B]
    # note: weights are scaled to sum to n
    if(!is.null(w.A)) ww.A <- data.A[, w.A]/sum(data.A[, w.A]) * nA
    else ww.A <- rep(1, nA)
    deff.uwA <- length(ww.A) * sum(ww.A^2)/(sum(ww.A)^2)
    
    if(!is.null(w.B)) ww.B <- data.B[, w.B]/sum(data.B[, w.B]) * nB
    else ww.B <- rep(1, nB)
    deff.uwB <- length(ww.B) * sum(ww.B^2)/(sum(ww.B)^2)
    
    ############################################
    # histograms & density plot comparing two distributions
    if(type=="hist" | type=="density"){
        if(ref) {
            bks <- hist.bks(x=xB, w=ww.B, n=min(nA/deff.uwA, nB/deff.uwB))
            if(bks[1] > min(xA)) bks[1] <- min(xA)
            if(bks[length(bks)] < max(xA)) bks[length(bks)] <- max(xA)
        }
        else    bks <- hist.bks(x=c(xA, xB), w=c(ww.A, ww.B), min(nA/deff.uwA, nB/deff.uwB))
         #cat(bks, fill=T)
        
        midp <- 0.5 * (bks[1:(length(bks)-1)] + bks[2:length(bks)])
        #cat(midp, fill=T)
        cxA <- cut(x = xA, breaks = bks, include.lowest = TRUE)
        cxB <- cut(x = xB, breaks = bks, include.lowest = TRUE)
        # tables
        tA <- prop.table(table(cxA))
        if(!is.null(w.A)) tA <- prop.table(xtabs(data.A[,w.A]~cxA))
        
        tB <- prop.table(table(cxB))
        if(!is.null(w.B)) tB <- prop.table(xtabs(data.B[,w.B]~cxB))
        
        if(type=="hist"){
            Var1 <- Freq <- NULL
            dfA <- cbind(data.frame(tA), sample="A")
            dfB <- cbind(data.frame(tB), sample="B")
            colnames(dfA) <- colnames(dfB) <- c("Var1", "Freq", "sample")
            df <- rbind(dfA, dfB)
            
            colnames(df) <- c("Var1", "Freq", "sample")
            
            # calculates tot. var dist.
            tvd <- 100*(1/2)*sum(abs(tA-tB))
            labtvd <- paste("tvd ", 
                            paste(round(tvd, 2), "%", sep=" "), 
                            sep="= ")
            btmlab <- paste(pr.lab, labtvd, sep=", ")
            
            out <- ggplot(data = df, 
                          aes(x = Var1, y = Freq, fill = sample)) +
                geom_bar(stat = "identity", position = "dodge") +
                labs(x = btmlab, y="rel freq")
        
        }
        if(type=="density"){
            x <- w <- NULL
            dfxA = data.frame(x=midp, sample="A")
            dfxB = data.frame(x=midp, sample="B")
            xx <- rbind(dfxA, dfxB)
            colnames(xx) <- c(pr.lab, "sample")
            xx$w <- c(tA, tB)
            out <- ggplot2::ggplot(xx, aes(xx[ ,pr.lab], weight=w, fill=sample, colour=sample)) +
                geom_density(alpha=0.4, lwd=0.8, adjust=0.5) +
                labs(x = pr.lab)
        }    
        
    }

    ###########################################
    
    if(type=="ecdf"){
        if(ref) usx <- unique(sort(xB))
        else usx <- unique(sort(c(xA ,xB)))
        
        k <- length(usx)
        ecdf.xA <- ecdf.xB <- numeric(k)
        for(i in 1:k){
            ecdf.xA[i] <- sum(ww.A[xA <= usx[i]])/sum(ww.A)
            ecdf.xB[i] <- sum(ww.B[xB <= usx[i]])/sum(ww.B)
        }
        
        dfxA = data.frame(x=usx, ecdf=ecdf.xA, sample="A")
        dfxB = data.frame(x=usx, ecdf=ecdf.xB, sample="B")
        xx <- rbind(dfxA, dfxB)
        
        out <- ggplot2::ggplot(xx, aes(x=x, y=ecdf, color=sample)) +
            geom_point(size=0.5) +
            geom_line() +
            xlab(label = xlab.A) +
            ylab(label = "Empirical cumulative density function")
            
    }
    ##############################################
    # comparison of empirical quantiles 
    if(type=="qqplot" | type=="qqshift"){
        # decision concerning quantiles
        if(n<20) pctp <- c(0.25, 0.50, 0.75)
        if(n>=20 & n<=30) pctp <- seq(from = 0.1,to = 0.9,by = 0.1)
        if(n>30) pctp <- seq(from = 0.05,to = 0.95,by = 0.05)
        # derive empirical quantiles
        qA <- quantile(x = xA, probs = pctp)
        if(!is.null(w.A)) {
            qA <- Hmisc::wtd.quantile(x=xA, 
                                      weights=data.A[, w.A], 
                                      probs=pctp)
        }
        qB <- quantile(x = xB, probs = pctp)
        if(!is.null(w.B)) {
            qB <- Hmisc::wtd.quantile(x=xB, 
                                      weights=data.B[, w.B], 
                                      probs=pctp)
        }
        if(type=="qqplot"){
            qq <- data.frame(qA=qA, qB=qB)
            #colnames(qq) <- paste(xlab, c("A","B"), sep=".")
            out <- ggplot2::ggplot(qq, aes(x=qB, y=qA)) +
                geom_point(shape=19, color="orangered1", size=3) +
                geom_abline(intercept = 0, slope = 1, color="blue", size=1.5, 
                        linetype="dashed") +
                xlab(label = paste("Quantiles", pr.lab, "B", sep=" ")) +
                ylab(label = paste("Quantiles", pr.lab, "A", sep=" ")) 
        }
        #############################################################
        # shift between two empirical distributions
        # i.e. scatterplot of difference between quantiles (Y)
        # and quantiles in A (X)
        if(type=="qqshift"){
            
            qq <- data.frame(qA=qB, qB=qA-qB)
            #colnames(qq) <- paste(xlab, c("A","B"), sep=".")
            out <- ggplot(qq, aes(x=qA, y=qB)) +
                geom_point(shape=19, color="orangered1", size=3) +
                geom_hline(yintercept = 0, color="blue", size=1.5, 
                           linetype="dashed") +
                xlab(label = paste("Quantiles", pr.lab, "B", sep=" ")) +
                ylab(label = "Quantiles.A - Quantiles.B") 
        }
        
    }
 ###################################################################
#####################################################################
 
    plot(out)

}