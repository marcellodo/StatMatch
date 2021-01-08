plotCont <- function(data.A, data.B, xlab.A, xlab.B=NULL, w.A=NULL, w.B=NULL,
                     type="density"){
################################################
    w <- Var1 <- Freq <- x <- Sample <- NULL
###############################    
    nA <- nrow(data.A)
    nB <- nrow(data.B)
    n <- min(nA, nB)
    if(is.null(xlab.B)) xlab.B <- xlab.A
    
    if(xlab.A==xlab.B) pr.lab <- xlab.A
    else pr.lab <- "x"
    
    ############################################
    # density plot comparing two distributions
    if(type=="density"){
        xA = data.frame(x=data.A[, xlab.A], sample="A")
        xB = data.frame(x=data.B[, xlab.B], sample="B")
        xx <- rbind(xA, xB)
        
        colnames(xx) <- c(pr.lab, "sample")

        if(!is.null(w.A)) ww.A <- data.A[, w.A]/sum(data.A[, w.A])
        else ww.A <- rep(1, nA)/nA
        
        if(!is.null(w.B)) ww.B <- data.B[, w.B]/sum(data.B[, w.B])
        else ww.B <- rep(1, nB)/nB
        
        xx$w <- c(ww.A, ww.B)
        out <- ggplot(xx, aes(xx[ ,pr.lab], weight=w, fill=sample, colour=sample)) +
               geom_density(alpha=0.4, lwd=0.8, adjust=0.5) +
               labs(x = pr.lab)
    }
    ##############################################
    # comparison of empirical quantiles 
    if(type=="qqplot"){
        # decision concerning quantiles
        if(n<20) pctp <- c(0.25, 0.50, 0.75)
        if(n>=20 & n<=30) pctp <- seq(from = 0.1,to = 0.9,by = 0.1)
        if(n>30) pctp <- seq(from = 0.05,to = 0.95,by = 0.05)
        # derive empirical quantiles
        qA <- quantile(x = data.A[, xlab.A], probs = pctp)
        if(!is.null(w.A)) {
            qA <- Hmisc::wtd.quantile(x=data.A[, xlab.A], 
                                      weights=data.A[, w.A], 
                                      probs=pctp)
        }
        qB <- quantile(x = data.B[, xlab.B], probs = pctp)
        if(!is.null(w.A)) {
            qB <- Hmisc::wtd.quantile(x=data.B[, xlab.B], 
                                      weights=data.B[, w.B], 
                                      probs=pctp)
        }
        
        qq <- data.frame(qA=qA, qB=qB)
        #colnames(qq) <- paste(xlab, c("A","B"), sep=".")
        out <- ggplot(qq, aes(x=qA, y=qB)) +
            geom_point(shape=19, color="orangered1", size=3) +
            geom_abline(intercept = 0, slope = 1, color="blue", size=1.5, 
                        linetype="dashed") +
            xlab(label = paste("Quantiles", pr.lab, "A", sep=" ")) +
            ylab(label = paste("Quantiles", pr.lab, "B", sep=" ")) 
        
    }
    #############################################################
    # shift between two empirical distributions
    # i.e. scatterplot of difference between quantiles (Y)
    # and quantiles in A (X)
    if(type=="qqshift"){
        if(n<20) pctp <- c(0.25, 0.50, 0.75)
        if(n>=20 & n<30) pctp <- seq(from = 0.1,to = 0.9,by = 0.1)
        if(n>30) pctp <- seq(from = 0.05,to = 0.95,by = 0.05)
        #
        qA <- quantile(x = data.A[, xlab.A], probs = pctp)
        if(!is.null(w.A)) {
            qA <- Hmisc::wtd.quantile(x=data.A[, xlab.A], 
                                      weights=data.A[, w.A], 
                                      probs=pctp)
        }
        qB <- quantile(x = data.B[, xlab.B], probs = pctp)
        if(!is.null(w.A)) {
                qB <- Hmisc::wtd.quantile(x=data.B[, xlab.B], 
                                      weights=data.B[, w.B], 
                                      probs=pctp)
        }
        
        qq <- data.frame(qA=qA, qB=qA-qB)
        #colnames(qq) <- paste(xlab, c("A","B"), sep=".")
        out <- ggplot(qq, aes(x=qA, y=qB)) +
            geom_point(shape=19, color="orangered1", size=3) +
            geom_hline(yintercept = 0, color="blue", size=1.5, 
                        linetype="dashed") +
            xlab(label = paste("Quantiles", pr.lab, "A", sep=" ")) +
            ylab(label = "Quantiles.A - Quantiles.B") 
    }
    ###################################################################
    # discretize continuous 
    # histograms with no. of classes = Frriedman-Diaconis
    if(type=="hist"){
        xx <- c(data.A[, xlab.A], data.B[, xlab.B])
        span <- 2 * IQR(xx)/n^(1/3)
        cat(span, fill=TRUE)
        cat(ceiling((max(xx)-min(xx))/span), fill=TRUE)
        
        brk <- seq(from = min(xx), to = max(xx), by = span)
        cxA <- cut(x = data.A[, xlab.A], breaks = brk, include.lowest = TRUE)
        cxB <- cut(x = data.B[, xlab.B], breaks = brk, include.lowest = TRUE)
        # tables
        tA <- prop.table(table(cxA))
        if(!is.null(w.A)) tA <- prop.table(xtabs(data.A[,w.A]~cxA))
        
        tB <- prop.table(table(cxB))
        if(!is.null(w.B)) tB <- prop.table(xtabs(data.B[,w.B]~cxB))
        # calculates tot. var dist.
        tvd <- 100*(1/2)*sum(abs(tA-tB))
        labtvd <- paste("tvd ", 
                        paste(round(tvd, 2), "%", sep=" "), 
                        sep="= ")
        btmlab <- paste(pr.lab, labtvd, sep=", ")
        
        dfA <- cbind(data.frame(tA), sample="A")
        dfB <- cbind(data.frame(tB), sample="B")
        colnames(dfA) <- colnames(dfB) <- c("Var1", "Freq", "sample")
        df <- rbind(dfA, dfB)
        
        colnames(df) <- c("Var1", "Freq", "sample")
        
        out <- ggplot(data = df, 
               aes(x = Var1, y = Freq, fill = sample)) +
            geom_bar(stat = "identity", position = "dodge") +
            labs(x = btmlab, y="rel freq")
        
        
    }
    
    
    plot(out)

}