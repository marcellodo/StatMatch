comp.cont <- function (data.A, data.B, xlab.A, xlab.B = NULL, w.A = NULL, 
          w.B = NULL, ref = FALSE) 
{
    ###################################################
    hist.bks <- function(x, w=NULL, n=NULL){
        
        xx <- x[!is.na(x)]
        
        if(is.null(w)) {
            # iqr.x <- IQR(xx)
            qq <- quantile(x = xx, probs = c(0.25, 0.5, 0.75))
            deff.uw <- 1
        }
        else{
            ww <- w[!is.na(x)]
            qq <- Hmisc::wtd.quantile(x = xx, weights = ww, 
                                      probs = c(0.25, 0.5, 0.75))
            deff.uw <- length(ww) * sum(ww^2)/(sum(ww)^2)
        }
        iqr.x <- c(qq[3] - qq[1])
        
        if(is.null(n)) nn <- length(xx) / deff.uw
        else nn <- n
        wid.FD <- 2 * iqr.x* 1/(nn^(1/3))
        # robust boxplot fences
        low <- max(min(xx), qq[1] - 1.5*2*(qq[2]-qq[1]) )
        up <-  min(max(xx), qq[3] + 1.5*2*(qq[3]-qq[2]) ) 
        
        bins <- (up - low)/wid.FD
        bins <- floor(bins)+1
        if(bins==1){
            bins <- 2
            wid.FD <- (max(xx) - min(xx))/bins
        }
        # adjust for low and up
        span <- bins*wid.FD
        dd <- span - (up - low)
        low <- low - dd/2
        up <- up + dd/2
        # final breaks
        seq(from=low, to=up, by=wid.FD)
    }
#################################################################    
    # preparation
    nA <- nrow(data.A)
    nB <- nrow(data.B)
    if (is.null(xlab.B)) 
        xlab.B <- xlab.A
    xA <- data.A[,xlab.A]
    xB <- data.B[,xlab.B]
    
    if(!is.null(w.A)) wA <- data.A[, w.A]/sum(data.A[, w.A]) * nA
    else wA <- rep(1, nA)
    deff.uwA <- length(wA) * sum(wA^2)/(sum(wA)^2)
    
    if(!is.null(w.B)) wB <- data.B[, w.B]/sum(data.B[, w.B]) * nB
    else wB <- rep(1, nB)
    deff.uwB <- length(wB) * sum(wB^2)/(sum(wB)^2)
    
    #############################################
    # empirical cumulative distance 
    if(ref) usx <- unique(sort(xB))
    else usx <- unique(sort(c(xA ,xB)))
    
    k <- length(usx)
    ecdf.xA <- ecdf.xB <- numeric(k)
    for(i in 1:k){
        ecdf.xA[i] <- sum(wA[xA <= usx[i]])
        ecdf.xB[i] <- sum(wB[xB <= usx[i]])
    }
    d.KS <- max(abs(ecdf.xA - ecdf.xB))
    d.Kui <- max(ecdf.xA - ecdf.xB) + max(ecdf.xB - ecdf.xA)
    d.Wass <- mean(abs(ecdf.xA - ecdf.xB))
    # reld.Wass <- (1/k) * sum(abs(ecdfA - ecdfB)/(ecdfB*(1-ecdfB)))
    
    d.ecdf <- c(KSdist=d.KS, Kuiper.dist=d.Kui, avAbsDiff=d.Wass)
    
    ############################################
    # quantiles
    n <- min(nA, nB)
    if(n <= 50) pctp <- c(0.25, 0.50, 0.75)
    if(n > 50 & n <= 150) pctp <- seq(from = 0.2, to = 0.8, by = 0.2)
    if(n > 150 & n <= 250) pctp <- seq(from = 0.1, to = 0.9, by = 0.1)
    if(n > 250) pctp <- seq(from = 0.05, to = 0.95, by = 0.05)
    
    qA <- quantile(x = xA, probs = pctp)
    if (!is.null(w.A)) {
        qA <- Hmisc::wtd.quantile(x = xA, weights = wA, probs = pctp)
    }
    
    qB <- quantile(x = xB, probs = pctp)
    if (!is.null(w.B)) {
        qB <- Hmisc::wtd.quantile(x = xB, weights = wB, probs = pctp)
    }
    
    dQa <- mean(abs(qA-qB))
    dQ2 <- mean((qA-qB)^2)
    
    #################################################
    # basic summaries
    if(is.null(w.A)) smryA <- c(summary(xA), sd=sd(xA, na.rm = TRUE))
    else{
        mA <- weighted.mean(x = xA, w = wA, na.rm=TRUE)
        sdA <- sqrt(sum(wA * ((xA - mA)^2))/(sum(wA) - 1))
        QA <- Hmisc::wtd.quantile(x = xA, weights = wA, probs = c(0, 0.25, 0.50, 0.75, 1))
        smryA <- c(QA[1:3], mA, QA[4:5], sdA)
        names(smryA) <- NULL
        names(smryA) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "sd")
    }
    if(is.null(w.B)) smryB <- c(summary(xB), sd=sd(xB, na.rm = TRUE))
    else{
        mB <- weighted.mean(x = xB, w = wB, na.rm=TRUE)
        sdB <- sqrt(sum(wB * ((xB - mB)^2))/(sum(wB) - 1))
        QB <- Hmisc::wtd.quantile(x = xB, weights = wB, probs = c(0, 0.25, 0.50, 0.75, 1))
        smryB <- c(QB[1:3], mB, QB[4:5], sdB)
        names(smryB) <- NULL
        names(smryB) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "sd")
    }
    smrys <- as.data.frame(dplyr::bind_rows(A=smryA, B=smryB))
    
    ###############################################
    # total variation distance (discretize continuous)
    if(ref)  {
        bks <- hist.bks(x=xB, w=wB, n=min(nA/deff.uwA, nB/deff.uwB))
        if(bks[1] > min(xA)) bks[1] <- min(xA)
        if(bks[length(bks)] < max(xA)) bks[length(bks)] <- max(xA)
    }
    else {
        bks <- hist.bks(x=c(xA, xB), w=c(wA, wB), min(nA/deff.uwA, nB/deff.uwB))
    } 
    
    cxA <- cut(x = xA, breaks = bks, include.lowest = TRUE)
    cxB <- cut(x = xB, breaks = bks, include.lowest = TRUE)
    # tables
    pA <- prop.table(table(cxA))
    if(!is.null(w.A)) tA <- prop.table(xtabs(data.A[,w.A]~cxA))
    
    pB <- prop.table(table(cxB))
    if(!is.null(w.B)) tB <- prop.table(xtabs(data.B[,w.B]~cxB))
    
    tvd <- 1/2 * sum(abs(pA-pB))
    bhatt <- sum(sqrt(pA * pB))
    
    
    ##########################
    # output
    list(summary = smrys,
         diff.Qs = c(avAbsD=dQa, avSqrtSqD=sqrt(dQ2)),
         dist.ecdf = d.ecdf,
         dist.discr = c(tvd=tvd, overlap=1-tvd, Hellinger = sqrt(1 - bhatt)) 
    )
}  