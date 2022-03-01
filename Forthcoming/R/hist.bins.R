hist.bins <- function(x, w=NULL, deff.c=1){
    
    xx <- x[!is.na(x)]
    n <- length(xx)
    
    if(is.null(w)){
        
        sd.x <- sd(xx)
        iqr.x <- IQR(xx)
        deff.uw <- 1
    }
    else{
        ww <- w[!is.na(x)]
        
        sd.x <- sqrt(Hmisc::wtd.var(x = xx, weights = ww))
        wQ <- Hmisc::wtd.quantile(x = xx, weights = ww, probs = c(0.25, 0.75))
        iqr.x <- c(wQ[2] - wQ[1])
        deff.uw <- n * sum(ww^2)/(sum(ww)^2)
        
    }
    nn <- n / (deff.c * deff.uw)
    
    A <- min(sd.x, iqr.x/1.349)
    
    wid=c(Scott= 3.49 * sd.x *1/(nn^(1/3)),
          FD= 2 * iqr.x* 1/(nn^(1/3)),
          adj = 2.59 * A * 1/(nn^(1/3)))
    bin <- (max(xx)-min(xx))/wid
    szs <- c(n=n)
    if(!is.null(w)) szs <- c(szs, deff.c=deff.c, deff.uw=deff.uw, eff.n =nn)
    
    list(size=szs, stats=c(min=min(xx), max=max(xx), IQR=iqr.x, sd=sd.x), width=wid,
         bins=floor(bin)+1)
    
}
