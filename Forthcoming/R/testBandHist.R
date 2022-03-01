h.fn <- function(x, w=NULL){
    
    n <- length(x)
    rng <- max(x) - min(x)
    
    if(is.null(w)){
        iqr.x <- IQR(x)
        n.eff <- n
    }
    else{
        ww <- w/sum(w)*n
        qq <- Hmisc::wtd.quantile(x=x, weights = ww, probs = c(0.25, 0.75))
        iqr.x <- c(qq[2] - qq[1])
        w.eff <- n * sum(ww^2)/(sum(ww)^2)
        n.eff <- n/w.eff
    } 
    h.FD <- 2 * iqr.x / (n.eff^(1/3))
    cat(h.FD, fill=T)
    k <- ceiling(rng/h.FD)
    cat(k, fill=TRUE)
   
    min.k <- k-30
    if(min.k<0) min.k <- 2
    max.k <- k+30
    if(max.k>n.eff) max.k <- n.eff-2
    
    vet.k <- min.k:max.k
    cost <- numeric()
    for(j in 1:length(vet.k)){
        hh <- rng/vet.k[j]
        bks <- seq(from=min(x), to=max(x), by=hh)
        cx <- cut(x, breaks = bks, include.lowest = TRUE)
        if(is.null(w)) fq <- xtabs(~cx)
        else fq <- xtabs(ww~cx)
        m.fq <- mean(fq)
        v.fq <- mean((fq - m.fq)^2)
        cost[j] <- (2*m.fq - v.fq)/(hh^2) 
    }
    vet.k[which.min(cost)]
    
}

data(samp.A, package = "StatMatch")
h.fn(samp.A$n.income)
h.fn(samp.A$n.income, w = samp.A$ww)


h.rot <- function(x, w=NULL){
    n <- length(x)
    if(is.null(w)){
        iqr.x <- IQR(x)
        sd.x <- sd(x)
        N <- n.eff <- n
    }
    else{
        ww <- w/sum(w)*n
        qq <- Hmisc::wtd.quantile(x=x, weights = ww, probs = c(0.25, 0.75))
        iqr.x <- c(qq[2] - qq[1])
        sd.x <- sqrt(Hmisc::wtd.var(x=x, weights = ww))
        w.eff <- n * sum(ww^2)/(sum(ww)^2)
        n.eff <- n/w.eff
        N <- sum(w)
    }
    
    nn <- c(N=N, n=n, n.eff=n.eff)
    h <- 1.06 * min(sd.x, iqr.x/1.34) * 1/(nn^(1/5))
    rbind(n=nn, h=h)
    
    
}
h.rot(x = samp.A$n.income, w=samp.A$ww)
