`pw.assoc` <- 
function(formula, data, weights=NULL, out.df=FALSE)
{
#####################################################################
## code for Cramer's V (standard and bias-corrected expression)
    V <- function(tab, n=NULL){
        if(is.null(n)) n <- sum(tab)
        if(n==1) warning("n = 1 will provide an erroneous estimate of \n 
                     bias-corrected Cramer's V")
        nr <- nrow(tab)
        nc <- ncol(tab)
        chi <- chisq.test(tab)
        mm <- min(nr - 1, nc - 1)
        stdV <- sqrt(chi$statistic/(n*mm)) # Cramer's V
        names(stdV) <- NULL
        bb <- chi$statistic/n - 1/(n-1)*(nr - 1)*(nc - 1)
        bb <- max(0, bb)
        avr <- nr - 1/(n-1)*(nr-1)^2
        avc <- nc - 1/(n-1)*(nc-1)^2
        bcV <- sqrt(bb/min(avr-1, avc-1)) # bias-corrected Cramer's V
        out <- c(V=stdV, bcV=bcV)
    }
#####################################################################
### function for computing: 
### (a) mutual information based measures
### (b) measures of proportional reduction of the variance Row|Column

    prv.rc <- function(tab){
        tab <- tab/sum(tab)
        rS <- rowSums(tab)
        cS <- colSums(tab)
        ## code for Goodman & Kruskal lambda(Row|Column)    
        V.r <- 1 - max(rS)
        EV.rgc <- 1-sum(apply(tab, 2, max))
        lambda <- (V.r - EV.rgc)/V.r
        ## code for Goodman & Kruskal tau(Row|Column) 
        V.r <- 1 - sum(rS^2)
        a <- colSums(tab^2)
        EV.rgc <- 1 - sum(a/cS)
        tau <- (V.r - EV.rgc)/V.r
        ## code for Theil's Uncertainty(Row|Column)    
        H.fcn <- function(x){
            x <- c(x)
            xx <- x[x>0]
            pp <- xx/sum(xx)
            -1 * sum(pp*log(pp)) 
        }
        # V.r <- (-1)*sum(rS*log(rS))
        # cS.mat <- matrix(cS, nrow=nrow(tab), ncol=ncol(tab), byrow=TRUE)
        # EV.rgc <- sum(tab *log(tab/cS.mat))
        # u <- (V.r + EV.rgc)/V.r
        
        H.r <- H.fcn(rS)
        H.c <- H.fcn(cS)
        mI <- (H.r + H.c - H.fcn(tab)) 
        u <- mI/H.r
        nI <- mI/min(H.r, H.c)
        
## output
        c(lambda.rc=lambda, tau.rc=tau, U.rc=u, I=mI, nI=nI)
    }
######################################################################
## function for computing AIC and BIC of Row|col    
    IC.based <- function(x, n=NULL){
        
        #####    
        if(is.null(n)) n <- sum(x)
        rS <- rowSums(x)
        cS <- colSums(x)
        # marginal distr. response variable (in row)
        IC.0 <- sum(rS * log(rS/n), na.rm = TRUE)   
        df.0 <- length(rS) - 1
        if(any(rS==0)) df.0 <- sum(rS>0) - 1
        
        # response (Row) conditional on Col
        IC.1 <- sum(x * log(prop.table(x=x, margin=2)), 
                    na.rm = TRUE)
        df.1 <- df.0 * length(cS)
        if(any(cS==0)) df.1 <- df.0 * sum(cS>0)
        # df.1a <- sum(tab>0) - df.0    
        
        # computes AIC & BIC
        AIC.0 <- (-2) * IC.0 + 2 * df.0
        BIC.0 <- (-2) * IC.0 + log(n) * df.0
        
        AIC.1 <- (-2) * IC.1 + 2 * df.1
        BIC.1 <- (-2) * IC.1 + log(n) * df.1
        # AIC.1a <- (-2) * IC.1 + 2 * df.1a
        # BIC.1a <- (-2) * IC.1 + log(n) * df.1a
        
        
        v0 <- c(df=df.0, AIC=AIC.0, BIC=BIC.0)
        v1 <- c(df=df.1, AIC=AIC.1, BIC=BIC.1)
        out <- rbind(v0, v1, v1-v0)
        row.names(out) <- c("row", "row|col", "diff")
        out
        
    }
    
###################################################
###################################################
    n <- nrow(data)
    if(is.null(weights)) ww <- rep(1, nrow(data))
    else{
        ww <- data[,weights]
        ww <- ww/sum(ww)*n
        data <- data[,setdiff(colnames(data), weights)]
    }
    # n <- sum(ww)
    
    df <- model.frame(formula=formula, data=data)
    lab <- colnames(df)
    p.x <- length(lab) - 1
    vV <- vbcV <- vlambda <- vtau <- vU  <- vI <- vnI <- vAIC <- vBIC <- vdf <- numeric(p.x)
    for(i in 1:p.x){
        pos <- i+1
        form <- paste(lab[1], lab[pos], sep="+")
        form <- paste("ww", form, sep="~")
        tab <- xtabs(as.formula(form), data=df)
        Vs <- V(tab)
        vV[i] <- Vs[1] #standard Cramer's V
        vbcV[i] <- Vs[2] # bias corrected Cramer's V
        
        appo <- prv.rc(tab)
        vI[i] <- appo["I"] # mutual info
        vnI[i] <- appo["nI"] # normalized mutual info
        vlambda[i] <- appo["lambda.rc"]
        vtau[i] <- appo["tau.rc"]
        vU[i] <- appo["U.rc"]
        
        ics <- IC.based(x=tab, n=n)
        vAIC[i] <- ics["row|col", "AIC"]
        vBIC[i] <- ics["row|col", "BIC"]
        vdf[i] <- ics["row|col", "df"]
    }
    lab.am.assoc <- paste(lab[1], lab[-1], sep="_" )
    lab.am.cond <- paste(lab[1], lab[-1], sep="|" )
    names(vV) <- names(vbcV) <- names(vI) <- names(vnI) <-  lab.am.assoc
    names(vlambda) <- names(vtau) <- names(vU) <- names(vAIC) <- names(vBIC) <- names(vdf) <- lab.am.cond
    
    fine <- list(V=vV, bcV=vbcV, mi=vI, norm.mi=vnI, 
                 lambda=vlambda, tau=vtau, U=vU, 
                 AIC=vAIC, BIC=vBIC, npar=vdf)
    if(out.df){
        fine <- do.call("cbind", fine)
        colnames(fine) <- c("V", "bcV", "mi", "norm.mi",
                            "lambda", "tau", "U", "AIC", "BIC", "npar")
        fine <- as.data.frame(fine)
    }
    fine
}
