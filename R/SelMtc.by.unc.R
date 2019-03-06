selMtc.by.unc <- 
function (tab.x, tab.xy, tab.xz, corr.d=2, 
          nA=NULL, nB=NULL, align.margins=FALSE) 
{
#######veq.fcn computes sparseness measure
    veq.fcn <- function(tab){
        if(sum(tab)!=1) tab <- prop.table(tab)
        nc <- length(tab)
        sqrt(
            sum((tab - 1/nc)^2 / (1/nc) )
        )
    }
    
################################################    
# vlist.fcn ----------------------------
# all possible combinations of
# x with (all.x != x)
    vlist.fcn <- function(x, all.x){
        pos <- match(x, all.x)
        y <- all.x[-pos]
        k <- length(y)
        out <- as.list(numeric(k))
        for(i in 1:k){
            out[[i]] <- c(x, y[i])
        }
        out
    }
############# end vlist.fcn
##################################################################    
# findbest.fcn ---------------------------
# to identify best X variable to add to the existing combination
# of X variables
    
findbest.fcn <- function(pp.x, pp.xy, pp.xz, x.vars, corr=0, 
                         N.xy=NULL, N.xz=NULL, align=FALSE){
        
        nc.xx.all <- length(pp.x)
        lab.x <- names(dimnames(pp.x))
        
        lab.xy <- names(dimnames(pp.xy))
        lab.y <- setdiff(lab.xy, lab.x)
        p.y <- match(lab.y, lab.xy)
        c.y <- dim(pp.xy)[length(dim(pp.xy))]
        
        lab.xz <- names(dimnames(pp.xz))
        lab.z <- setdiff(lab.xz, lab.x)
        p.z <- match(lab.z, lab.xz)
        c.z <- dim(pp.xz)[length(dim(pp.xz))]
        
        H <- length(x.vars)
        
        av <- pen <- av.pen <- numeric()
        cell.info <- matrix(0, nrow=H, ncol=12)
        
        for (h in 1:H) {
            lab <- x.vars[[h]]
            p.x <- match(lab, lab.x)
            xx <- margin.table(pp.x, p.x)
            nc.xx <- length(xx)
            nc0.xx <- sum(xx==0)
            av.xx <- mean(xx)
            
            p.xx <- prop.table(xx)
            veq.xx <- veq.fcn(p.xx)
            
            p.xy <- match(c(lab, lab.y), lab.xy)
            xy <- margin.table(pp.xy, p.xy)
            nc.xy <- length(xy)
            nc0.xy <- sum(xy==0)
            av.xy <- mean(xy)
            veq.xy <- veq.fcn(xy)
            
            p.xz <- match(c(lab, lab.z), lab.xz)
            xz <- margin.table(pp.xz, p.xz)
            nc.xz <- length(xz)
            nc0.xz <- sum(xz==0)
            av.xz <- mean(xz)
            veq.xz <- veq.fcn(xz)
            
            cell.info[h,] <- c(nc.xx, nc0.xx, av.xx, veq.xx, 
                               nc.xy, nc0.xy, av.xy, veq.xy,
                               nc.xz, nc0.xz, av.xz, veq.xz)
            
            #if( mean(xx) < 5 ) xx <- pBayes(xx, method=p.meth)$pseudoB
            #if( mean(xy) < 5 ) xy <- pBayes(xy, method=p.meth)$pseudoB
            #if( mean(xz) < 5 ) xz <- pBayes(xz, method=p.meth)$pseudoB
            
            # measure sparseness 
            if(is.null(N.xy)) N.xy <- sum(xy)
            if(is.null(N.xz)) N.xz <- sum(xz)
            
            av.n <- min(N.xy/nc.xy, N.xz/nc.xz) 
            if(av.n>1){
                fb <- Frechet.bounds.cat(xx, xy, xz, print.f = "tables", 
                                         align.margins=align)
                
                if(corr==0) {
                    cc <- 0
                }
                if(corr==1){
                    cc <- log(1 + nc.xx/nc.xx.all) 
                    if(is.na(cc)) cc <- 0
                }
                if(corr==2){
                    cc <- max(1/(nA-nc.xy), 1/(nB-nc.xz)) 
                    if(is.na(cc)) cc <- 0
                }
                av[h] <- fb$uncertainty[2]
                pen[h] <- cc
                av.pen[h] <- fb$uncertainty[2] + cc    
            }
            else {
                av[h] <- NA
                pen[h] <- NA
                av.pen[h] <- NA
            }
        }
        
        if(all(is.na(av.pen))){
            out <- list(ord=NA,
                        bst=NA, 
                        av=NA,
                        pen=NA,
                        av.pen=NA,
                        c.info=NA)
        }
        else{
            pos <- which.min(av.pen)
            ord.av <- order(av.pen, decreasing=FALSE, na.last = TRUE)
            out <- list(ord=ord.av,
                        bst=x.vars[[pos]], 
                        av=av[ord.av],
                        pen=pen[ord.av],
                        av.pen=av.pen[ord.av],
                        c.info=cell.info[ord.av,])
        }   
        out
    }
###################### end findbest.fcn ######################
##############################################################    
# computations start here    

    if(!(corr.d %in% 0:2)) stop("argument corr.d should be 0, 1, or 2" ) 
#sample sizes
    if(is.null(nA)) nA <- sum(tab.xy)
    if(is.null(nB)) nB <- sum(tab.xz)
    if(nA==1 | nB==1) stop("sample sizes of samples A and B are needed")
    
# names of variables 
    lab.x <- names(dimnames(tab.x))
    if (all(nchar(lab.x) == 0)) 
        lab.x <- paste("x", 1:length(lab.x), sep = "")
    names(attr(tab.x, "dimnames")) <- lab.x

    lab.xy <- names(dimnames(tab.xy))
    if (all(nchar(lab.xy) == 0)) 
        lab.xy <- c(lab.x, "y")
    names(attr(tab.xy, "dimnames")) <- lab.xy
    
    lab.y <- setdiff(lab.xy, lab.x)
    p.y <- match(lab.y, lab.xy)
    
    lab.xz <- names(dimnames(tab.xz))
    if (all(nchar(lab.xz) == 0)) 
        lab.xz <- c(lab.x, "z")
    names(attr(tab.xz, "dimnames")) <- lab.xz
    
    lab.z <- setdiff(lab.xz, lab.x)
    p.z <- match(lab.z, lab.xz)
    
    # fbc.0 = Frechet.bounds.cat(tab.x=NULL, tab.xy = margin.table(tab.xy, p.y),
    #                            tab.xz = margin.table(tab.xz, p.z),
    #                            print.f = "tables")
# summary info    
    n.x <- length(lab.x)
    
# step 0) order available X vars by small av.width
    fnd.0 <- findbest.fcn(pp.x=tab.x, pp.xy=tab.xy, pp.xz=tab.xz, 
                          x.vars=as.list(lab.x), corr=corr.d, 
                          N.xy = nA, N.xz = nB,
                          align = align.margins)
    lab.x <- lab.x[fnd.0$ord]
    avw0 <- fnd.0$av.pen
    names(avw0) <- lab.x
    
    list.vars.x <- as.list(numeric(n.x))

    best.x <- list.vars.x[[1]] <- lab.x[1]
    av.w <- pen <- av.w.pen <- rep(NA, n.x)
    c.info <- matrix(NA, nrow=n.x, ncol=12)

    # best model with a single X
    av.w[1] <- fnd.0$av[1]
    pen[1] <- fnd.0$pen[1]
    av.w.pen[1] <- fnd.0$av.pen[1]
    c.info[1, ] <- fnd.0$c.info[1,]
    # i <- 2
    # chk <- 1
    #while(!is.na(chk) & i < (n.x-1)){
# start searching best model with increasing combinations of Xs
# until the stopping rule is satisfied
    # ref.d <- fnd.0$av.pen[1]
    for(i in 2:(n.x-1)){
        test.x <- vlist.fcn(x=best.x, all.x=lab.x)
        fnd <- findbest.fcn(pp.x=tab.x, pp.xy=tab.xy, pp.xz=tab.xz, 
                            x.vars=test.x, N.xy = nA, N.xz = nB,
                            corr = corr.d, 
                            align = align.margins)
        chk <- fnd$bst
        
        # if(is.na(chk) | fnd$av.pen[1]>ref.d) break
        if(is.na(chk)) break
        else {
            # ref.d <- fnd$av.pen[1]
            best.x <- fnd$bst
            list.vars.x[[i]] <- fnd$bst
            av.w[i] <- fnd$av[1]
            pen[i] <- fnd$pen[1]
            av.w.pen[i] <- fnd$av.pen[1]
            #cat(fnd$c.info, fill=T)
            c.info[i, ] <- fnd$c.info[1,]
        }    
    }

### full model (all x included)
    if(min(nA/length(tab.xy), nB/length(tab.xz)) > 1 ){
        fnd <-  Frechet.bounds.cat(tab.x=tab.x, tab.xy=tab.xy, tab.xz=tab.xz,
                                   print.f = "tables",
                                   align.margins = align.margins)
        
        av.w[n.x] <- fnd$uncertainty[2]
        if(corr.d==0) {
            pen[n.x] <- 0
            av.w.pen[n.x] <- fnd$uncertainty[2]
        }
        if(corr.d==1) {
            pen[n.x] <- log(2)
            av.w.pen[n.x] <- fnd$uncertainty[2] + log(2)
        }    
        if(corr.d==2) {
            pen[n.x] <- max(1/(nA-length(tab.xy)),1/(nB-length(tab.xz)))
            av.w.pen[n.x] <- fnd$uncertainty[2] + max(1/(nA-length(tab.xy)),1/(nB-length(tab.xz)))
        }
        c.info[n.x, ] <- c(length(tab.x), sum(tab.x==0),  mean(tab.x), 
                           veq.fcn(tab.x),
                           length(tab.xy), sum(tab.xy==0), mean(tab.xy), 
                           veq.fcn(tab.xy), 
                           length(tab.xz), sum(tab.xz==0), mean(tab.xz),
                           veq.fcn(tab.xz)
        )
    }
# final output
    colnames(c.info) <- c("nc.x", "nc0.x", "av.crf.x", 'veq.x',
                          "nc.xy", "nc0.xy", "av.crf.xy",'veq.xy',
                          "nc.xz", "nc0.xz", "av.crf.xz", 'veq.xz')
    c.info <- cbind(c.info, 
                    min.av=apply(c.info[, c("av.crf.xy","av.crf.xz")], 1, min))
    list.vars.x[[n.x]] <- lab.x

    vv <- unlist(lapply(list.vars.x, paste, collapse="*"))
    nv <- unlist(lapply(list.vars.x, length))
    
    aa <- data.frame( avw=av.w, penalty=pen, avw.pen=av.w.pen)
    av.df <- cbind(x.vars=vv, nxv=nv, c.info, aa)
    
    ty <- margin.table(tab.xy, p.y)
    tz <- margin.table(tab.xz, p.z)
    
# case of unconditional uncertainty (no Xs)    
    fbc0 <- Frechet.bounds.cat(tab.x=NULL, tab.xy = ty, tab.xz=tz)
   
    l0 = data.frame(x.vars=NA, nxv=0, 
                    nc.x=NA, nc0.x=NA, av.crf.x=NA, veq.x=NA,
                    nc.xy= nrow(fbc0$low.u), nc0.xy=sum(ty==0), 
                    av.crf.xy=mean(ty), veq.xy=veq.fcn(ty),
                    nc.xz= ncol(fbc0$low.u), nc0.xz=sum(tz==0), 
                    av.crf.xz=mean(tz), veq.xz=veq.fcn(tz),
                    min.av=min(mean(ty), mean(tz)), 
                    avw=fbc0$uncertainty[1],
                    penalty=0, avw.pen=fbc0$uncertainty[1]
                    )
    #cat(dim(l0), fill=T)
# final output
    av.df <- rbind(l0, av.df)
    tst <- !is.na(av.df$avw) 
    av.df.ok <- subset(av.df, tst) # fulfill sparseness rule
    # av.df.no <- subset(av.df, tst) # NOT fulfill sparseness rule
    pos <- which.min(av.df.ok$avw.pen)
    list(ini.ord=avw0, list.xs=list.vars.x[tst[-1]], 
         av.df=av.df.ok[1:pos, ])
}
