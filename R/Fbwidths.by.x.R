Fbwidths.by.x <-
function (tab.x, tab.xy, tab.xz, 
          deal.sparse="discard", 
          nA=NULL, nB=NULL, ...) 
{
###########################################
# Cohen effect size wrt uniform distribution a
# as a measure of sparseness
# 
    
    cohen.ef <- function(tab){
        nc <- length(tab)
        eqpr <- 1/nc
        relfreq <- prop.table(tab)
        sqrt(sum((relfreq - eqpr)^2/eqpr))
    }
#########################################
# relative frequencies    
    
    N <- sum(tab.xy) + sum(tab.xz)
    prop.x <- prop.table(tab.x)
    prop.xy <- prop.table(tab.xy)
    prop.xz <- prop.table(tab.xz)
###    check variable names in the table
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

# creates a list with all possible combinations of 
# 1, 2, ... of X variables
    n.x <- length(lab.x)
    appo.var <- as.list(lab.x)
    for (k in 2:n.x) {
        b <- combn(lab.x, k)
        b <- data.frame(b, stringsAsFactors = FALSE)
        appo.var <- c(appo.var, as.list(b))
    }

# sample sizes (needed to compute sparseness)    
    if(is.null(nA)) nA <- round(sum(tab.xy))
    if(is.null(nB)) nB <- round(sum(tab.xz))
    max.x.cells <- length(tab.x)
# computating average conditional bounds width for each possible
# combination of the Xs
    H <- length(appo.var)
    out.rng <- as.list(as.numeric(H))
    #av.rng <- matrix(NA, H, 8)
    av.rng <- matrix(NA, H, 13)
    
    for (h in 1:H) {
        lab <- appo.var[[h]]
        p.x <- match(lab, lab.x)
       
        xx <- margin.table(prop.x, p.x)
        av.rng[h, 1] <- length(xx) # cells in tables crossing Xs
        av.rng[h, 2] <- sum(xx == 0) # 0s in tables crossing Xs
        p.xy <- match(c(lab, lab.y), lab.xy)
        xy <- margin.table(prop.xy, p.xy)
        av.rng[h, 3] <- length(xy) # cells in X*Y
        av.rng[h, 4] <- sum(xy == 0) # empty cells in X*Y
        p.xz <- match(c(lab, lab.z), lab.xz)
        xz <- margin.table(prop.xz, p.xz)
        av.rng[h, 5] <- length(xz) # cells in X*Z
        av.rng[h, 6] <- sum(xz == 0) # empty cells in X*Z
        
        avrg.n <- min(nA/length(xy), nB/length(xz)) # min average size
        av.rng[h, 7] <- avrg.n 
        av.rng[h, 8] <- cohen.ef(xx) # Cohen's effect size wrt unif
        
        penalty1 <- log(1 + length(xx)/max.x.cells ) # 1st penalty of av.width
        penalty2 <- max(1/(nA - length(xy)), 1/(nB - length(xz))) # 2nd penalty of av.width
        
        if(deal.sparse=='relfreq'){ 
            #estimates uncertainty without caring of sparseness
            
            fb <- Frechet.bounds.cat(xx, xy, xz, print.f = "tables", 
                                       ...)
            appo <- data.frame(fb$low.cx)
            out.rng[[h]] <- data.frame(appo[, 1:2], lower = c(fb$low.cx), 
                                       upper = c(fb$up.cx), 
                                       width = c(fb$up.cx - fb$low.cx),
                                       CIA = c(fb$CIA))
            
            av.rng[h, 9] <- fb$uncertainty[2]
            av.rng[h, 10] <- penalty1
            av.rng[h, 11] <- fb$uncertainty[2] + penalty1
            av.rng[h, 12] <- penalty2
            av.rng[h, 13] <- fb$uncertainty[2] + penalty2    
                
        }
        if(deal.sparse=='discard'){
            # if table is not sparse procede with estimation 
            if(avrg.n > 1){
                fb <- Frechet.bounds.cat(xx, xy, xz, print.f = "tables", 
                                           ...)
                appo <- data.frame(fb$low.cx)
                out.rng[[h]] <- data.frame(appo[, 1:2], lower = c(fb$low.cx), 
                                           upper = c(fb$up.cx), 
                                           width = c(fb$up.cx - fb$low.cx),
                                           CIA = c(fb$CIA))
                
                av.rng[h, 9] <- fb$uncertainty[2]
                av.rng[h, 10] <- penalty1
                av.rng[h, 11] <- fb$uncertainty[2] + penalty1
                av.rng[h, 12] <- penalty2
                av.rng[h, 13] <- fb$uncertainty[2] + penalty2
            }
            else{
            # if table is sparse no sestimation is performed    
                out.rng[[h]] <- NULL
                av.rng[h, 9] <- NA
                av.rng[h, 10] <- penalty1
                av.rng[h, 11] <- NA
                av.rng[h, 12] <- penalty2
                av.rng[h, 13] <- NA
            }
        }
        # if(deal.sparse=='pseudoB'){
        # # estimation of rel. frequencies is performed using 
        # # pseudoBayes estimator
        #     if(length(p.x)>1) pb.xx <- pBayes(x=xx, method="m.ind")
        #     else pb.xx <- xx
        #     pb.xy <- pBayes(x=xy, method="m.ind")
        #     pb.xz <- pBayes(x=xz, method="m.ind")
        #     fb <- Frechet.bounds.cat(pb.xx, pb.xy, pb.xz, print.f = "tables", 
        #                              ...)
        #     appo <- data.frame(fb$low.cx)
        #     out.rng[[h]] <- data.frame(appo[, 1:2], lower = c(fb$low.cx), 
        #                                upper = c(fb$up.cx), 
        #                                width = c(fb$up.cx - fb$low.cx),
        #                                CIA = c(fb$CIA))
        #     
        #     av.rng[h, 9] <- fb$uncertainty[2]
        #     av.rng[h, 10] <- penalty1
        #     av.rng[h, 11] <- fb$uncertainty[2] + penalty1
        #     av.rng[h, 12] <- penalty2
        #     av.rng[h, 13] <- fb$uncertainty[2] + penalty2
        #     
        # }
            
    }
    lab.list <- paste("|", lapply(appo.var, paste, collapse = "+"), 
                      sep = "")
    n.vars <- lapply(appo.var, length)
    
    # summary information 
    av.rng <- data.frame(x.vars = unlist(n.vars), 
                         x.cells = av.rng[, 1], x.freq0 = av.rng[, 2], 
                         xy.cells = av.rng[, 3], xy.freq0 = av.rng[, 4],
                         xz.cells = av.rng[, 5], xz.freq0 = av.rng[, 6],
                         av.n = av.rng[, 7], cohen.ef=av.rng[, 8],
                         av.width = av.rng[, 9],
                         penalty1 = av.rng[, 10], av.width.pen1 = av.rng[, 11], 
                         penalty2 = av.rng[, 12], av.width.pen2 = av.rng[, 13])
    
    row.names(av.rng) <- paste("|", lapply(appo.var, paste, collapse = "*"), 
                               sep = "")
    # summary information in the unconditional case (no Xs)
    av.rng.0 <- c(x.vars=0, x.cells=NA, x.freq0=NA, 
                  xy.cells = NA, xy.freq0 = NA,
                  xz.cells = NA, xz.freq0 = NA,
                  av.n = 0, cohen.ef=NA, av.width = fb$uncertainty[1], 
                  penalty1 = 0, av.width.pen1 = fb$uncertainty[1], 
                  penalty2 = 0, av.width.pen2 = fb$uncertainty[1])
    
    # whole summary information
    av.rng <- rbind(unconditional=av.rng.0, av.rng)
    
    aa <- n.x - av.rng$x.vars
    ord.lab <- order(aa, av.rng$av.width, decreasing = TRUE)
    av.rng <- av.rng[ord.lab, ]

    # if(compress.sum){
    #     sp.av <- split(av.rng, av.rng$x.vars)
    #     G <- length(sp.av)
    #     sp.new <- as.list(G)
    #     sp.new[[1]] <- sp.av[[1]]
    #     sp.new[[2]] <- sp.av[[2]]
    #     for(g in 3:G){
    #         min.p <- min(sp.av[[(g-1)]][,"av.width"])
    #         tst <- sp.av[[g]][,"av.width"] <= min.p
    #         sp.new[[g]] <- sp.av[[g]][tst,]
    #     }
    #     av.rng <- do.call("rbind", sp.new)
    # }
    out.rng[[(H + 1)]] <- av.rng
    names(out.rng) <- c(lab.list, "sum.unc")
    out.rng
}