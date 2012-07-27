'Fbwidths.by.x' <-
function(tab.x, tab.xy, tab.xz)
{
    prop.x <- prop.table(tab.x)
    prop.xy <- prop.table(tab.xy)
    prop.xz <- prop.table(tab.xz)

    lab.x <- names(dimnames(tab.x))
    if(all(nchar(lab.x)==0)) lab.x <- paste("x",1:length(lab.x), sep="")
    names(attr(tab.x, "dimnames")) <- lab.x

    lab.xy <- names(dimnames(tab.xy))
    if(all(nchar(lab.xy)==0)) lab.xy <- c(lab.x, "y")
    names(attr(tab.xy, "dimnames")) <- lab.xy
    lab.y <- setdiff(lab.xy, lab.x)
    p.y <- match(lab.y, lab.xy)

    lab.xz <- names(dimnames(tab.xz))
    if(all(nchar(lab.xz)==0)) lab.xz <- c(lab.x, "z")
    names(attr(tab.xz, "dimnames")) <- lab.xz
    lab.z <- setdiff(lab.xz, lab.x)
    p.z <- match(lab.z, lab.xz)

##
    n.x <- length(lab.x)
    appo.var <- as.list(lab.x)
    for(k in 2:n.x){
        b <- combn(lab.x, k)
        b <- data.frame(b, stringsAsFactors=FALSE)
        appo.var <- c(appo.var, as.list(b))
    }

    H <- length(appo.var)
    out.rng <- as.list(as.numeric(H))
    av.rng <- matrix(NA, H, 4)
    for(h in 1:H){
        lab <- appo.var[[h]]
        p.x <- match(lab, lab.x)
        xx <- margin.table(prop.x, p.x)
        av.rng[h,1] <- length(xx)
        av.rng[h,2] <- sum(xx==0)
        
        p.xy <- match(c(lab,lab.y), lab.xy)
        xy <- margin.table(prop.xy, p.xy)
        
        p.xz <- match(c(lab, lab.z), lab.xz)
        xz <- margin.table(prop.xz, p.xz)

        fb <- Frechet.bounds.cat(xx, xy, xz, print.f="tables")
        appo <- data.frame(fb$low.cx)
        out.rng[[h]] <- data.frame(appo[,1:2], lower=c(fb$low.cx), upper=c(fb$up.cx), width=c(fb$up.cx-fb$low.cx))
        av.rng[h,3] <- mean( c(fb$up.cx-fb$low.cx))
        av.rng[h,4] <- fb$unc
   
    }
    lab.list <- paste("|", lapply(appo.var, paste, collapse="+"), sep="")
    n.vars <- lapply(appo.var, length)
    av.rng <- data.frame(x.vars=unlist(n.vars), x.cells=av.rng[,1], x.freq0=av.rng[,2], 
                         av.width=av.rng[,3], ov.unc=av.rng[,4])
    row.names(av.rng) <- paste("|", lapply(appo.var, paste, collapse="+"), sep="")

    aa <- n.x - av.rng$x.vars
    ord.lab <- order(aa, av.rng$ov, decreasing=TRUE)
    
    out.rng[[(H+1)]] <- av.rng[ord.lab,]
    names(out.rng) <- c(lab.list, "av.widths")
    out.rng
}
