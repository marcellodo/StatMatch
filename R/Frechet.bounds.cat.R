'Frechet.bounds.cat' <-
function(tab.x, tab.xy, tab.xz, print.f="tables", tol= 0.0001)
{
    lab.x <- names(dimnames(tab.x))
    if(all(nchar(lab.x)==0)) lab.x <- paste("x",1:length(lab.x), sep="")
    names(attr(tab.x, "dimnames")) <- lab.x

    lab.xy <- names(dimnames(tab.xy))
    if(all(nchar(lab.xy)==0)) lab.xy <- c(lab.x, "y")
    names(attr(tab.xy, "dimnames")) <- lab.xy

    lab.y <- setdiff(lab.xy, lab.x)
    pos.y <- match(lab.y, lab.xy)

    lab.xz <- names(dimnames(tab.xz))
    if(all(nchar(lab.xz)==0)) lab.xz <- c(lab.x, "z")
    names(attr(tab.xz, "dimnames")) <- lab.xz

    lab.z <- setdiff(lab.xz, lab.x)
    pos.z <- match(lab.z, lab.xz)

    p.x <- prop.table(tab.x)
    p.xy <- prop.table(tab.xy)
    p.y <- margin.table(p.xy, pos.y)

    p.xz <- prop.table(tab.xz)
    p.z <- margin.table(p.xz, pos.z)

    ndim <- 1:length(dim(p.xy))
    y.gx <- prop.table(p.xy, ndim[-pos.y])

    ndim <- 1:length(dim(p.xz))
    z.gx <- prop.table(p.xz, ndim[-pos.z])

# check marginal distribution of the X variables

    d1.x <- 1:length(dim(p.xy))
    m1.x <- margin.table(p.xy, d1.x[-pos.y])
    if(any(abs(m1.x-p.x)>tol) )
        warning("The marginal distr. of the X variables \n in tab.xy is not equal to tab.x")

    d2.x <- 1:length(dim(p.xz))
    m2.x <- margin.table(p.xz, d2.x[-pos.z])
    if(any(abs(m2.x-p.x)>tol) )
        warning("The marginal distr. of the X variables \n in tab.xz is not equal to tab.x")

    if(any(abs(m1.x-m2.x)>tol) )
        warning("The marginal distr. of the X variables \n in tab.xy and in tab.xz are not equal")

########################################################
# computes Fréchet bounds _without_ using X variables

    ll <- outer(p.y, p.z, FUN="+") - 1
    m0 <- matrix(0, nrow(ll), ncol(ll))
    low <- pmax(m0, ll)
    upper <- outer(p.y, p.z, FUN="pmin")

    dimnames(low) <- dimnames(upper) <- list(names(p.y), names(p.z))
    class(low) <- class(upper)  <- "table"
    res.0 <- list(low.u=low, up.u=upper)

#############################################
# computes Fréchet bounds using X variables

    dm.x <- data.frame(p.x)
    sdm.x <- split(dm.x, dm.x[,lab.x])

    ay.gx <- data.frame(y.gx)
    say.gx <- split(ay.gx, ay.gx[,lab.x])

    bz.gx <- data.frame(z.gx)
    sbz.gx <- split(bz.gx, bz.gx[,lab.x])

    H <- length(sdm.x)
    out.CIA <- as.list(numeric(H))
    out.low <- as.list(numeric(H))
    out.up <- as.list(numeric(H))
    unc <- as.list(numeric(H)) 

    for(h in 1:H){
        yy <- say.gx[[h]][,"Freq"]
        yy[is.nan(yy)] <- 0
        zz <- sbz.gx[[h]][,"Freq"]
        zz[is.nan(zz)] <- 0
        xx <- sdm.x[[h]][,"Freq"]
        xx[is.nan(xx)] <- 0
        out.CIA[[h]] <- outer(yy, zz, FUN="*") * xx

        thetas <- outer(yy, zz, FUN="+")-1
        m0 <- matrix(0, nrow=nrow(thetas), ncol=ncol(thetas))
        ll <- pmax(thetas,m0)
        uu <- outer(yy, zz, FUN="pmin")
        out.low[[h]] <- ll*xx
        out.up[[h]] <- uu*xx
        unc[[h]] <- sum((uu-ll) * outer(yy, zz, FUN="*") * xx)
    }

    aa.CIA <- array(unlist(out.CIA), dim=c(dim(out.CIA[[1]]),H) )
    fine.CIA <- apply(aa.CIA, c(1,2), sum)

    aa.low <- array(unlist(out.low), dim=c(dim(out.low[[1]]),H) )
    fine.low <- apply(aa.low, c(1,2), sum)

    aa.up <- array(unlist(out.up), dim=c(dim(out.up[[1]]),H) )
    fine.up <- apply(aa.up, c(1,2), sum)
    l.y <- dimnames(y.gx)
    p.y <- match(lab.x, names(l.y))
    l.y <- l.y[-p.y]

    l.z <- dimnames(z.gx)
    p.z <- match(lab.x, names(l.z))
    l.z <- l.z[-p.z]
    class(fine.CIA) <- class(fine.low) <- class(fine.up) <- "table"
    dimnames(fine.CIA) <- dimnames(fine.low) <-  dimnames(fine.up) <- c(l.y, l.z)
    res.1 <- list(CIA=fine.CIA, low.cx=fine.low, up.cx=fine.up, unc=sum(unlist(unc)))

    if(print.f=="tables"){
        out <- c(res.0, res.1)
    }
    else if(print.f=="data.frame"){
        dataf <- data.frame(res.0$low.u)
        labdf <- c(lab.y, lab.z, "low.u")
        colnames(dataf) <- labdf
        dataf$low.cx <- c(res.1$low.cx)
        dataf$CIA <- c(res.1$CIA)
        dataf$up.cx <- c(res.1$up.cx)
        dataf$up.u <- c(res.0$up.u)
        out <- list(bounds=dataf, unc=res.1$unc)
    }
    out
}
