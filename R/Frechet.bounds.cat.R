Frechet.bounds.cat <- function (tab.x, tab.xy, tab.xz, 
                                print.f = "tables", align.margins = FALSE,
                                tol = 0.001, warn=TRUE) 
{
# fucntion to compute uncertainty bounds in the unconditional case
# i.e. no X variables    
     
    fb.yz <- function(y, z, prn = "tables") {
        lab.y <- names(y)
        if (is.null(lab.y)) 
            lab.y <- paste("y", 1:length(y), sep = "")
        lab.z <- names(z)
        if (is.null(lab.z)) 
            lab.z <- paste("z", 1:length(z), sep = "")
        p.y <- prop.table(y)
        p.z <- prop.table(z)
        ll <- outer(p.y, p.z, FUN = "+") - 1
        m0 <- matrix(0, nrow(ll), ncol(ll))
        low <- pmax(m0, ll)
        upper <- outer(p.y, p.z, FUN = "pmin")
        ind <- outer(p.y, p.z, FUN = "*")
        dimnames(low) <- dimnames(upper) <- dimnames(ind) <- list(lab.y, 
                                                                  lab.z)
        class(low) <- class(upper) <- class(ind) <- "table"
        H.y <- sum(-1 * p.y * log(p.y), na.rm = TRUE)
        H.z <- sum(-1 * p.z * log(p.z), na.rm = TRUE)
        res.0 <- list(low.u = low, up.u = upper, IA = ind, H = c(H.y, 
                                                                 H.z), uncertainty = mean(upper - low))
        if (prn == "tables") {
            out <- res.0
        }
        else if (prn == "data.frame") {
            df <- data.frame(low)
            colnames(df) <- c("Y", "Z", "low.u")
            df$IA <- c(ind)
            df$up.u <- c(upper)
            out <- list(bounds = df, H = c(H.y, H.z), uncertainty = mean(upper - 
                                                                             low))
        }
        out
    }
########################### end fb.yz #################
############################################################    
#
# computation when Xs are NOT available
    if (is.null(tab.x)) 
        out <- fb.yz(y = tab.xy, z = tab.xz, prn = print.f)
# computation when Xs are available 
    else {
        # check variable names
        lab.x <- names(dimnames(tab.x))
        if (all(nchar(lab.x) == 0)) 
            lab.x <- paste("x", 1:length(lab.x), sep = "")
        names(attr(tab.x, "dimnames")) <- lab.x
        
        lab.xy <- names(dimnames(tab.xy))
        if (all(nchar(lab.xy) == 0)) 
            lab.xy <- c(lab.x, "y")
        names(attr(tab.xy, "dimnames")) <- lab.xy
        lab.y <- setdiff(lab.xy, lab.x)
        pos.y <- match(lab.y, lab.xy) #position of Y
        
        lab.xz <- names(dimnames(tab.xz))
        if (all(nchar(lab.xz) == 0)) 
            lab.xz <- c(lab.x, "z")
        names(attr(tab.xz, "dimnames")) <- lab.xz
        lab.z <- setdiff(lab.xz, lab.x)
        pos.z <- match(lab.z, lab.xz) # position of Z
        
        # tables with relative frequencies
        p.x <- prop.table(tab.x)
        p.xy <- prop.table(tab.xy)
        p.xz <- prop.table(tab.xz)
        p.y <- margin.table(p.xy, pos.y) # Y marginal
        p.z <- margin.table(p.xz, pos.z) # Z marginal
        
        # check coherence marginal of X estimated from tab.xy and 
        # tab.xy wrt to p.x 
        if(!align.margins){
            d1.x <- 1:length(dim(p.xy))
            m1.x <- margin.table(p.xy, d1.x[-pos.y]) 
            if (any(abs(m1.x - p.x) > tol) & warn) 
                warning("The marginal distr. of the X variables \n 
                        in tab.xy is not equal to tab.x")
            d2.x <- 1:length(dim(p.xz))
            m2.x <- margin.table(p.xz, d2.x[-pos.z])
            if (any(abs(m2.x - p.x) > tol) & warn) 
                warning("The marginal distr. of the X variables \n 
                        in tab.xz is not equal to tab.x")
            if (any(abs(m1.x - m2.x) > tol) & warn) 
                warning("The marginal distr. of the X variables \n 
                        in tab.xy and in tab.xz are not equal")    
        }
        #
        # derives bounds for cells in YxZ in the unconditional case
        ll <- outer(p.y, p.z, FUN = "+") - 1
        low <- pmax(matrix(0, nrow(ll), ncol(ll)), ll)
        upper <- outer(p.y, p.z, FUN = "pmin")
        dimnames(low) <- dimnames(upper) <- list(names(p.y), 
                                                 names(p.z))
        class(low) <- class(upper) <- "table"
        res.0 <- list(low.u = low, up.u = upper)
        
        # derives bounds for cells in YxZ in the conditioning on
        # the joint distr. of available Xs
        # tables are re-arranged in a two-way contingency table
        if (length(dim(p.x)) == 1) 
            xx.0 <- cbind(p.x)
        else xx.0 <- ftable(p.x, row.vars = 1:length(dim(p.x)))
        xx.y <- ftable(p.xy, col.vars = length(dim(p.xy)))
        xx.z <- ftable(p.xz, col.vars = length(dim(p.xz)))
        
        # margins of xx.y and xx.z aligned to have same marginal distr.
        # of xx.0 without changing that of Y and Z
        # IPF is used, empty cells can assume a positive value
        # to avoid incoherencies
        
        if(align.margins){
            tot.x <- rowSums(xx.0)
            tot.y <- colSums(xx.y)
            tot.z <- colSums(xx.z)
            
            xx.y[xx.y==0] <- 1e-06 # constant added to empty cells
            xx.z[xx.z==0] <- 1e-06 # constant added to empty cells
            axy <- mipfp::Estimate(seed = xx.y, 
                                   target.list = list(1,2), 
                                   target.data = list(tot.x, tot.y)
                                   )
            axz <- mipfp::Estimate(seed = xx.z, 
                                   target.list = list(1,2), 
                                   target.data = list(tot.x, tot.z)
                                   )
            xx.y <- axy$p.hat
            xx.z <- axz$p.hat
        }
        # conditional rel. freq Y|X and Z|X
        ygxx <- prop.table(xx.y, margin = 1)
        zgxx <- prop.table(xx.z, margin = 1)
        
        # start computing conditional bounds for cells YxZ 
        # for each category of joint distr. of Xs
        H <- nrow(xx.0)
        out.CIA <- out.low <- out.up <- array(0, dim = c(ncol(xx.y), 
                                                         ncol(xx.z), H))
        for (h in 1:H) {
            out.CIA[, , h] <- outer(ygxx[h, ], zgxx[h, ], FUN = "*") * 
                xx.0[h, ]
            thetas <- outer(ygxx[h, ], zgxx[h, ], FUN = "+") - 1
            # thetas[thetas<0] <- 0
            ll <- pmax(thetas, matrix(0, nrow = nrow(thetas), 
                                      ncol = ncol(thetas)))
            uu <- outer(ygxx[h, ], zgxx[h, ], FUN = "pmin")
            out.low[, , h] <- ll * xx.0[h, ]
            out.up[, , h] <- uu * xx.0[h, ]
        }
        # computes expected values for CIA and bounds
        fine.CIA <- apply(out.CIA, c(1, 2), sum, na.rm = TRUE)
        fine.low <- apply(out.low, c(1, 2), sum, na.rm = TRUE)
        fine.up <- apply(out.up, c(1, 2), sum, na.rm = TRUE)
        l.y <- attr(ygxx, "col.vars")[[1]]
        l.z <- attr(zgxx, "col.vars")[[1]]
        
        class(fine.CIA) <- class(fine.low) <- class(fine.up) <- "table"
        dimnames(fine.CIA) <- dimnames(fine.low) <- dimnames(fine.up) <- list(l.y, 
                                                                              l.z)
        # average width of conditional and unconditional uncertainty bounds
        vet.unc <- c(av.u = mean(c(upper - low)), av.cx = mean(c(fine.up - 
                                                                     fine.low)))
        # output
        res.1 <- list(CIA = fine.CIA, low.cx = fine.low, up.cx = fine.up)
        res.2 <- list(uncertainty = vet.unc)
        if (print.f == "tables") {
            out <- c(res.0, res.1, res.2)
        }
        else if (print.f == "data.frame") {
            dataf <- data.frame(res.0$low.u)
            labdf <- c(lab.y, lab.z, "low.u")
            colnames(dataf) <- labdf
            dataf$low.cx <- c(res.1$low.cx)
            dataf$CIA <- c(res.1$CIA)
            dataf$up.cx <- c(res.1$up.cx)
            dataf$up.u <- c(res.0$up.u)
            out <- c(list(bounds = dataf), res.2)
        }
    }
    out
}
