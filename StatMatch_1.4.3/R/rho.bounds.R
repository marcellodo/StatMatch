#' Title
#'
#' @param data.rec
#' @param data.don
#' @param match.vars
#' @param y.rec
#' @param z.don
#' @param w.rec
#' @param w.don
#'
#' @return
#' @export
#'
#' @examples
rho.bounds <-
    function (data.rec, data.don, match.vars, y.rec, z.don,
              w.rec=NULL, w.don=NULL)
    {
##############################################################
        # function to estimate bounds by rodgers & devol
        corr.bounds <- function(cxx, cxy, cxz){
            C <- sum( outer(X = c(cxy), Y = c(cxz), FUN = "*") * solve(cxx) )
            Dy <- 1 - sum( outer(X = c(cxy), Y = c(cxy), FUN = "*") * solve(cxx) )
            Dz <- 1 - sum( outer(X = c(cxz), Y = c(cxz), FUN = "*") * solve(cxx) )
            low.c0 <- C - sqrt(Dy * Dz)
            up.c0 <- C + sqrt(Dy * Dz)

            c(low = low.c0,
              midp =(low.c0 + up.c0) / 2,
              up = up.c0)
        }
############################################################
        n.A <- nrow(data.rec)
        n.B <- nrow(data.don)
        px <- length(match.vars)
        if(!is.null(w.rec) & !is.null(w.don)){
            w.A <- data.rec[, w.rec]
            w.B <- data.don[, w.don]
        }

        # treat predictors in rec
        xrec <- data.rec[ ,match.vars, drop = FALSE]
        xrec <- StatMatch::fact2dummy(data = xrec, all = FALSE)
        m.rec <- cbind(xrec, data.rec[, y.rec])
        colnames(m.rec) <- c(colnames(xrec), y.rec)

        # treat predictors in rec
        xdon <- data.don[ ,match.vars, drop = FALSE]
        xdon <- StatMatch::fact2dummy(data = xdon, all = FALSE)
        m.don <- cbind(xdon, data.don[, z.don])
        colnames(m.don) <- c(colnames(xdon), z.don)

        # info on Xs
        mtc.v <- colnames(xrec)
        p.x <- length(mtc.v)
        x.A <- matrix(c(m.rec[ ,mtc.v]), ncol = p.x)
        pos.x.A <- match(mtc.v, colnames(m.rec))
        x.B <- matrix(c(m.don[ ,mtc.v]), ncol = p.x)
        pos.x.B <- match(mtc.v, colnames(m.don))

        # preparing dependent variables Y and Z
        # (assumed to be continuous)
        y.A <- m.rec[ , y.rec]
        z.B <- m.don[ , z.don]

        # estimate parameters of the multivariate Gaussian distribution
        # with Moriarity and Scheuren method

        # estimates pieces of corr matrix without survey weights
        if(is.null(w.rec) & is.null(w.don)){
            # estimates for X
            if(px==1) S.x <- 1
            else S.x <- cor(rbind(x.B, x.A))

            # estimates for Y
            S.xy <- c(cor(x.A, y.A))

            # estimates for Z
            S.xz <- c(cor(x.B, z.B))
        }
        # estimates pieces of corr matrix with survey weights
        else{
            if(p.x==1) S.x <- 1
            else{
                # weighted estimation
                if(!is.null(w.rec) & !is.null(w.don)){
                    d.A <- 1 + (sd(w.A)/mean(w.A))^2 # deff weights rec
                    d.B <- 1 + (sd(w.B)/mean(w.B))^2 # deff weights don
                    l.A <- (n.A/d.A) / (n.A/d.A + n.B/d.B)
                    l.B <- 1 - l.A
                }

                oo <- cov.wt(x = rbind(x.A, x.B),
                             wt = c(l.A * w.A, l.B * w.B),
                             cor = TRUE)
                S.x <- oo$cor
            }

            # estimates for Y
            oo <- cov.wt(x = cbind(x.A, y.A), wt = w.A, cor = TRUE)
            S.xy <- c(oo$cor[nrow(oo$cor), -nrow(oo$cor)])

            # estimates for Z
            oo <- cov.wt(x = cbind(x.B, z.B), wt = w.B, cor = TRUE)
            S.xz <- c(oo$cor[nrow(oo$cor), -nrow(oo$cor)])

            rm(oo)
            gc()

        }

        # fills in the Var-Cov matrix
        p <- p.x+2
        pos.x <- 1:p.x
        pos.y <- p.x+1
        pos.z <- p.x+2

        cc <- diag(p)
        cc[pos.x, pos.x] <- S.x
        cc[pos.x, pos.y] <- cc[pos.y, pos.x] <- S.xy
        cc[pos.x, pos.z] <- cc[pos.z, pos.x] <- S.xz

        # estimation of S.yz
        # estimate rho,yz under CI
        if(p.x == 1){
            c.xy <- cc[pos.x, pos.y]
            c.xz <- cc[pos.x, pos.z]
            low.c <- c.xy * c.xz - sqrt((1 - c.xy^2)*(1 - c.xz^2))
            up.c <-  c.xy * c.xz + sqrt((1 - c.xy^2)*(1 - c.xz^2))

            out <- c(low = low.c, midp = c.xy * c.xz, up = up.c)
        }
        else{
            # apply rodgers & devol method
            bb1 <- corr.bounds(cxx=cc[pos.x, pos.x],
                               cxy=cc[pos.x, pos.y],
                               cxz=cc[pos.x, pos.z])

            # # estimate by grid search
            # eps <- 0.0001
            #
            # rr <- seq(-1, 1, eps)
            # k <- length(rr)
            # vdet <- rep(0,k)
            # for(i in 1:k){
            #     cc[pos.z, pos.y] <- cc[pos.y, pos.z] <- rr[i]
            #     vdet[i] <- det(cc)
            # }
            # cc.yz <- rr[vdet>=0]
            # low.c1 <- min(cc.yz)
            # up.c1 <- max(cc.yz)
            #
            # bbsim <- c(low=low.c1, midp=(low.c1+up.c1)/2, up=up.c1)
            #
            # out <- rbind(analytic=bb1, grid.search=bbsim)
            out <- bb1
        }
        out
    }
