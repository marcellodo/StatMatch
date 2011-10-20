'comb.samples' <-
function(svy.A, svy.B, svy.C=NULL, y.lab, z.lab, form.x, estimation=NULL, ...)
{

    require(survey)
    data.A <- svy.A$variables
    y.lev <- levels(data.A[, y.lab])
    levels(data.A[, y.lab]) <- 1:nlevels(data.A[, y.lab])
    svy.A$variables <- data.A

    n.A <- nrow(data.A)
    w.A <- weights(svy.A)

    data.B <- svy.B$variables
    z.lev <- levels(data.B[, z.lab])
    levels(data.B[, z.lab]) <- 1:nlevels(data.B[, z.lab])
    svy.B$variables <- data.B

    n.B <- nrow(data.B)
    w.B <- weights(svy.B)

    X.A <- model.matrix(form.x, data=data.A)
    X.B <- model.matrix(form.x, data=data.B)

    form.y <- as.formula(paste("~", y.lab, "-1", collapse=""))
    form.z <- as.formula(paste("~", z.lab, "-1", collapse=""))

    Y.A <- model.matrix(form.y, data=data.A)
    Z.B <- model.matrix(form.z, data=data.B)

    QR.A <- qr(X.A * sqrt(w.A))
    beta.yx.A <- qr.coef(QR.A, Y.A*sqrt(w.A))
    beta.yx.A[is.na(beta.yx.A)] <- 0

    QR.B <- qr(X.B * sqrt(w.B))
    beta.zx.B <- qr.coef(QR.B, Z.B*sqrt(w.B))
    beta.zx.B[is.na(beta.zx.B)] <- 0

    XX.wA <- t(X.A) %*% (X.A*w.A)
    XX.wB <- t(X.B) %*% (X.B*w.B)
    gamma.p <- n.A/(n.A+n.B)
    XX.pool <- gamma.p*XX.wA + (1-gamma.p)*XX.wB
    YZ.CIA <- t(beta.yx.A) %*% XX.pool %*% beta.zx.B
    dimnames(YZ.CIA) <- list(y.lev, z.lev)
    
###########################
# use of auxiliary sample "svy.C"
    if(!is.null(svy.C)){
        data.C <- svy.C$variables
        y.lev.C <- levels(data.C[, y.lab])
        z.lev.C <- levels(data.C[, z.lab])
        #check
        if( all.equal(y.lev, y.lev.C) ) levels(data.C[, y.lab]) <- 1:nlevels(data.C[, y.lab])
        else stop("The levels of y.lab in svy.A and in svy.C do not match")
        if( all.equal(z.lev, z.lev.C) ) levels(data.C[, z.lab]) <- 1:nlevels(data.C[, z.lab])
        else stop("The levels of z.lab in svy.B and in svy.C do not match")
        svy.C$variables <- data.C
        n.C <- nrow(data.C)
        w.C <- weights(svy.C)

# Incomplete two-way Stratification
        if(estimation=="ITWS" || estimation=="i2ws" || estimation=="incomplete"){
            tot.y.A <- colSums(Y.A * w.A)
            tot.z.B <- colSums(Z.B * w.B)
            tot.yz <- c(tot.y.A, tot.z.B[-1])
            form.yz <- as.formula(paste("~",paste(y.lab, z.lab, sep="+"), "- 1", sep=""))
            cal.C <- calibrate(design=svy.C, formula=form.yz,
                        population=tot.yz, ...)
        }

# Synthetic two-way Stratification
        if(estimation=="STWS" || estimation=="s2ws" || estimation=="synthetic" ){
            X.C <- model.matrix(form.x, data=data.C)
            Y.C <- model.matrix(form.y, data=data.C)
            Z.C <- model.matrix(form.z, data=data.C)
            resY.C <- Y.C - (X.C %*% beta.yx.A)
            resZ.C <- Z.C - (X.C %*% beta.zx.B)
            c.y <- ncol(Y.C)
            c.z <- ncol(Z.C)
            new.YZ <- matrix(NA, nrow=n.C, ncol=(c.y*c.z))
            for(i in 1:n.C){
                m1 <- cbind(Y.C[i,]) %*% rbind(Z.C[i,])
                m2 <- cbind(resY.C[i,]) %*% rbind(resZ.C[i,])
                new.YZ[i, ] <- c(m1)-c(m2)
            }
            lab1 <- rep(colnames(Y.C), c.z)
            lab2 <- rep(colnames(Z.C), each=c.y)
            lab <- paste(lab1, lab2, sep="_")
            colnames(new.YZ) <- lab
            orig.vars <- colnames(svy.C$variables)
            svy.C$variables <- data.frame(svy.C$variables, new.YZ)

            vec.tot <- c(YZ.CIA)
            names(vec.tot) <- lab
            form.yz <- as.formula(paste("~",paste(lab, collapse="+"), "- 1", sep=""))
            cal.C <- calibrate(design=svy.C, formula=form.yz,
                            population=vec.tot, ...)
            cal.C$variables <-  cal.C$variables[,orig.vars]
        }
        ww.C <- weights(cal.C)
        f.yz <- paste("ww.C", paste(y.lab, z.lab, sep="+"), sep="~")
        YZ.noCIA <- xtabs(as.formula(f.yz), data=data.C)
        dimnames(YZ.noCIA) <- list(y.lev, z.lev)
        out <- list(yz.CIA=YZ.CIA, cal.C=cal.C, yz.est=YZ.noCIA, call=match.call())
    }
# output
    else out <- list(yz.CIA=YZ.CIA, call=match.call())
    out
}
