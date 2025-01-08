rho.bounds.pred <-
    function (data.rec, data.don, match.vars, y.rec, z.don,
              pred="lm", w.rec=NULL, w.don=NULL, 
              out.pred = FALSE, ...)
    {

        n.A <- nrow(data.rec)
        n.B <- nrow(data.don)

        data.rec <- data.frame(data.rec)
        data.don <- data.frame(data.don)

        fx <- paste(match.vars, collapse = "+")

        # fy <- paste(y.lab, ".", sep="~")
        # fz <- paste(z.lab, ".", sep="~")
        fy <- paste(y.rec, fx, sep = "~")
        fz <- paste(z.don, fx, sep = "~")
        #
        # linear regression predictions
        if(pred == "lm"){
            fit.y <- lm(as.formula(fy), data = data.rec )
            py.rec <- predict(fit.y)
            py.don <- predict(fit.y, newdata = data.don)

            fit.z <- lm(as.formula(fz), data = data.don)
            pz.don <- predict(fit.z)
            pz.rec <- predict(fit.z, newdata = data.rec)
        }
        #
        # prediction with robust regression
        if(pred == "roblm"){
            fit.y <- MASS::rlm(as.formula(fy), data = data.rec, ...)
            py.rec <- predict(fit.y)
            py.don <- predict(fit.y, newdata = data.don)

            fit.z <- MASS::rlm(as.formula(fz), data = data.don, ...)
            pz.don <- predict(fit.z)
            pz.rec <- predict(fit.z, newdata = data.rec)
        }
        #
        # lasso feature selection and regression predictions
        if(pred == "lasso"){
            x.rec <- data.rec[, match.vars, drop=FALSE]
            x.don <- data.don[, match.vars, drop=FALSE]
            xtype <- sapply(x.rec, class)
            if(any(xtype == "factor")){
                x.rec <- glmnet::makeX(x.rec)
                x.don <- glmnet::makeX(x.don)
            }
            # y get lambda.1se by cross-validation with 10 folds
            # and use it for predictions
            fit.y <- glmnet::cv.glmnet(x = x.rec, y = data.rec[, y.rec], ...)
            py.rec <- c(predict(object = fit.y, newx = x.rec, s = "lambda.1se"))
            py.don <- c(predict(object = fit.y, newx = x.don, s = "lambda.1se"))

            fit.z <- glmnet::cv.glmnet(x = x.don, y = data.don[, z.don], ...)
            pz.don <- c(predict(object = fit.z, newx = x.don, s = "lambda.1se"))
            pz.rec <- c(predict(object = fit.z, newx = x.rec, s = "lambda.1se"))

        }
        #
        # random forest predictions
        if(pred=="randomforest" | pred=="randomForest" | pred=="rf" | pred=="rF"){
            fit.y <- randomForest::randomForest(as.formula(fy), data=data.rec, ...)
            py.rec <- predict(object = fit.y, newdata = data.rec)
            py.don <- predict(object = fit.y, newdata = data.don)

            fit.z <- randomForest::randomForest(as.formula(fz), data=data.don, ...)
            pz.don <- predict(object = fit.z, newdata = data.don)
            pz.rec <- predict(object = fit.z, newdata = data.rec)
        }
        ###############################################################
        # estimate correlations between predicted and obs (sqrt(R^2))
        c.yy <- cor(py.rec, data.rec[ , y.rec])
        c.zz <- cor(pz.don, data.don[ , z.don])

        # puts Y, Z and the corresponding predictions in data.frames
        labpyz <- paste("pred", c(y.rec, z.don), sep=".")

        # new recipient
        if(is.null(w.rec))  {
            new.rec <- data.frame(py=c(py.rec), pz=c(pz.rec), y.rec=data.rec[,y.rec])
            colnames(new.rec) <- c(labpyz, y.rec)
        }
        else{
            new.rec <- data.frame(py=c(py.rec), pz=c(pz.rec), data.rec[,c(y.rec, w.rec)])
            colnames(new.rec) <- c(labpyz, y.rec, w.rec)
        }
        # new donor

        if(is.null(w.don))  {
            new.don <- data.frame(py=c(py.don), pz=c(pz.don), data.don[ , z.don])
            colnames(new.don) <- c(labpyz, z.don)
        }
        else{
            new.don <- cbind(py=c(py.don), pz=c(pz.don), data.don[ , c(z.don, w.don)])
            colnames(new.don) <- c(labpyz, z.don, w.don)
        }
        # cat(sapply(new.rec, class), fill=T)
        # cat(sapply(new.don, class), fill=T)
        #
        # applies function rho.bounds
        if( is.null(w.rec) & is.null(w.don)) {
            out <- rho.bounds(data.rec=new.rec, data.don=new.don,
                         match.vars = labpyz, y.rec=y.rec, z.don=z.don)
        }
        else{
            out <- rho.bounds(data.rec=new.rec, data.don=new.don,
                         match.vars = labpyz, y.rec=y.rec, z.don=z.don,
                         w.rec=w.rec, w.don=w.don)

        }
        ####################################################
        # output
        fine <- list(corr=c(py.y=c.yy, pz.z=c.zz), bounds=out)
        
        if(out.pred){
           up.rec <- data.frame(data.rec[ ,c(match.vars, y.rec, w.rec)], 
                                c(py.rec), c(pz.rec))
           colnames(up.rec) <- c(match.vars, y.rec, w.rec, 
                                        paste("pred", c(y.rec, z.don), sep="."))
           
           up.don <- data.frame(data.don[ ,c(match.vars, z.don, w.don)], 
                                c(py.don), pz=c(pz.don))
           colnames(up.don) <- c(match.vars, z.don, w.don, 
                                        paste("pred", c(y.rec, z.don), sep="."))
           ld <- list(up.rec = up.rec, up.don = up.don)
           fine <- c(ld, fine)
        }
        
        fine
        
    }

