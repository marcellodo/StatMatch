Fbounds.pred <-
    function (data.rec, data.don, match.vars, y.rec, z.don, 
              pred="multinom", w.rec=NULL, w.don=NULL, 
              type.pred = "random", out.pred = FALSE,
              ...)
    {
#######################################################
        acc.fcn <- function(tt){
            sum(diag(prop.table(tt)))
        }
#######################################################
        get.class <- function(pp, how="rnd"){
            n <- length(pp)
            if(n == 1){
                n <- 2
                pp <- c(pp, 1-pp)
            }
            if(how == "rnd") out <- sample(x = 1:n, size = 1, prob = pp)
            if(how == "mv") out <- which.max(pp)    

            out
        }
#######################################################
        link2prob <- function(x){
            exp(x)/sum(exp(x))
        }
#######################################################
        # ck.cross <- function(xf, probs, p.xf, p.rob){
        #     x <- as.integer(xf)
        #     J <- length(p.xf)
        #     for(j in 1:J){
        #         probs[x == p.xf[j], p.rob[j]] <- 0
        #     }
        #     # rescale probs
        #     t(apply(probs, 1, function(x){x/sum(x)}))
        # }
########################################################
        n.A <- nrow(data.rec)
        n.B <- nrow(data.don)
        
        data.rec <- data.frame(data.rec)
        data.don <- data.frame(data.don)
        #####################################################
        # 0) get conditional probs of both Y and Z
        fx <- paste(match.vars, collapse="+")
        fy <- paste(y.rec, fx, sep="~")
        fz <- paste(z.don, fx, sep="~")
        
        # 0.a) predictions with multinomial regression
        if(pred=="multinom" | pred=="multinomial"){
            
            fit.y <- nnet::multinom(as.formula(fy), data = data.rec)
            ppy.rec <- predict(fit.y, type = "prob") 
            ppy.don <- predict(fit.y, newdata = data.don, type = "prob")
            
            fit.z <- nnet::multinom(as.formula(fz), data = data.don)
            ppz.don <- predict(fit.z, type = "prob") 
            ppz.rec <- predict(fit.z, newdata = data.rec, type = "prob")
            
            if(is.null(dim(ppy.rec))){
                ppy.rec <- cbind(ppy.rec, 1 - ppy.rec)
                ppy.don <- cbind(ppy.don, 1 - ppy.don)
            }
            
            if(is.null(dim(ppz.don))){
                ppz.rec <- cbind(ppz.rec, 1 - ppz.rec)
                ppz.don <- cbind(ppz.don, 1 - ppz.don)
            }
        }
# 
        # 0.b) predictions with multinomial regression with lasso selection
        # of predictors
        if(pred=="lasso"){
            x.rec <- data.rec[, match.vars, drop = FALSE]
            x.don <- data.don[, match.vars, drop = FALSE]
            xtype <- sapply(x.rec, class)
            if(any(xtype == "factor")){
                x.rec <- glmnet::makeX(x.rec)
                x.don <- glmnet::makeX(x.don)
            }
            # y get lambda.1se by cross-validation with 10 folds 
            # and use it for predictions
            fit.y <- glmnet::cv.glmnet(x = x.rec, y = data.rec[, y.rec], 
                                       family = "multinomial", 
                                       type.multinomial = "grouped")
            py.rec <- predict(object = fit.y, newx = x.rec, 
                              s = "lambda.1se",
                              type="link")
            ppy.rec <- t(apply(py.rec[ , ,1], 1, link2prob))
            
            py.don <- predict(object = fit.y, newx = x.don, 
                              s = "lambda.1se",
                              type="link")
            ppy.don <- t(apply(py.don[ , ,1], 1, link2prob))
            
            #
            fit.z <- glmnet::cv.glmnet(x = x.don, y = data.don[, z.don], 
                                       family = "multinomial",
                                       type.multinomial = "grouped")
            pz.don <- predict(object = fit.z, newx = x.don, 
                              s = "lambda.1se",
                              type="link")
            ppz.don <- t(apply(pz.don[ , ,1], 1, link2prob))
            
            pz.rec <- predict(object = fit.z, newx = x.rec, s = "lambda.1se", 
                                type="link")
            ppz.rec <- t(apply(pz.rec[ , ,1], 1, link2prob))
            
        }
        #
        # 0.c) predictions with random forest
        if(pred=="randomforest" | pred=="randomForest" | pred=="rf" | pred=="rF"){
            fit.y <- randomForest::randomForest(y = data.rec[ ,y.rec ], 
                                                x = data.rec[ ,match.vars], 
                                                replace = FALSE)
            ppy.rec <- predict(object = fit.y, newdata = data.rec, type = "prob")
            ppy.don <- predict(object = fit.y, newdata = data.don, type = "prob")
            
            fit.z <- randomForest::randomForest(y = data.don[ ,z.don], 
                                                x = data.don[ ,match.vars],
                                                replace = FALSE)
            ppz.don <- predict(object = fit.z, newdata = data.don, type = "prob")
            ppz.rec <- predict(object = fit.z, newdata = data.rec, type = "prob")
        }
        #
        # 0.d) predictions with naive bayes
        if(pred=="naivebayes" | pred=="naiveBayes" | pred=="nb" | pred=="nB"){
            # x.rec <- data.rec[, match.vars, drop = FALSE]
            # xtype <- sapply(x.rec, class)
            # if(!all(xtype == "factor")) stop("the match.vars are not all factors \n
            #                                  naive Bayes classifier requires factor predictors")
            fit.y <- naivebayes::naive_bayes(as.formula(fy), data=data.rec, 
                                             usekernel = TRUE, ...)
            ppy.rec <- predict(object = fit.y, newdata = data.rec, type = "prob")
            ppy.don <- predict(object = fit.y, newdata = data.don, type = "prob")
            
            fit.z <- naivebayes::naive_bayes(as.formula(fz), data=data.don, 
                                             usekernel = TRUE, ...)
            ppz.don <- predict(object = fit.z, newdata = data.don, type = "prob")
            ppz.rec <- predict(object = fit.z, newdata = data.rec, type = "prob")
        }
        # END : # get estimated conditional probs (scores)
        # ###################################################################
        #
        # from estimated probs to predicted class 
        yy <- data.rec[ , y.rec]
        if(is.factor(yy))  lev.y <- levels(yy)
        else lev.y <- sort(unique(yy)) # integer
        
        nl.y <- length(lev.y)
        
        zz <- data.don[ , z.don]
        if(is.factor(zz)) lev.z <- levels(zz)
        else lev.z <- sort(unique(zz))
        nl.z <- length(lev.z)
        
        if(type.pred == "random" | type.pred == "rnd"){
            py.rec <- factor(apply(ppy.rec, 1, get.class, how="rnd"), 
                             levels = 1:nl.y)    
            pz.don <- factor(apply(ppz.don, 1, get.class, how="rnd"), 
                             levels = 1:nl.z)
            pz.rec <- factor(apply(ppz.rec, 1, get.class, how="rnd"), 
                             levels = 1:nl.z)
            py.don <- factor(apply(ppy.don, 1, get.class, how="rnd"), 
                             levels = 1:nl.y)
        }
        if(type.pred == "mostvoted" | type.pred == "mv"){
            py.rec <- factor(apply(ppy.rec, 1, get.class, how="mv"), 
                             levels = 1:nl.y)    
            pz.don <- factor(apply(ppz.don, 1, get.class, how="mv"), 
                             levels = 1:nl.z)
            pz.rec <- factor(apply(ppz.rec, 1, get.class, how="mv"), 
                             levels = 1:nl.z)
            py.don <- factor(apply(ppy.don, 1, get.class, how="mv"), 
                             levels = 1:nl.y)
        }
        
        
        # if(is.null(ini.yz)){
        #     pz.rec <- factor(apply(ppz.rec, 1, rnd.class), levels = 1:nl.z)
        #     py.don <- factor(apply(ppy.don, 1, rnd.class), levels = 1:nl.y)
        # }
        # else{
        #     pos.0 <- which(ini.yz == 0, arr.ind = TRUE)
        #     pps.rec <- ck.cross(xf = yy, probs = ppz.rec, 
        #                         p.xf = pos.0[ ,1], p.rob = pos.0[ ,2])
        #     pz.rec <- factor(apply(pps.rec, 1, rnd.class), levels = 1:nl.z)
        #     
        #     pps.don <- ck.cross(xf = zz, probs = ppy.don, 
        #                         p.xf = pos.0[ ,2], p.rob = pos.0[ ,1])
        #     py.don <- factor(apply(pps.don, 1, rnd.class), levels = 1:nl.y)
        # }
        # 
        # puts Y, Z and the corresponding predictions in data.frames
        # new recipient 
        labpyz <- paste("pred", c(y.rec, z.don), sep=".")
        if(is.null(w.rec))  {
            new.rec <- data.frame(py = c(py.rec), 
                                  pz = c(pz.rec), 
                                  y.rec = yy)
            colnames(new.rec) <- c(labpyz, y.rec)
        }
        else{
            new.rec <- data.frame(py = c(py.rec), 
                                  pz = c(pz.rec), 
                                  data.rec[ ,c(y.rec, w.rec)])
            colnames(new.rec) <- c(labpyz, y.rec, w.rec)
        }
        
        # new donor
        if(is.null(w.don))  {
            new.don <- data.frame(py = c(py.don), 
                                  pz = c(pz.don), 
                                  z.don = zz)
            colnames(new.don) <- c(labpyz, z.don)
        }
        else{
            new.don <- cbind(py = c(py.don), 
                             pz = c(pz.don), 
                             data.don[ , c(z.don, w.don)])
            colnames(new.don) <- c(labpyz, z.don, w.don)
        }
        
        # estimate tables to get the Frechet.bounds.cat
        # pooling the samples to estimate X distr
        fxx <- paste("~", paste(labpyz, collapse="+"))
        fxy <- paste(fxx, y.rec, sep="+")
        fxz <- paste(fxx, z.don, sep="+")
        
        if( is.null(w.rec) & is.null(w.don)) {
            l.A <- n.A / (n.A + n.B)
            # estimate input tables for function Frechet.bounds.cat
            tx.rec <- xtabs(as.formula(fxx), data=new.rec)
            tx.don <- xtabs(as.formula(fxx), data=new.don)
            txx <- tx.rec * l.A + tx.don * (1 - l.A)

            txy <- xtabs(as.formula(paste(fxx, y.rec, sep="+")), data=new.rec)
            txz <- xtabs(as.formula(paste(fxx, z.don, sep="+")), data=new.don)
        }
        else{
            # rescale weights by pooling the samples
           
            w.A <- data.rec[, w.rec]
            w.B <- data.don[, w.don]

            d.A <- 1 + (sd(w.A)/mean(w.A))^2
            d.B <- 1 + (sd(w.B)/mean(w.B))^2
            l.A <- (n.A/d.A) / (n.A/d.A + n.B/d.B)
            # l.B <- 1 - l.A
            new.rec[ , w.rec] <- l.A * w.A
            new.don[ , w.don] <- (1 - l.A) * w.B

            ## estimate input tables for function Frechet.bounds.cat
            tx.rec <- xtabs(as.formula(paste(w.rec, fxx)), data=new.rec)
            tx.don <- xtabs(as.formula(paste(w.don, fxx)), data=new.don)
            txx <- tx.rec + tx.don
            
            fxy <- paste(w.rec, fxy, sep="")
            txy <- xtabs(as.formula(fxy), data=new.rec)
            
            fxz <- paste(w.don, fxz, sep="")
            txz <- xtabs(as.formula(fxz), data=new.don)
        }

        # applies function Frechet.bounds.cat

        ######################################################
        # diagnostic: accuracy
        acc.y <- acc.fcn(margin.table(txy, margin = c(1, 3)))
        acc.z <- acc.fcn(margin.table(txz, margin = c(2, 3)))
        
        #  get joint distribution of predictions
       
        out <- Frechet.bounds.cat(tab.x = prop.table(txx), 
                                  tab.xy = prop.table(txy), 
                                  tab.xz = prop.table(txz), 
                                  align.margins = TRUE,
                                  print.f = "data.frame", ...)
        
        # prepare output
        if(out.pred){
           
            out.rec <- data.frame(ppy.rec, ppz.rec, py.rec, pz.rec,
                                  data.rec[, c(y.rec, match.vars, w.rec)])
            colnames(out.rec) <- c(paste("pr", y.rec, lev.y, sep = "."),
                                   paste("pr", z.don, lev.z, sep = "."), 
                                   paste("pred", c(y.rec, z.don), sep="."),
                                   y.rec, match.vars, w.rec)
              
            out.don <- data.frame(ppy.don, ppz.don, py.don, pz.don,
                                  data.don[, c(z.don, match.vars, w.don)])
            colnames(out.don) <- c(paste("prob", y.rec, lev.y, sep = "."),
                                   paste("prob", z.don, lev.z, sep = "."), 
                                   paste("pred", c(y.rec, z.don), sep="."),
                                   z.don, match.vars, w.don)
              
            linp <- list(up.rec = out.rec, up.don = out.don,
                           p.xx.ini = prop.table(txx), 
                           p.xy.ini = prop.table(txy), 
                           p.xz.ini = prop.table(txz),
                           accuracy=c(y = acc.y, z = acc.z))
        }
        else{
            linp <- list(p.xx.ini = prop.table(txx), 
                         p.xy.ini = prop.table(txy), 
                         p.xz.ini = prop.table(txz),
                         accuracy=c(y = acc.y, z = acc.z))
        }
        
        # final output
        c(linp, out)
    }

