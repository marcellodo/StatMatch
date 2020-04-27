kproto.dist <- function(x, y=NULL, lambda=NULL){
######################
    varCat1.fun <- function(x){
        px <- prop.table(table(x))
        1 - sum(px^2)
    }
##########################    
    #
    dx.num <- sapply(x, is.numeric)
    dx.cat <- sapply(x, is.factor)
    
    # checks
    if(!any(dx.num)) stop ("There are no numeric variables in the input dataframe")
    if(!any(dx.cat)) stop ("There are no categorical variables in the input dataframe")
    
    nx <- nrow(x)
    # separate variables
    xx.num <- x[,dx.num]
    xx.cat <- x[,dx.cat]
    if(is.null(y)){
        xx.num.sc <- yy.num.sc <- scale(x=as.matrix(xx.num),
                                  center = FALSE, 
                                  scale = TRUE)
        yy.cat <- xx.cat
    }
    
    if(!is.null(y)){
        yy.num <- y[, sapply(y, is.numeric)]
        xxyy <- as.matrix(rbind(xx.num, yy.num))
        xxyy.sc <- scale(xxyy, center = FALSE,
                         scale = TRUE)
        xx.num.sc <- xxyy.sc[1:nx, ]
        yy.num.sc <- xxyy.sc[-(1:nx), ]

        yy.cat <- y[, sapply(y, is.factor)]
    }
    

    # if(is.null(lambda)){
    #         var.cat <- sapply(rbind(xx.cat, yy.cat), varCat1.fun)
    #     }
    #     else{
    #         var.num <- sapply(xx.num, var)
    #         var.cat <- sapply(xx.cat, varCat1.fun)
    #     }
    #     #lambda <- mean(var.num)/mean(var.cat)
    #     
    # if( lambda==0 | is.infinite(lambda)) stop("")
    cat("The value of lambda is: ", lambda, fill=T)
    #

    # calculate distance
    dd.num <- proxy::dist(x = xx.num.sc, y = yy.num.sc, 
                          method = "Euclidean")
    xx.cat.dum <- fact2dummy(data=xx.cat, all=TRUE)
    yy.cat.dum <- fact2dummy(data=yy.cat, all=TRUE)
    dd.cat <- proxy::dist(x = xx.cat.dum, y = yy.cat.dum, 
                          method = "Simple matching")
    dd.num^2 + lambda * dd.cat
    
    
    
}