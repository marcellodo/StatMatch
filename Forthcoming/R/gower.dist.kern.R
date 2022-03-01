gower.dist<- function(data.x, data.y = data.x, rngs = NULL, KR.corr = TRUE,
                      cat.cont=TRUE, var.weights = NULL, kern="kde0"){

#########################
# core function: Gower with vectors    
    
        gower.fcn <- function(x, y, rng = NULL, KR.corr = TRUE, catc=TRUE, krn="kde0") {
        nx <- length(x)
        ny <- length(y)
        cx <- class(x)
        cy <- class(y)
        delta <- matrix(1, nx, ny)
        if (!identical(cx, cy)) 
            stop("the x and y object are of different type")
        if (is.logical(x)) {
            dd <- abs(outer(X = x, Y = y, FUN = "-"))
            delta[outer(x == FALSE, y == FALSE, FUN = "&")] <- 0
            delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0
        }
        else if (is.character(x) || (is.factor(x) && !is.ordered(x))) {
            if (is.factor(x) && !identical(levels(x), levels(y))) 
                stop("x and y have different levels")
            dd <- 1 - outer(x, y, FUN = "==")
            delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0
        }
        else if (is.ordered(x)) {
            if (KR.corr) {
                x <- as.numeric(x)
                y <- as.numeric(y)
                if (is.null(rng) || is.na(rng)) 
                    rng <- max(x, y, na.rm = TRUE) - 1
                if(rng==0) {
                    dd <- matrix(0, nx, ny)
                    delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0
                }    
                else{
                    zx <- (x - 1)/rng
                    zy <- (y - 1)/rng
                    dd <- abs(outer(X = zx, Y = zy, FUN = "-"))/(max(zx, 
                                                                 zy, na.rm=TRUE) - min(zx, zy, na.rm=TRUE))
                    delta[outer(is.na(zx), is.na(zy), FUN = "|")] <- 0
                }
            }
            else {
                x <- as.numeric(x)
                y <- as.numeric(y)
                if (is.null(rng) || is.na(rng)) 
                    rng <- max(x, y, na.rm=TRUE) - 1
                if(rng==0) dd <- matrix(0, nx, ny)
                else dd <- abs(outer(X = x, Y = y, FUN = "-"))/rng
                delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0
            }
        }
        else {
            # categorization of continuous variables
            if(catc){
                dd <- abs(outer(X = x, Y = y, FUN = "-"))
                sd.y <- sd(y)
                
                A <- min(sd.y, IQR(y)/1.34)
                
                if(tolower(kern)=="kde0"){ 
                    h <- 1.06 * sd.y / (ny^(1/5))
                }
                else if(tolower(kern)=="kde1") {
                    h <- 1.06 * A / (ny^(1/5))
                }
                else if(tolower(kern)=="kde2"){ 
                    h <- 0.9 * A / (ny^(1/5))
                }
                dd[dd < h] <- 0
                dd[dd >= h] <- 1
                delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0
            }
            
            
            # standard treatment of continuous variables
            # Manhattan scaled by the range
            else{
                if (is.null(rng) || is.na(rng)) rng <- max(x, y, na.rm=TRUE) - min(x, y, na.rm=TRUE)
                
                if(rng==0) dd <- matrix(0, nx, ny)
                else dd <- abs(outer(X = x, Y = y, FUN = "-"))/rng
                delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0    
            }
            
            
        }
        list(dist = dd, delta = delta)
    }

######## END  gower.fcn() ###################       

################################################
# input data are single variables (vectors)
    if (is.null(dim(data.x)) && is.null(dim(data.y))) {
        out.gow <- gower.fcn(x = data.x, y = data.y, rng = rngs, 
                             KR.corr = KR.corr, catc=cat.cont, krn=kern)
        out <- (out.gow$dist * out.gow$delta)/out.gow$delta
    }
    # data.x is a single obs while data.y has more obs.  
    else if (is.null(dim(data.x)) && !is.null(dim(data.y))) {
        p <- ncol(data.y)
        if (length(data.x) != p) 
            stop("data.x should be of the same length of the no. of cols of data.y")
        num <- array(0, c(1, nrow(data.y)))
        den <- array(0, c(1, nrow(data.y)))
        if(is.null(var.weights)) var.weights <- rep(1, p)
        for (k in 1:p) {
            
            if (is.null(rngs)) 
                rng.k <- NULL
            else rng.k <- rngs[k]
            w.k <- var.weights[k]
            out.gow <- gower.fcn(x = data.x[, k], y = data.y[,k],
              rng = rng.k, KR.corr = KR.corr, catc=cat.cont, krn=kern)
            n <- out.gow$dist * out.gow$delta * w.k
            
            n[is.na(n)] <- 0
            num <- num + n
            d <- out.gow$delta * w.k
            d[is.na(d)] <- 0
            den <- den + d
        }
        out <- num/den
    }
# data.x and data.y have bot more than one observation        
    else {
        p <- ncol(data.y)
        if (ncol(data.x) != p) 
            stop("data.x and data.y must have the same no. of cols")
        num <- array(0, c(nrow(data.x), nrow(data.y)))
        den <- array(0, c(nrow(data.x), nrow(data.y)))
        if(is.null(var.weights)) var.weights <- rep(1, p)
        for (k in 1:p) {
            
            if (is.null(rngs)) 
                rng.k <- NULL
            else rng.k <- rngs[k]
            w.k <- var.weights[k]
            out.gow <- gower.fcn(x = data.x[, k], y = data.y[, k], rng = rng.k,
              KR.corr = KR.corr, catc=cat.cont, krn=kern)
            n <- out.gow$dist * out.gow$delta * w.k
            
            n[is.na(n)] <- 0
            num <- num + n
            d <- out.gow$delta * w.k
            d[is.na(d)] <- 0
            den <- den + d
        }
        out <- num/den
    }
    out
}    
