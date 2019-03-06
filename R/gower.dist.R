gower.dist<- function(data.x, data.y = data.x, rngs = NULL, KR.corr = TRUE,
                      var.weights = NULL){

####################    
        gower.fcn <- function(x, y, rng = NULL, KR.corr = TRUE) {
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
            if (is.null(rng) || is.na(rng)) rng <- max(x, y, na.rm=TRUE) - min(x, y, na.rm=TRUE)
            
            if(rng==0) dd <- matrix(0, nx, ny)
            else dd <- abs(outer(X = x, Y = y, FUN = "-"))/rng
            delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0
        }
        list(dist = dd, delta = delta)
    }

######## END  gower.fcn() ###################       

    if (is.null(dim(data.x)) && is.null(dim(data.y))) {
        out.gow <- gower.fcn(x = data.x, y = data.y, rng = rngs, 
                             KR.corr = KR.corr)
        out <- (out.gow$dist * out.gow$delta)/out.gow$delta
    }
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
              rng = rng.k, KR.corr = KR.corr)
            n <- out.gow$dist * out.gow$delta * w.k
            
            n[is.na(n)] <- 0
            num <- num + n
            d <- out.gow$delta * w.k
            d[is.na(d)] <- 0
            den <- den + d
        }
        out <- num/den
    }
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
              KR.corr = KR.corr)
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
