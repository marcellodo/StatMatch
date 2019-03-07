NND.hotdeck <- function (data.rec, data.don, match.vars, don.class = NULL, 
          dist.fun = "Manhattan", constrained = FALSE, constr.alg = "Hungarian", 
          k=1, keep.t = FALSE, ...) 
{
    p <- length(match.vars)
    if (!is.null(dim(data.rec))) {
        nr <- nrow(data.rec)
        r.lab <- row.names(data.rec)
    }
    else {
        nr <- length(data.rec)
        r.lab <- names(data.rec)
    }
    if (!is.null(dim(data.don))) {
        nd <- nrow(data.don)
        d.lab <- row.names(data.don)
    }
    else {
        nd <- length(data.don)
        d.lab <- names(data.don)
    }
    if (is.null(r.lab)) 
        r.lab <- paste("rec", 1:nr, sep = "=")
    else r.lab <- paste("rec", r.lab, sep = "=")
    row.names(data.rec) <- r.lab
    if (is.null(d.lab)) 
        d.lab <- paste("don", 1:nd, sep = "=")
    else d.lab <- paste("don", d.lab, sep = "=")
    row.names(data.don) <- d.lab
    if (!is.null(match.vars)) {
        if (dist.fun == "Euclidean" || dist.fun == "euclidean" || 
            dist.fun == "Manhattan" || dist.fun == "Mahalanobis" || 
            dist.fun == "mahalanobis" || dist.fun == "manhattan" || 
            dist.fun == "minimax" || dist.fun == "MiniMax" || 
            dist.fun == "Minimax") {
            cat("Warning: The ", dist.fun, " distance is being used", 
                fill = TRUE)
            cat("All the categorical matching variables in rec and don \n data.frames, if present are recoded into dummies", 
                fill = TRUE)
        }
        if (dist.fun == "exact" || dist.fun == "exact matching") {
            cat("Warning: the exact matching distance is being used", 
                fill = TRUE)
            cat("all the matching variables in rec and don are converted \n to character variables and are treated as categorical nominal", 
                fill = TRUE)
        }
    }
    NND.hd <- function(rec, don, dfun = "Manhattan", constr = FALSE, 
                       c.alg = NULL, ...) {
        x.rec <- rec
        x.don <- don
        p <- ncol(rec)
        nr <- nrow(x.rec)
        nd <- nrow(x.don)
        if (nr > nd) 
            cat("Warning: the number of donors is less than \n the number of recipients", 
                fill = TRUE)
        r.lab <- rownames(x.rec)
        if (is.null(r.lab)) 
            r.lab <- paste("rec", 1:nr, sep = "=")
        d.lab <- rownames(x.don)
        if (is.null(d.lab)) 
            d.lab <- paste("don", 1:nr, sep = "=")
        if (dfun == "Euclidean" || dfun == "euclidean" || dfun == 
            "Manhattan" || dfun == "manhattan") {
            if (is.data.frame(x.rec)) 
                x.rec <- fact2dummy(x.rec, all = TRUE)
            if (is.data.frame(x.don)) 
                x.don <- fact2dummy(x.don, all = TRUE)
            mdist <- proxy::dist(x = x.rec, y = x.don, method = dfun, 
                          ...)
        }
        else if (dfun == "Mahalanobis" || dfun == "mahalanobis") {
            if (is.data.frame(x.rec)) 
                x.rec <- fact2dummy(x.rec, all = TRUE)
            if (is.data.frame(x.don)) 
                x.don <- fact2dummy(x.don, all = TRUE)
            mdist <- mahalanobis.dist(data.x = x.rec, data.y = x.don, 
                                      ...)
        }
        else if (dfun == "minimax" || dfun == "MiniMax" || dfun == 
                 "Minimax") {
            if (is.data.frame(x.rec)) 
                x.rec <- fact2dummy(x.rec, all = TRUE)
            if (is.data.frame(x.don)) 
                x.don <- fact2dummy(x.don, all = TRUE)
            mdist <- maximum.dist(data.x = x.rec, data.y = x.don, 
                                  ...)
        }
        else if (dfun == "exact" || dfun == "exact matching") {
            dxr <- dim(x.rec)
            x.rec <- as.character(as.matrix(x.rec))
            dim(x.rec) <- dxr
            dxd <- dim(x.don)
            x.don <- as.character(as.matrix(x.don))
            dim(x.don) <- dxd
            xx <- data.frame(rbind(x.rec, x.don))
            x.rec <- xx[1:nr, ]
            x.don <- xx[-(1:nr), ]
            mdist <- gower.dist(data.x = x.rec, data.y = x.don, 
                                ...)
        }
        else if (dfun == "Gower" || dfun == "gower") {
            mdist <- gower.dist(data.x = x.rec, data.y = x.don, 
                                ...)
            mdist[is.nan(mdist)] <- 1
            mdist[is.na(mdist)] <- 1
        }
        else {
            mdist <- proxy::dist(x = x.rec, y = x.don, method = dfun, 
                          ...)
        }
        dimnames(mdist) <- list(r.lab, d.lab)
        if (!constr) {
            dist.rd <- numeric(nr)
            nad <- rep(NA, nr)
            don.lab <- numeric(nr)
            for (i in 1:nr) {
                vd <- mdist[i, ]
                if (all(is.na(vd))) 
                    stop("All the distances are NA \n this is due to the presence of NA in the \n matching variables")
                min.d <- min(vd)
                dist.rd[i] <- min.d
                appo <- d.lab[vd == min.d]
                nad[i] <- length(appo)
                if (length(appo) == 1) 
                    don.lab[i] <- appo
                else don.lab[i] <- sample(appo, 1)
            }
            rec.lab <- r.lab
        }
        if (constr && (c.alg == "lpSolve" || c.alg == "lpsolve")) {
            if (nr == nd){
                # appo <- lp.assign(cost.mat = mdist)
                r.sig <- rep("==", nr)
                r.rhs <- rep(1, nr)
                c.sig <- rep("<=", nd)
                c.rhs <- rep(k, nd)
                appo <- lpSolve::lp.transport(cost.mat = mdist, row.signs = r.sig, 
                                     row.rhs = r.rhs, col.signs = c.sig, col.rhs = c.rhs)
            } 
                
            else if (nr < nd) {
                r.sig <- rep("==", nr)
                r.rhs <- rep(1, nr)
                c.sig <- rep("<=", nd)
                c.rhs <- rep(k, nd)
                appo <- lpSolve::lp.transport(cost.mat = mdist, row.signs = r.sig, 
                                     row.rhs = r.rhs, col.signs = c.sig, col.rhs = c.rhs)
            }
            else if (nr > nd) {
                warning("There are more recipients than donors!")
                if(nr/nd > k) {
                    cat("The ratio (no. dons)/(no. recs) is ", nr/nd, fill=TRUE)
                    cat("it should be <= k")
                    stop()
                } 
                r.sig <- rep("==", nr)
                r.rhs <- rep(1, nr)
                c.sig <- rep("<=", nd)
                c.rhs <- rep(k, nd)
                appo <- lpSolve::lp.transport(cost.mat = mdist, row.signs = r.sig, 
                                     row.rhs = r.rhs, col.signs = c.sig, col.rhs = c.rhs)
            }
            sol <- appo$solution
            ss <- c(t(sol))
            cc <- c(t(col(sol)))
            dist.rd <- mdist[cbind(1:nr, cc[as.logical(ss)])]
            rec.lab <- r.lab
            don.lab <- d.lab[c(cc[as.logical(ss)])]
        }
        if (constr && (c.alg == "Hungarian" || c.alg == "hungarian")) {
            if (nr > nd) 
                stop("When using the Hungarian algorithm \n the no. of donors must be greater \n or equal than the no. of recipients")
            sol <- solve_LSAP(x = mdist, maximum = FALSE)
            rec.lab <- r.lab
            don.lab <- d.lab[as.integer(sol)]
            dist.rd <- mdist[cbind(rec.lab, don.lab)]
        }
        mtc.ids <- cbind(rec.id = rec.lab, don.id = don.lab)
        if (constr) 
            fine <- list(mtc.ids = mtc.ids, dist.rd = dist.rd, 
                         call = match.call())
        else fine <- list(mtc.ids = mtc.ids, dist.rd = dist.rd, 
                          noad = nad, call = match.call())
        fine
    }
    if (is.null(don.class)) {
        out <- NND.hd(rec = data.rec[, match.vars, drop = FALSE], 
                      don = data.don[, match.vars, drop = FALSE], dfun = dist.fun, 
                      constr = constrained, c.alg = constr.alg)
        mmm <- out$mtc.ids
        mmm <- substring(mmm, 5)
        if (is.null(rownames(data.rec)) && is.null(rownames(data.don))) 
            mtc.ids <- matrix(as.numeric(mmm), ncol = 2, byrow = TRUE)
        else mtc.ids <- mmm
        dimnames(mtc.ids) <- list(NULL, c("rec.id", "don.id"))
        dist.rd <- out$dist.rd
        if (!constrained) 
            noad <- out$noad
    }
    else {
        if (length(don.class) == 1) {
            l.rec <- split(data.rec[, match.vars, drop = FALSE], 
                           f = data.rec[, don.class], drop = TRUE)
            l.don <- split(data.don[, match.vars, drop = FALSE], 
                           f = data.don[, don.class], drop = TRUE)
        }
        else {
            l.rec <- split(data.rec[, match.vars, drop = FALSE], 
                           f = as.list(data.rec[, don.class]), drop = TRUE)
            l.don <- split(data.don[, match.vars, drop = FALSE], 
                           f = as.list(data.don[, don.class]), drop = TRUE)
        }
        tst <- which(!(names(l.rec) %in% names(l.don)))
        if (length(tst)) {
            cat("The following donation classes:", fill = TRUE)
            cat(names(l.rec)[tst], fill = TRUE)
            cat("in data.don are empty, i.e. there are no donors", 
                fill = TRUE)
            stop()
        }
        H <- length(l.rec)
        mtc.ids <- as.list(numeric(H))
        dist.rd <- as.list(numeric(H))
        if (!constrained) 
            noad <- as.list(numeric(H))
        for (h in 1:H) {
            lab.h <- names(l.rec)[h]
            if (keep.t) 
                cat("Selecting donors for donation class: ", 
                    lab.h, fill = TRUE)
            out <- NND.hd(rec = l.rec[[lab.h]], don = l.don[[lab.h]], 
                          dfun = dist.fun, constr = constrained, c.alg = constr.alg)
            mtc.ids[[h]] <- out$mtc.ids
            dist.rd[[h]] <- out$dist.rd
            if (!constrained) 
                noad[[h]] <- out$noad
        }
        mmm <- unlist(lapply(mtc.ids, t))
        mmm <- substring(mmm, 5)
        mtc.ids <- matrix(mmm, ncol = 2, byrow = TRUE)
        if (is.null(rownames(data.rec)) && is.null(rownames(data.don))) 
            mtc.ids <- matrix(as.numeric(mmm), ncol = 2, byrow = TRUE)
        dimnames(mtc.ids) <- list(NULL, c("rec.id", "don.id"))
        dist.rd <- unlist(dist.rd)
        if (!constrained) 
            noad <- unlist(noad)
    }
    if (constrained) 
        end.out <- list(mtc.ids = mtc.ids, dist.rd = dist.rd, 
                        call = match.call())
    else end.out <- list(mtc.ids = mtc.ids, dist.rd = dist.rd, 
                         noad = noad, call = match.call())
    end.out
}
