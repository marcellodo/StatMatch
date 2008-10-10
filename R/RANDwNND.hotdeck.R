`RANDwNND.hotdeck` <-
function (data.rec, data.don, match.vars=NULL, don.class=NULL, dist.fun="Euclidean", cut.don="rot", k=NULL) 
{
	if(dist.fun!="Gower" || dist.fun!="exact" || dist.fun!="exact matching"){
		require(proxy)
	}
	
	p <- length(match.vars)
	if(!is.null(dim(data.rec))){
		nr <- nrow(data.rec)
		r.lab <- row.names(data.rec)
	}
	else{
		nr <- length(data.rec)
		r.lab <- names(data.rec)
	}
	if(!is.null(dim(data.don))){
		nd <- nrow(data.don)
		d.lab <- row.names(data.don)
	}
	else{
		nd <- length(data.don)
		d.lab <- names(data.don)
	}
		
	if(is.null(r.lab)) r.lab <- paste("rec", 1:nr, sep="=")
	else r.lab <- paste("rec", r.lab, sep="=")
	row.names(data.rec) <- r.lab
	
	if(is.null(d.lab)) d.lab <- paste("don", 1:nd, sep="=")
	else d.lab <- paste("don", d.lab, sep="=")
	row.names(data.don) <- d.lab

################
RANDwNND.hd <- function (rec, don, dfun="Euclidean", cut.don="rot", k=NULL)
{ 
	p <- ncol(rec)
	x.rec <- rec
	x.don <- don
		
	nr <- nrow(rec)
	nd <- nrow(don)
	if(nr>nd) cat("Warning: the no. do donors is less than the no.f recipients", fill=TRUE)

	r.lab <- rownames(rec)
	if(is.null(r.lab)) r.lab <- 1:nr
	d.lab <- rownames(don)
	if(is.null(d.lab)) d.lab <- 1:nd

# compute matrix of distances between obs. in x.don and obs. in x.rec
# function dist() in package "proxy" is used! 

    if(dfun=="Euclidean" || dfun=="Manhattan"){
        cat("The ", dfun, " distance is being used", fill=TRUE)
        cat("Warning: all the categorical variables in rec and don data.frame are recoded into dummies", fill=TRUE)
        x.rec <- fact2dummy(x.rec, all=FALSE)
        x.don <- fact2dummy(x.don, all=FALSE)
        mdist <- dist(x=x.rec, y=x.don, method=dfun)
    }
    else if(dfun=="exact" || dfun=="exact matching"){
        cat("The exact matching distance is being used", fill=TRUE)
        cat(" Warning: all the variables in rec and don data frame are converted to character variables and are treated as categorical nominal", fill=TRUE)
        dxr <- dim(x.rec)
        x.rec <- as.character(as.matrix(x.rec))
        dim(x.rec) <- dxr
        dxd <- dim(x.don)
        x.don <- as.character(as.matrix(x.don))
        dim(x.don) <- dxd
		    mdist <- gower.dist(data.x=x.rec, data.y=x.don)
    }
    else if(dfun=="Gower"){
        if(p==1 && is.factor(x.rec)) x.rec <- list(x.rec)
        if(p==1 && is.factor(x.don)) x.don <- list(x.don)
        mdist <- gower.dist(data.x=x.rec, data.y=x.don)
        mdist[is.nan(mdist)] <- 1 # NaN can occur when p=1 and x.rec and x.don is of type logical
        mdist[is.na(mdist)] <- 1 # NA can occur when p=1 and x.rec and x.don is of type logical
    }
    else mdist <- dist(x=x.rec, y=x.don, method=dfun)
    dimnames(mdist) <- list(r.lab, d.lab)
    
    if(cut.don=="rot"){
        k <- ceiling(sqrt(nd))
    		if(k==0) stop("k=sqrt(no. of dons) is equal to 0")  
    }
    else if(cut.don=="span"){ 
        if(k==0 || k>1) stop("When cut.don=span, k should be such that  0 < k <= 1") 
		    k <- ceiling(nd*k)
		    if(k>nd) k <- nd
		    if(k==0) stop("k=round(no. of dons x k) is equal to 0")
	  } 
    else if(cut.don=="exact") {
	    if(k==0 || k>nd) stop("When cut.don=exact, k should be such that  1 < k <= no. of dons") 
    }
	
    min.d <- numeric(nr)
	  max.d <- numeric(nr)
    sd.d <- numeric(nr)
    cut.d <- numeric(nr)
    dist.rd <- numeric(nr)
    nad <- rep(NA, nr)

    don.lab <- numeric(nr)
	  for(i in 1:nr){
		    vd <- mdist[i,]
        min.dist <- min(vd, na.rm=TRUE) # smallest distance recipient-donor
        min.d[i] <- min.dist
        max.d[i] <- max(vd, na.rm=TRUE)
		    sd.d[i] <- sd(vd, na.rm=TRUE)
        if(cut.don=="min"){
        short.vd <- vd[(vd==min.dist) & !is.na(vd)]
        appo <- d.lab[(vd==min.dist) & !is.na(vd)]
        dist.rd[i] <- min.dist
        cut.d[i] <- min.dist
		}	
        else if(cut.don=="k.dist"){
			if(k<min.dist) {
				cat("Warning: the value of k,", k, fill=TRUE)
				cat("is smaller than the minimum distance:", min.d, fill=TRUE)
			}			
			appo <- d.lab[(vd<=k) & !is.na(vd)]
			short.vd <- vd[(vd<=k) & !is.na(vd)]
            cut.d[i] <- k
        }
        else {
            appo <- d.lab[order(vd, na.last=NA)]
			if(length(appo)<k) kk <- length(appo)
			else kk <- k
			appo <- appo[1:kk]
			short.vd <- vd[order(vd, na.last=NA)][1:kk]
            cut.d[i] <- short.vd[kk]
        }

        nad[i] <- length(appo) # number of availabe donors
        if(length(appo)==0) {
            nad[i] <- 0
			don.lab[i] <- NA
			dist.rd[i] <- NA
            cat("Warning: there are no available donors for the", d.lab[i], "recipient!", fill=TRUE)
        }       
        else if(length(appo)==1){
			nad[i] <- 1
			don.lab[i] <- appo 
			dist.rd[i] <- short.vd
		}	
        else{
			nn.dd <- length(appo)
			nad[i] <- nn.dd
			choi <- sample(1:nn.dd, 1)
			don.lab[i] <- appo[choi]
			dist.rd[i] <- short.vd[choi]
		}	
    }       
    rec.lab <- r.lab

# output
    mtc.ids <- cbind(rec.id=rec.lab, don.id=don.lab)
	sum.dist <- cbind(min=min.d, max=max.d, sd=sd.d, cut=cut.d, dist.rd=dist.rd)
	list(mtc.ids=mtc.ids, sum.dist=sum.dist, noad=nad, call=match.call())
}
#########################
# RANDwNND.hd ends here
############################################	

	if(is.null(don.class)){ 
		if(is.null(match.vars)){
			don.lab <- sample(d.lab, nr, replace=TRUE)
			mtc.ids <- cbind(rec.id=r.lab, don.id=don.lab)
			noad <- rep(nd, nr)
			sss.dist <- NULL
		}
		else {
			if(p==1){
				REC <- data.frame(data.rec[, match.vars])
				names(REC) <- match.vars
				row.names(REC) <- r.lab
				DON <- data.frame(data.don[, match.vars])
				names(DON) <- match.vars
				row.names(DON) <-d.lab
			}
			else{
				REC <- data.rec[, match.vars]
				DON <- data.don[, match.vars]
			}
			out <- RANDwNND.hd(rec=REC, don=DON, dfun=dist.fun, cut.don=cut.don, k=k)
			mtc.ids <- out$mtc.ids
			noad <- out$noad
			sss.dist <- out$sum.dist
		}
		mmm <- substring(c(mtc.ids), 5)
		mtc.ids <- matrix(mmm, ncol=2)
		if(is.null(rownames(data.rec)) && is.null(rownames(data.don)))  mtc.ids <- matrix(as.numeric(mmm), ncol=2)
	}
	else{
		if(length(don.class)==1){
			l.r.lab <- split(r.lab, f=data.rec[ ,don.class])
			l.rec <- split(data.rec[ ,match.vars], f=data.rec[ ,don.class])
			l.d.lab <- split(d.lab, f=data.don[ ,don.class])
			l.don <- split(data.don[ ,match.vars], f=data.don[ ,don.class])
		}
		else{
			l.r.lab <- split(r.lab, f=as.list(data.rec[ ,don.class]))
			l.rec <- split(data.rec[ ,match.vars], f=as.list(data.rec[ ,don.class]))
			l.d.lab <- split(d.lab, f=as.list(data.don[ ,don.class]))
			l.don <- split(data.don[ ,match.vars], f=as.list(data.don[ ,don.class]))
		}
		if(length(l.rec)!=length(l.don)){
			cat("The no. of donation classes in recipient data is not equal to the no. of donation classes in donor data", fill=TRUE)
			stop("Possible reason: the variables used to classify units are not defined as factors or are factors with different levels") 	
		}
		
		nn.r <- unlist(lapply(l.r.lab, length))
		nn.d <- unlist(lapply(l.d.lab, length))
				
		if(sum( (nn.r>0) & (nn.d==0))>0) {
			stop("In some donation classes there are NO available donors. Please modify the definition of the donation classes")
		}	
		l.r.lab <- l.r.lab[nn.r>0]
		l.d.lab <- l.d.lab[nn.r>0]
		l.rec <- l.rec[nn.r>0]
		l.don <- l.don[nn.r>0]
				
		H <- length(l.rec)
		mtc.ids <- as.list(numeric(H))
		sum.dist <- as.list(numeric(H))
		noad <- as.list(numeric(H))
		
		for(h in 1:H){
			if(is.null(match.vars)){
				don.lab <- sample(l.d.lab[[h]], nn.r[[h]], replace=TRUE)
				mtc.ids[[h]] <- cbind(rec.id=l.r.lab[[h]], don.id=don.lab)
				sum.dist[[h]] <- NA
				noad[[h]] <- rep(nn.d[[h]], nn.r[[h]])
			}
			else{
				if(p==1){
					REC <- data.frame(l.rec[[h]])
					row.names(REC) <- l.r.lab[[h]]
					DON <- data.frame(l.don[[h]])
					row.names(DON) <- l.d.lab[[h]]
				}
				else{
					REC <- l.rec[[h]]
					DON <- l.don[[h]]
				}
				out <- RANDwNND.hd(rec=REC, don=DON, dfun=dist.fun, cut.don=cut.don, k=k)
				mtc.ids[[h]] <- out$mtc.ids
				sum.dist[[h]] <- out$sum.dist
				noad[[h]] <- out$noad
			}	
		}
		
		mmm <- unlist(lapply(mtc.ids, t))
		mmm <- substring(mmm, 5)
		mtc.ids <- matrix(mmm, ncol=2, byrow=TRUE)
		if(is.null(rownames(data.rec)) && is.null(rownames(data.don)))  mtc.ids <- matrix(as.numeric(mmm), ncol=2, byrow=TRUE)
		sss.dist <- NA
		if(!is.null(match.vars)){
			
			sss <- unlist(lapply(sum.dist, t))
			sss.dist <- matrix(sss, ncol=5, byrow=TRUE)
			colnames(sss.dist)<- colnames(out$sum.dist)
		}
		noad <- unlist(noad)
	}
	dimnames(mtc.ids) <- list(NULL, c("rec.id", "don.id"))
	list(mtc.ids=mtc.ids, sum.dist=sss.dist, noad=noad, call=match.call())
}

