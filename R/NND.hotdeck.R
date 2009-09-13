`NND.hotdeck` <-
function (data.rec, data.don, match.vars, don.class=NULL, dist.fun="Euclidean", constrained=FALSE, constr.alg=NULL) 
{

    if(dist.fun!="Gower" || dist.fun!="gower" || dist.fun!="exact" || dist.fun!="exact matching"){
		require(proxy)
	}
	if(constrained && constr.alg=="relax"){
		require(optmatch)
	}
	if(constrained && (constr.alg=="lpSolve" || constr.alg=="lpsolve")){
		require(lpSolve)
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
	
########################
NND.hd <- function (rec, don, dfun="Euclidean", constr=FALSE, c.alg=NULL)
{ 
    if(is.null(dim(rec))) x.rec <- data.frame(rec)
	else x.rec <- rec
    if(is.null(dim(don))) x.don <- data.frame(don)
	else x.don <- don
    p <- ncol(rec)
	nr <- nrow(x.rec)
	nd <- nrow(x.don)
	if(nr>nd) cat("Warning: the no. do donors is less than the no.f recipients", fill=TRUE)

    r.lab <- rownames(x.rec)
	if(is.null(r.lab)) r.lab <- 1:nr
	d.lab <- rownames(x.don)
	if(is.null(d.lab)) d.lab <- 1:nd

# compute matrix of distances between obs. in x.don and obs. in x.rec
# function dist() in package "proxy" is used! 

	if(dfun=="Euclidean" || dfun=="Manhattan"){
		cat("Warning:", dfun, "distance is being used", fill=TRUE)
        cat("Warning: all the categorical variables in rec and don data.frame are recoded into dummies", fill=TRUE)
        x.rec <- fact2dummy(x.rec, all=FALSE)
        x.don <- fact2dummy(x.don, all=FALSE)
        mdist <- dist(x=x.rec, y=x.don, method=dfun)
	}
	else if(dfun=="exact" || dfun=="exact matching"){
		cat("Warning: exact matching distance is being used", fill=TRUE)
		cat("all the matching variables in data.rec and data.don are converted to character variables and treated as categorical nominal", fill=TRUE)
		dxr <- dim(x.rec)
		x.rec <- as.character(as.matrix(x.rec))
		dim(x.rec) <- dxr
		dxd <- dim(x.don)
		x.don <- as.character(as.matrix(x.don))
		dim(x.don) <- dxd
		xx <- data.frame(rbind(x.rec, x.don))
		x.rec <- xx[1:nr,]
		x.don <- xx[-(1:nr),]
		mdist <- gower.dist(data.x=x.rec, data.y=x.don)
	}
	else if(dfun=="Gower" || dfun=="gower"){
		# if(p==1 && is.factor(x.rec)) x.rec <- list(x.rec)
		# if(p==1 && is.factor(x.don)) x.don <- list(x.don)
		mdist <- gower.dist(data.x=x.rec, data.y=x.don)
		mdist[is.nan(mdist)] <- 1 # NaN can occur when p=1 and x.rec and x.don is of type logical
		mdist[is.na(mdist)] <- 1 # NA can occur when p=1 and x.rec and x.don is of type logical
	}
	else{
		mdist <- dist(x=x.rec, y=x.don, method=dfun)
	}
	dimnames(mdist) <- list(r.lab, d.lab)

# UNCONSTRAINED nearest neighbour matching

	if(!constr){
		dist.rd <- numeric(nr)
		nad <- rep(NA, nr)
		don.lab <- numeric(nr)
		for(i in 1:nr){
			vd <- mdist[i,]
			min.d <- min(vd) # smallest distance recipient-donor
			dist.rd[i] <- min.d
			appo <- d.lab[vd==min.d]
			nad[i] <- length(appo) # number of availabe donors
			if(length(appo)==1) don.lab[i] <- appo 
			else don.lab[i] <- sample(appo, 1)
		}
		rec.lab <- r.lab
	}

# CONSTRAINED nearest neighbour matching.
# the functions in library lpSolve are used

	if(constr && (c.alg=="lpSolve" || c.alg=="lpsolve")){
		if(nr==nd) appo <- lp.assign(cost.mat=mdist)
		else if(nr<nd){
                r.sig <- rep("==", nr)
                r.rhs <- rep(1, nr)
                c.sig <- rep("<=", nd)
                c.rhs <- rep(1, nd)
                appo <- lp.transport(cost.mat=mdist, row.signs=r.sig, row.rhs=r.rhs, col.signs=c.sig, col.rhs=c.rhs)
		}   
		else if(nr > nd){
			stop("There more recipients than donors!")
			cat("some donors will be used more than once", fill=TRUE)
			r.sig <- rep("==", nr)
			r.rhs <- rep(1, nr)
			c.sig <- rep(">=", nd)
			c.rhs <- rep(1, nd)
			appo <- lp.transport(cost.mat=mdist, row.signs=r.sig, row.rhs=r.rhs, col.signs=c.sig, col.rhs=c.rhs)
		}
		sol <- appo$solution
		ss <- c(t(sol))
		cc <- c(t(col(sol)))
		dist.rd <- mdist[cbind(1:nr, cc[as.logical(ss)] )]
		rec.lab <- r.lab
		don.lab <- d.lab[c(cc[as.logical(ss)])]
	}

# the function pairmatch() in library optMatch are used
	if(constr && c.alg=="relax"){
		if(nr > nd) stop("pairmatch() function in library optmatch requires the no. 
							of donors to be greater or equal than the no. of recipients")
		out.pr <- pairmatch(mdist)
		labs <- names(out.pr)
		tt <- labs %in% r.lab
		df1 <- data.frame(id.rec=labs[tt], mm=out.pr[tt], stringsAsFactors=FALSE)
		df2 <- data.frame(id.don=labs[!tt], mm=out.pr[!tt], stringsAsFactors=FALSE)
		df <- merge(df1, df2, by="mm")
		rmd <- mdist[df$id.rec, ]
		rmd <- rmd[ ,df$id.don]
		dist.rd <- diag(rmd[1:nr, 1:nr])
		aa <- data.frame(id.rec=r.lab, pos=1:nr, stringsAsFactors=FALSE)
		bb <- merge(aa, df, by="id.rec")
		don.lab <- bb[order(bb$pos) , "id.don"]
		rec.lab <- r.lab
		dist.rd <- dist.rd[order(bb$pos)]
	}
# output
    mtc.ids <- cbind(rec.id=rec.lab, don.id=don.lab)
	if(constr) fine <- list(mtc.ids=mtc.ids, dist.rd=dist.rd, call=match.call())
	else fine <- list(mtc.ids=mtc.ids, dist.rd=dist.rd, noad=nad, call=match.call())
	fine
}
################ NND.hd ends here #############################
	
	if(is.null(don.class)){ 
		out <- NND.hd(rec=data.rec[,match.vars], don=data.don[,match.vars], dfun=dist.fun, constr=constrained, c.alg=constr.alg )
		mmm <- out$mtc.ids
		mmm <- substring(mmm, 5)
		if(is.null(rownames(data.rec)) && is.null(rownames(data.don)))  mtc.ids <- matrix(as.numeric(mmm), ncol=2, byrow=TRUE)
		else mtc.ids <- mmm
		dimnames(mtc.ids) <- list(NULL, c("rec.id", "don.id"))
		dist.rd <- out$dist.rd
		if(!constrained) noad <- out$noad
	}
	else{
		if(length(don.class)==1){
			l.rec <- split(data.rec[ ,match.vars], f=data.rec[ ,don.class])
			l.don <- split(data.don[ ,match.vars], f=data.don[ ,don.class])
		}
		else{
			l.rec <- split(data.rec[ ,match.vars], f=as.list(data.rec[ ,don.class]))
			l.don <- split(data.don[ ,match.vars], f=as.list(data.don[ ,don.class]))
		}
		if(length(l.rec)!=length(l.don)){
			cat("The no. of donation classes in recipient data is not equal to the no. of donation classes in donor data", fill=TRUE)
			stop("Possible reason: the variables used to classify units are not defined as factors or are factors with different levels") 	
		}
		if(!identical(names(l.rec), names(l.don)))
			cat("Warning: the donation classes seem built using different factors with differnt levels")
		if(p==1){
            nn.r <- unlist(lapply(l.rec, length))
            nn.d <- unlist(lapply(l.don, length))
        }
		else {
            nn.r <- unlist(lapply(l.rec, nrow))
            nn.d <- unlist(lapply(l.don, nrow))
        }
        l.rec <- l.rec[nn.r>0]
        l.don <- l.don[nn.r>0]
        nn.r <- nn.r[nn.r>0]
        nn.d <- nn.d[nn.r>0]
        if(any(nn.d==0)) {
			stop("For some donation classes there are NO donors available. Please modify the definition of the donation classes")
		}	
		H <- length(l.rec)
		mtc.ids <- as.list(numeric(H))
		dist.rd <- as.list(numeric(H))
		if(!constrained) noad <- as.list(numeric(H))
		for(h in 1:H){
        out <- NND.hd(rec=l.rec[[h]], don=l.don[[h]], dfun=dist.fun, constr=constrained, c.alg=constr.alg )
        mtc.ids[[h]] <- out$mtc.ids
			  dist.rd[[h]] <- out$dist.rd
			  if(!constrained) noad[[h]] <- out$noad
		}
        mmm <- unlist(lapply(mtc.ids, t))
		mmm <- substring(mmm, 5)
		mtc.ids <- matrix(mmm, ncol=2, byrow=TRUE)
		if(is.null(rownames(data.rec)) && is.null(rownames(data.don)))  mtc.ids <- matrix(as.numeric(mmm), ncol=2, byrow=TRUE)
		dimnames(mtc.ids) <- list(NULL, c("rec.id", "don.id"))
		dist.rd <- unlist(dist.rd)
		if(!constrained) noad <- unlist(noad)
	}
	if(constrained) end.out <- list(mtc.ids=mtc.ids, dist.rd=dist.rd, call=match.call())
	else end.out <- list(mtc.ids=mtc.ids, dist.rd=dist.rd, noad=noad, call=match.call())
	end.out
}

