`fact2dummy` <-
function (data, all=TRUE) 
{
	dum.fcn <- function(x, all=TRUE){
		n <- length(x)
		k <- nlevels(x)
		lev <- levels(x)
		mat <- matrix(lev, n, k, byrow=TRUE)
		fine <- (mat==x)*1
		dim(fine) <- c(n,k)
		colnames(fine) <- lev
		if(!all) {
			fine <- fine[,-k]
			fine <- matrix(fine, n, (k-1))
			colnames(fine) <- lev[-k]
		}
		fine
	} 
	
	if(is.null(dim(data))){
		if(class(data)[1]=="numeric" || class(data)[1]=="integer" || class(data)[1]=="logical") oo <- cbind(1*data)
		else{
			oo <- dum.fcn(data, all=all)
			colnames(oo) <- paste("x", colnames(oo), sep="=")
		}
	}
	else{
		p <- length(data)
		n <- nrow(data)
		vtype <- lapply(data, class)
		out <- as.list(rep(NA,p))
		
		for(i in 1:p){
			if(vtype[[i]][1]=="numeric" || vtype[[i]][1]=="integer" || vtype[[i]][1]=="logical") {
				aa <- matrix(1*data[,i])
				colnames(aa) <- names(data)[i]
				out[[i]] <- aa
			}	
			else {
				aa <- dum.fcn(data[,i], all=all)
				colnames(aa) <- paste(names(data)[i], colnames(aa), sep="=")
				out[[i]] <- aa
			}
		}
		oo <- unlist(out)
		oo <- matrix(oo, nrow=n)
		dimnames(oo) <- list(row.names(data), unlist(lapply(out, colnames)))
	}	
	oo
}

