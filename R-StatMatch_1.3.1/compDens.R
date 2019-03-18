compDens <- function(data.A, data.B, xlab, wA=NULL, wB=wA){
    
    if(is.null(wA)){
        xA = data.frame(x=data.A[, xlab], sample="A")
        xB = data.frame(x=data.B[, xlab], sample="B")
        xx <- rbind(xA, xB)
        colnames(xx) <- c(xlab, "sample")
        
        ggplot(xx, aes(xx[,xlab], fill=sample, colour=sample)) +
            geom_density(alpha=0.4, lwd=0.8, adjust=0.5) +
            labs(x = xlab)
    }
    else{
        wwA <- data.A[, wA]/sum(data.A[, wA])
        xA = data.frame(x=data.A[, xlab], w=wwA, sample="A")
        
        wwB <- data.B[, wB]/sum(data.B[, wB])
        xB = data.frame(x=data.B[, xlab], w=wwB, sample="B")
        xx <- rbind(xA, xB)
        colnames(xx) <- c(xlab, "w", "sample")
        
        ggplot(xx, aes(xx[,xlab], weight=w, fill=sample, colour=sample)) +
            geom_density(alpha=0.4, lwd=0.8, adjust=0.5) +
            labs(x = xlab)
    }
}