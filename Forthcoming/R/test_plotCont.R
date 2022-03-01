data(samp.A, package = "StatMatch")

set.seed(23357)
pos <- sample(3009, 1500)

aa <- samp.A[pos,]
bb <- samp.A[-pos,]
library(ggplot2)
plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A=NULL, w.B=NULL,
                     type="density", ref=FALSE)
plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A="ww", w.B="ww",
         type="density", ref=FALSE)

plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A=NULL, w.B=NULL,
         type="density", ref=T)
plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A="ww", w.B="ww",
         type="density", ref=T)



plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A=NULL, w.B=NULL,
         type="hist", ref=FALSE)
plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A="ww", w.B="ww",
         type="hist", ref=FALSE)
plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A=NULL, w.B=NULL,
         type="hist", ref=T)
plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A="ww", w.B="ww",
         type="hist", ref=T)

plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A=NULL, w.B=NULL,
         type="ecdf", ref=FALSE)
plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A="ww", w.B="ww",
         type="ecdf", ref=FALSE)
plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A=NULL, w.B=NULL,
         type="ecdf", ref=T)
plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A="ww", w.B="ww",
         type="ecdf", ref=T)

plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A=NULL, w.B=NULL,
         type="qqplot", ref=FALSE)
plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A="ww", w.B="ww",
         type="qqplot", ref=FALSE)
plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A=NULL, w.B=NULL,
         type="qqplot", ref=T)
plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A="ww", w.B="ww",
         type="qqplot", ref=T)

plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A=NULL, w.B=NULL,
         type="qqshift", ref=FALSE)
plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A="ww", w.B="ww",
         type="qqshift", ref=FALSE)
plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A=NULL, w.B=NULL,
         type="qqshift", ref=T)
plotCont(data.A=aa, data.B=bb, xlab.A="n.income", xlab.B=NULL, w.A="ww", w.B="ww",
         type="qqshift", ref=T)


