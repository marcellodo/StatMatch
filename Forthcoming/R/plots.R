data(samp.A, package = "StatMatch")
data(samp.B, package = "StatMatch")



## ggplot
library(ggplot2)



compTab(data.A = samp.A, data.B = samp.B, xlab.A="edu7", w.A = "ww")
compTab(data.A = samp.A, data.B = samp.B, xlab.A=c("urb", "sex"), w.A = "ww", w.B="ww")

all(c("urb", "sex", "marital")==c("urb", "sex", "marital"))

compCont(data.A = samp.A, data.B = samp.B, xlab.A="age", w.A = "ww", w.B="ww")
compCont(data.A = samp.A, data.B = samp.B, xlab.A="age", type='qqplot')
compCont(data.A = samp.A, data.B = samp.B, xlab.A="age", w.A = "ww", type='qqplot')
compCont(data.A = samp.A, data.B = samp.B, xlab.A="age", type='qqshift')
compCont(data.A = samp.A, data.B = samp.B, xlab.A="age", type='hist')
compCont(data.A = samp.A, data.B = samp.B, xlab.A="age", type='hist',
         w.A="ww", w.B='ww')

############################
# frechet bounds
library(StatMatch)

# compute the tables required by Frechet.bounds.cat()
freq.xA <- xtabs(~sex+c.age, data=samp.A)
freq.xB <- xtabs(~sex+c.age, data=samp.B)
freq.xy <- xtabs(~sex+c.age+c.neti, data=samp.A)
freq.xz <- xtabs(~sex+c.age+labour5, data=samp.B)

# apply Frechet.bounds.cat()
library(StatMatch)
bounds.yz <- Frechet.bounds.cat(tab.x=freq.xA+freq.xB, tab.xy=freq.xy,
                                tab.xz=freq.xz, print.f="data.frame")

df <- bounds.yz[["bounds"]]
n <- nrow(df)
xs <- 1:n

plot(x = xs, y=df[,"CIA"], pch=18, 
     xlab = "", ylab="", cex=1,
     ylim = c(min(df[,"low.u"])-0.1, max(df[,"up.u"]))+0.1,
     xaxt="n")

for(i in 1:n){
  lines(x = c(i,i), y=c(df[i,"low.u"], df[i,"up.u"]), lwd=1, lty=3)
  lines(x = c(i,i), y=c(df[i,"low.cx"], df[i,"up.cx"]), lwd=2, lty=1)
}

lab1 <- paste(names(df)[1], df[,1], sep=" = ")
lab2 <- paste(names(df)[2], df[,2], sep=" = ")
lab <- paste(lab1, lab2, sep=", ")
lab
text(x=(1:n)-0.2, y=df[,"up.u"], labels = lab, pos =3 , cex=0.8, srt=90)


##########################################################
mdf <- data.matrix(df)
mdf
eps <- tapply(mdf[,"up.u"]-mdf[,"low.u"], mdf[,2], max)
aa <- round(max(eps),1) + 0.05
aa <- aa*(1:max(mdf[,2]))
xs <- rep(aa, each=max(mdf[,1]))
xxs <- xs + mdf[, -(1:2)]

yys <- max(mdf[,1])+1-mdf[,1]

par(mar=c(2.1, 7.1, 4.1, 1.1))

plot(x = xxs[,"CIA"], y=yys, pch=18,
     xlab="",ylab="", xaxt="n", yaxt="n",axes = F,
     xlim = c(min(xxs)-0.1, max(xxs)+0.1),
     ylim = c(min(yys)-1, max(yys)+0.5)
     
)
for(i in 1:n){
  lines(x = c(xxs[i, "low.u"], xxs[i, "up.u"]), y=c(yys[i], yys[i]), lty=3)
  lines(x = c(xxs[i, "low.cx"], xxs[i, "up.cx"]), y=c(yys[i], yys[i]), lty=1)
}

lab1 <- paste(names(df)[1], unique(df[,1]), sep=" = ")
lab2 <- paste(names(df)[2], unique(df[,2]), sep=" = ")

posCIA <- tapply(xxs[,"CIA"], df[,2], mean)
axis(side = 3, at = posCIA, labels = lab2)
axis(side = 2, at = unique(yys), labels = lab1, las=2)

wdtu <- df$up.u - df$low.u
wdtcx <- df$up.cx - df$low.cx

labwdt <- paste(round(wdtcx,3), " (", round(wdtu,3), ")", sep="")
labwdt
text(xxs[,"CIA"], y=yys, labels = round(df$CIA,3), pos = 3, cex=0.8) 
text(xxs[,"CIA"], y=yys, labels = labwdt, pos = 1, cex=0.8) 