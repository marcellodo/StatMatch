# compute the tables required by Frechet.bounds.cat()
freq.xy <- xtabs(~c.neti, data=samp.A)
freq.xz <- xtabs(~labour5, data=samp.B)

# apply Frechet.bounds.cat()
library(StatMatch)
bounds.yz <- Frechet.bounds.cat(tab.x=NULL, tab.xy=freq.xy,
                                tab.xz=freq.xz, print.f="data.frame")

plot.bounds(bounds.yz)

bounds.yz <- Frechet.bounds.cat(tab.x=NULL, tab.xy=freq.xy,
                                tab.xz=freq.xz)

plot.bounds(bounds.yz)


###########################################################à
# # compute the tables required by Frechet.bounds.cat()
# freq.xA <- xtabs(~sex+c.age, data=samp.A)
# freq.xB <- xtabs(~sex+c.age, data=samp.B)
# freq.xy <- xtabs(~sex+c.age+c.neti, data=samp.A)
# freq.xz <- xtabs(~sex+c.age+labour5, data=samp.B)
# 
# # apply Frechet.bounds.cat()
# bounds.yz <- Frechet.bounds.cat(tab.x=freq.xA+freq.xB, tab.xy=freq.xy,
#                                 tab.xz=freq.xz, print.f="data.frame")
# 
# 
# plot.bounds(bounds.yz)

bounds.yz <- Frechet.bounds.cat(tab.x=freq.xA+freq.xB, tab.xy=freq.xy,
                                tab.xz=freq.xz)


plot.bounds(bounds.yz)


##########################################à
freq.xA <- xtabs(~sex+c.age+edu7, data=samp.A)
freq.xB <- xtabs(~sex+c.age+edu7, data=samp.B)
freq.xy <- xtabs(~sex+c.age+edu7+c.neti, data=samp.A)
freq.xz <- xtabs(~sex+c.age+edu7+labour5, data=samp.B)

# apply Frechet.bounds.cat()
bounds.yz <- Frechet.bounds.cat(tab.x=freq.xA+freq.xB, tab.xy=freq.xy,
                                tab.xz=freq.xz, print.f="data.frame")


plot.bounds(bounds.yz)

bounds.yz <- Frechet.bounds.cat(tab.x=freq.xA+freq.xB, tab.xy=freq.xy,
                                tab.xz=freq.xz)


plot.bounds(bounds.yz)


