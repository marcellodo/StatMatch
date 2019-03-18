data(samp.A, package = "StatMatch")
data(samp.B, package = "StatMatch")


lab <- "age"
# density compare
xA <- samp.A[, lab]
xB <- samp.B[, lab]

dta_A <- density(xA, na.rm = TRUE)
dta_B <- density(xB, na.rm = TRUE)

plot(dta_A, col = "blue", main = paste0("Dens ", lab), 
     ylim = c(0, max(dta_A$y,dta_B$y)))
lines(dta_B, col = "red")
legend(, paste0(lab, c(" A"," B")),
       lty = c(1,1), col = c("blue","red"))


xrng = range(xA, xB)

kdeA = density(xA, from = xrng[1L], to = xrng[2L])
kdeB = density(xB, from = xrng[1L], to = xrng[2L])

matplot(kdeA$x, cbind(kdeA$y, kdeB$y))

## ggplot
library(ggplot2)




data.A <- samp.A
data.B <- samp.B
xlab <- "urb"
compTab(data.A = samp.A, data.B = samp.B, xlab="urb")
compTab(data.A = samp.A, data.B = samp.B, xlab="urb", wA = "ww")

compDens(data.A = samp.A, data.B = samp.B, xlab="age")
compDens(data.A = samp.A, data.B = samp.B, xlab="age", wA = "ww")
    
    