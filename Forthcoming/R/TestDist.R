library(StatMatch)
data("samp.A")
data("samp.B")

lab <- c("hsize5", "age", "sex")

gower.dist(data.x = samp.A[1:10,lab], data.y = samp.B[1:10,lab], cat.cont = F)

gower.dist(data.x = samp.A[1:10,lab], data.y = samp.B[1:10,lab], cat.cont = T, kern = "kde0")
gower.dist(data.x = samp.A[1:10,lab], data.y = samp.B[1:10,lab], cat.cont = T, kern = "kde1")
gower.dist(data.x = samp.A[1:10,lab], data.y = samp.B[1:10,lab], cat.cont = T, kern = "kde2")


gower.dist(data.x = samp.A[1:10,"age"], data.y = samp.B[1:10,"age"], cat.cont = F)
gower.dist(data.x = samp.A[1:10,"age"], data.y = samp.B[1:10,"age"], cat.cont = T)
