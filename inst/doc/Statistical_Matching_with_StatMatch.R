### R code from vignette source 'Statistical_Matching_with_StatMatch.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: Statistical_Matching_with_StatMatch.Rnw:108-112
###################################################
options(useFancyQuotes="UTF-8")
#options(useFancyQuotes=FALSE)
options(width=66)
options(warn=-1)


###################################################
### code chunk number 2: Statistical_Matching_with_StatMatch.Rnw:116-119
###################################################
library(simPopulation, warn.conflicts=FALSE) #loads pkg simPopulation
data(eusilcS) 
str(eusilcS)


###################################################
### code chunk number 3: Statistical_Matching_with_StatMatch.Rnw:124-144
###################################################
# discard units with age<16
silc.16 <- subset(eusilcS, age>15) # units 
nrow(silc.16)
# categorize age
silc.16$c.age <- cut(silc.16$age, c(16,24,49,64,100), include.lowest=T)
#
# truncate hsize
aa <- as.numeric(silc.16$hsize)
aa[aa>6] <- 6
silc.16$hsize6 <- factor(aa, ordered=T)
#
# recode personal economic status
aa <- as.numeric(silc.16$pl030)
aa[aa<3] <- 1
aa[aa>1] <- 2
silc.16$work <- factor(aa, levels=1:2, labels=c("working","not working"))
#
# categorize personal net income
silc.16$c.netI <- cut(silc.16$net/1000,
                      breaks=c(-6,0,5,10,15,20,25,30,40,50,200))


###################################################
### code chunk number 4: Statistical_Matching_with_StatMatch.Rnw:149-169
###################################################
# simulate samples
set.seed(123456)
obs.A <- sample(nrow(silc.16), 4000, replace=F)

X.vars <- c("hsize","hsize6","db040","age","c.age",
            "rb090","pb220a","rb050")
y.var <- c("pl030","work")
z.var <- c("netIncome", "c.netI")

rec.A <- silc.16[obs.A, c(X.vars, y.var)]
don.B <- silc.16[-obs.A, c(X.vars, z.var)]

#
# determine a rough weighting 
# compute N, the est. size of pop(age>16)
N <- round(sum(silc.16$rb050)) 
N
#rescale origin weights
rec.A$wwA <- rec.A$rb050/sum(rec.A$rb050)*N
don.B$wwB <- don.B$rb050/sum(don.B$rb050)*N


###################################################
### code chunk number 5: Statistical_Matching_with_StatMatch.Rnw:183-189
###################################################
# analyses on A
library(StatMatch) #loads StatMatch
# response is pl030
pw.assoc(pl030~db040+hsize6+c.age+rb090+pb220a, data=rec.A)
#response is work (aggregated pl030)
pw.assoc(work~db040+hsize6+c.age+rb090+pb220a, data=rec.A)


###################################################
### code chunk number 6: Statistical_Matching_with_StatMatch.Rnw:196-199
###################################################
# analyses on B
require(Hmisc)
spearman2(netIncome~db040+hsize+age+rb090+pb220a, data=don.B)


###################################################
### code chunk number 7: Statistical_Matching_with_StatMatch.Rnw:237-246
###################################################
xx <- xtabs(~db040+hsize6+c.age+rb090+pb220a, data=rec.A)
xy <- xtabs(~db040+hsize6+c.age+rb090+pb220a+work, data=rec.A)
xz <- xtabs(~db040+hsize6+c.age+rb090+pb220a+c.netI, data=don.B)

library(StatMatch) #loads StatMatch
out.fbw <-  Fbwidths.by.x(tab.x=xx, tab.xy=xy, tab.xz=xz)
# sort according to overall uncertainty
sort.ov.unc <- out.fbw$sum.unc[order(out.fbw$sum.unc$ov.unc),]
head(sort.ov.unc) # best 6 models


###################################################
### code chunk number 8: Statistical_Matching_with_StatMatch.Rnw:267-271
###################################################
group.v <- c("rb090","db040")
X.mtc <- "age" 
out.nnd <- NND.hotdeck(data.rec=rec.A, data.don=don.B,
                         match.vars=X.mtc, don.class=group.v)


###################################################
### code chunk number 9: Statistical_Matching_with_StatMatch.Rnw:276-279
###################################################
summary(out.nnd$dist.rd) # summary distances rec-don
summary(out.nnd$noad) # summary available donors at min. dist.
table(out.nnd$noad)


###################################################
### code chunk number 10: Statistical_Matching_with_StatMatch.Rnw:284-289
###################################################
head(out.nnd$mtc.ids)
fA.nnd <- create.fused(data.rec=rec.A, data.don=don.B,
                       mtc.ids=out.nnd$mtc.ids,
                       z.vars=c("netIncome","c.netI"))
head(fA.nnd) #first 6 obs.


###################################################
### code chunk number 11: Statistical_Matching_with_StatMatch.Rnw:296-305
###################################################
group.v <- c("rb090","db040")
X.mtc <- "age"
out.nnd.c <- NND.hotdeck(data.rec=rec.A, data.don=don.B, 
                           match.vars=X.mtc, don.class=group.v, 
                           dist.fun="Manhattan", constrained=TRUE, 
                           constr.alg="Hungarian")
fA.nnd.c <- create.fused(data.rec=rec.A, data.don=don.B,
                    mtc.ids=out.nnd.c$mtc.ids,
                    z.vars=c("netIncome","c.netI"))


###################################################
### code chunk number 12: Statistical_Matching_with_StatMatch.Rnw:310-313
###################################################
#comparing distances
sum(out.nnd$dist.rd) # unconstrained
sum(out.nnd.c$dist.rd) # constrained


###################################################
### code chunk number 13: Statistical_Matching_with_StatMatch.Rnw:318-326
###################################################
# estimating marginal distribution of C.netI
tt0 <- xtabs(~c.netI, data=don.B) # reference distr.
tt <- xtabs(~c.netI, data=fA.nnd)  # synt unconstr.
ttc <- xtabs(~c.netI, data=fA.nnd.c) #synt. constr.
#
# comparing marginal distributions
comp.prop(p1=tt, p2=tt0, n1=nrow(fA.nnd), n2=NULL, ref=TRUE)
comp.prop(p1=ttc, p2=tt0, n1=nrow(fA.nnd), n2=NULL, ref=TRUE)


###################################################
### code chunk number 14: Statistical_Matching_with_StatMatch.Rnw:336-342
###################################################
group.v <- c("db040","rb090")
rnd.1 <- RANDwNND.hotdeck(data.rec=rec.A, data.don=don.B, 
                           match.vars=NULL, don.class=group.v)
fA.rnd <- create.fused(data.rec=rec.A, data.don=don.B,
                        mtc.ids=rnd.1$mtc.ids, 
                        z.vars=c("netIncome", "c.netI"))


###################################################
### code chunk number 15: Statistical_Matching_with_StatMatch.Rnw:349-359
###################################################
# random choiches of a donor among the closest k=20 wrt age
group.v <- c("db040","rb090")
X.mtc <- "age"
rnd.2 <- RANDwNND.hotdeck(data.rec=rec.A, data.don=don.B, 
                           match.vars=X.mtc, don.class=group.v, 
                           dist.fun="Manhattan", 
                            cut.don="exact", k=20)
fA.knnd <- create.fused(data.rec=rec.A, data.don=don.B,
                        mtc.ids=rnd.2$mtc.ids, 
                        z.vars=c("netIncome", "c.netI"))


###################################################
### code chunk number 16: Statistical_Matching_with_StatMatch.Rnw:364-365
###################################################
head(rnd.2$sum.dist)


###################################################
### code chunk number 17: Statistical_Matching_with_StatMatch.Rnw:383-391
###################################################
rnk.1 <- rankNND.hotdeck(data.rec=rec.A, data.don=don.B, 
                          var.rec="age", var.don="age")
#create the synthetic data set
fA.rnk <- create.fused(data.rec=rec.A, data.don=don.B,
                        mtc.ids=rnk.1$mtc.ids, 
                        z.vars=c("netIncome", "c.netI"), 
                        dup.x=TRUE, match.vars="age")
head(fA.rnk)


###################################################
### code chunk number 18: Statistical_Matching_with_StatMatch.Rnw:396-404
###################################################
rnk.2 <- rankNND.hotdeck(data.rec=rec.A, data.don=don.B, var.rec="age",
                        var.don="age", don.class="rb090",
                         constrained=TRUE, constr.alg="Hungarian")
fA.grnk <- create.fused(data.rec=rec.A, data.don=don.B,
                        mtc.ids=rnk.2$mtc.ids, 
                        z.vars=c("netIncome", "c.netI"),
                        dup.x=TRUE, match.vars="age")
head(fA.grnk)


###################################################
### code chunk number 19: Statistical_Matching_with_StatMatch.Rnw:429-456
###################################################
# step 0) introduce missing values in iris
set.seed(1324)
miss <- rbinom(150, 1, 0.30) #generates randomly missing
data(iris, package="datasets")
iris.miss <- iris
iris.miss$Petal.Length[miss==1] <- NA
summary(iris.miss$Petal.L)
#
# step 1) separate units in two data sets
rec <- subset(iris.miss, is.na(Petal.Length), select=-Petal.Length)
don <- subset(iris.miss, !is.na(Petal.Length))
#
# step 2) search for closest donors
X.mtc <- c("Sepal.Length", "Sepal.Width", "Petal.Width")
nnd <- NND.hotdeck(data.rec=rec, data.don=don,
                         match.vars=X.mtc, don.class="Species",
                         dist.fun="Manhattan")
# fills rec
imp.rec <- create.fused(data.rec=rec, data.don=don,
                        mtc.ids=nnd$mtc.ids, z.vars="Petal.Length")
imp.rec$imp.PL <- 1 # flag for imputed
#
# step 3) re-aggregate data sets
don$imp.PL <- 0
imp.iris <- rbind(imp.rec, don)
#summary stat of imputed and non imputed Petal.Length
tapply(imp.iris$Petal.Length, imp.iris$imp.PL, summary)


###################################################
### code chunk number 20: Statistical_Matching_with_StatMatch.Rnw:472-489
###################################################
# uses iris data set
iris.A <- iris[101:150, 1:3]
iris.B <- iris[1:100, c(1:2,4)]

X.mtc <- c("Sepal.Length","Sepal.Width") # matching variables

# parameters estimated using ML
mix.1 <- mixed.mtc(data.rec=iris.A, data.don=iris.B, match.vars=X.mtc,
                    y.rec="Petal.Length", z.don="Petal.Width", 
                    method="ML", rho.yz=0, 
                    micro=TRUE, constr.alg="Hungarian")

mix.1$mu #estimated means
mix.1$cor #estimated cor. matrix

head(mix.1$filled.rec) # A filled in with Z
cor(mix.1$filled.rec)


###################################################
### code chunk number 21: Statistical_Matching_with_StatMatch.Rnw:496-503
###################################################
# parameters estimated using ML and rho_YZ|X=0.85
mix.2 <- mixed.mtc(data.rec=iris.A, data.don=iris.B, match.vars=X.mtc,
                    y.rec="Petal.Length", z.don="Petal.Width", 
                    method="ML", rho.yz=0.85, 
                    micro=TRUE, constr.alg="Hungarian")
mix.2$cor
head(mix.2$filled.rec)


###################################################
### code chunk number 22: Statistical_Matching_with_StatMatch.Rnw:508-514
###################################################
mix.3 <- mixed.mtc(data.rec=iris.A, data.don=iris.B, match.vars=X.mtc,
                    y.rec="Petal.Length", z.don="Petal.Width", 
                    method="MS", rho.yz=0.75, 
                    micro=TRUE, constr.alg="Hungarian")

mix.3$rho.yz


###################################################
### code chunk number 23: Statistical_Matching_with_StatMatch.Rnw:529-554
###################################################
# summary info on the weights
sum(rec.A$wwA) # estimated pop size from A
sum(don.B$wwB) # estimated pop size from B
summary(rec.A$wwA)
summary(don.B$wwB)

# NND constrained hot deck
group.v <- c("rb090","db040")
out.nnd <- NND.hotdeck(data.rec=rec.A, data.don=don.B,
                         match.vars="age", don.class=group.v,
                         dist.fun="Manhattan",
                         constrained=TRUE, constr.alg="Hungarian")

fA.nnd.m <- create.fused(data.rec=rec.A, data.don=don.B,
                    mtc.ids=out.nnd$mtc.ids,
                    z.vars=c("netIncome","c.netI"))

# estimating average net income
weighted.mean(fA.nnd.m$netIncome, fA.nnd.m$wwA) # imputed in A
weighted.mean(don.B$netIncome, don.B$wwB) # ref. estimate in B

# comparing marginal distribution of C.netI using weights
tt.0w <- xtabs(wwB~c.netI, data=don.B)
tt.fw <- xtabs(wwA~c.netI, data=fA.nnd.m)
comp.prop(p1=tt.fw, p2=tt.0w, n1=nrow(fA.nnd.m), ref=TRUE)


###################################################
### code chunk number 24: Statistical_Matching_with_StatMatch.Rnw:561-577
###################################################
group.v <- c("rb090","db040")
X.mtc <- "age"
rnd.2 <- RANDwNND.hotdeck(data.rec=rec.A, data.don=don.B, 
                           match.vars=NULL, don.class=group.v, 
                           weight.don="wwB")
fA.wrnd <- create.fused(data.rec=rec.A, data.don=don.B, 
                         mtc.ids=rnd.2$mtc.ids,
                         z.vars=c("netIncome","c.netI"))

weighted.mean(fA.wrnd$netIncome, fA.wrnd$wwA) # imputed in A
weighted.mean(don.B$netIncome, don.B$wwB) # ref. estimate in B

# comparing marginal distribution of C.netI using weights
tt.0w <- xtabs(wwB~c.netI, data=don.B)
tt.fw <- xtabs(wwA~c.netI, data=fA.wrnd)
comp.prop(p1=tt.fw, p2=tt.0w, n1=nrow(fA.nnd.m), ref=TRUE)


###################################################
### code chunk number 25: Statistical_Matching_with_StatMatch.Rnw:588-607
###################################################
rnk.w <- rankNND.hotdeck(data.rec=rec.A, data.don=don.B, 
                          don.class="db040", var.rec="age", 
                          var.don="age", weight.rec="wwA",
                          weight.don="wwB", constrained=TRUE,
                          constr.alg="Hungarian")
#
#create the synthetic data set
fA.wrnk <- create.fused(data.rec=rec.A, data.don=don.B,
                        mtc.ids=rnk.w$mtc.ids, 
                        z.vars=c("netIncome", "c.netI"), 
                        dup.x=TRUE, match.vars="age")
#
weighted.mean(fA.wrnk$netIncome, fA.wrnk$wwA) # imputed in A
weighted.mean(don.B$netIncome, don.B$wwB) # ref. estimate in B

# comparing marginal distribution of C.netI using weights
tt.0w <- xtabs(wwB~c.netI, data=don.B)
tt.fw <- xtabs(wwA~c.netI, data=fA.wrnk)
comp.prop(p1=tt.fw, p2=tt.0w, n1=nrow(fA.nnd.m), ref=TRUE)


###################################################
### code chunk number 26: Statistical_Matching_with_StatMatch.Rnw:630-652
###################################################
tt.A <- xtabs(wwA~rb090+c.age, data=rec.A)
tt.B <- xtabs(wwB~rb090+c.age, data=don.B)
(prop.table(tt.A)-prop.table(tt.B))*100
comp.prop(p1=tt.A, p2=tt.B, n1=nrow(rec.A),
          n2=nrow(don.B), ref=FALSE)

library(survey, warn.conflicts=FALSE) # loads survey
# creates svydesign objects
svy.rec.A <- svydesign(~1, weights=~wwA, data=rec.A)
svy.don.B <- svydesign(~1, weights=~wwB, data=don.B)
#
# harmonizes wrt to joint distr. of gender vs. c.age
out.hz <- harmonize.x(svy.A=svy.rec.A, svy.B=svy.don.B,
                form.x=~c.age:rb090-1)
#
summary(out.hz$weights.A) # new calibrated weights for A
summary(out.hz$weights.B) # new calibrated weights for B

tt.A <- xtabs(out.hz$weights.A~rb090+c.age, data=rec.A)
tt.B <- xtabs(out.hz$weights.B~rb090+c.age, data=don.B)
comp.prop(p1=tt.A, p2=tt.B, n1=nrow(rec.A),
          n2=nrow(don.B), ref=FALSE)


###################################################
### code chunk number 27: Statistical_Matching_with_StatMatch.Rnw:668-674
###################################################
# estimating c.pl030 vs. c.netI under the CI assumption
out <- comb.samples(svy.A=out.hz$cal.A, svy.B=out.hz$cal.B,
            svy.C=NULL, y.lab="work", z.lab="c.netI",
            form.x=~c.age:rb090-1)
#
addmargins(t(out$yz.CIA))  # table estimated under the CIA


###################################################
### code chunk number 28: Statistical_Matching_with_StatMatch.Rnw:679-697
###################################################
# generating artificial sample C
set.seed(43210)
obs.C <- sample(nrow(silc.16), 200, replace=F)
#
X.vars <- c("hsize","hsize6","db040","age","c.age",
            "rb090","pb220a", "rb050")
y.var <- c("pl030","work")
z.var <- c("netIncome","c.netI")
#
aux.C <- silc.16[obs.C, c(X.vars, y.var, z.var)]
aux.C$wwC <- aux.C$rb050/sum(aux.C$rb050)*round(sum(silc.16$rb050)) # rough w
svy.aux.C <- svydesign(~1, weights=~wwC, data=aux.C) 
#
# incomplete two-way estimation
out.inc <- comb.samples(svy.A=out.hz$cal.A, svy.B=out.hz$cal.B,
                svy.C=svy.aux.C, y.lab="work", z.lab="c.netI",
                form.x=~c.age:rb090-1, estimation="incomplete")
addmargins(t(out.inc$yz.est))           


###################################################
### code chunk number 29: Statistical_Matching_with_StatMatch.Rnw:702-722
###################################################
new.wwC <- weights(out.inc$cal.C) #new cal. weights for C 
#
# marginal distributions of work
m.work.cA <- xtabs(out.hz$weights.A~work, data=rec.A)
m.work.cC <- xtabs(new.wwC~work, data=aux.C)
m.work.cA-m.work.cC
#
# marginal distributions of c.netI
m.cnetI.cB <- xtabs(out.hz$weights.B~c.netI, data=don.B)
m.cnetI.cC <- xtabs(new.wwC~c.netI, data=aux.C)
m.cnetI.cB-m.cnetI.cC 

# joint distribution of the matching variables
tt.A <- xtabs(out.hz$weights.A~rb090+c.age, data=rec.A)
tt.B <- xtabs(out.hz$weights.B~rb090+c.age, data=don.B)
tt.C <- xtabs(new.wwC~rb090+c.age, data=aux.C)
comp.prop(p1=tt.A, p2=tt.B, n1=nrow(rec.A),
          n2=nrow(don.B), ref=FALSE)
comp.prop(p1=tt.C, p2=tt.A, n1=nrow(aux.C),
          n2=nrow(rec.A), ref=FALSE)


###################################################
### code chunk number 30: Statistical_Matching_with_StatMatch.Rnw:727-733
###################################################
# synthetic two-way estimation
out.synt <- comb.samples(svy.A=out.hz$cal.A, svy.B=out.hz$cal.B,
                svy.C=svy.aux.C, y.lab="work", z.lab="c.netI",
                form.x=~c.age:rb090-1, estimation="synthetic")
#
addmargins(t(out.synt$yz.est))           


###################################################
### code chunk number 31: Statistical_Matching_with_StatMatch.Rnw:738-758
###################################################
new.wwC <- weights(out.synt$cal.C) #new cal. weights for C 
#
# marginal distributions of work
m.work.cA <- xtabs(out.hz$weights.A~work, data=rec.A)
m.work.cC <- xtabs(new.wwC~work, data=aux.C)
m.work.cA-m.work.cC

# marginal distributions of c.netI
m.cnetI.cB <- xtabs(out.hz$weights.B~c.netI, data=don.B)
m.cnetI.cC <- xtabs(new.wwC~c.netI, data=aux.C)
m.cnetI.cB-m.cnetI.cC 

# joint distribution of the matching variables
tt.A <- xtabs(out.hz$weights.A~rb090+c.age, data=rec.A)
tt.B <- xtabs(out.hz$weights.B~rb090+c.age, data=don.B)
tt.C <- xtabs(new.wwC~rb090+c.age, data=aux.C)
comp.prop(p1=tt.A, p2=tt.B, n1=nrow(rec.A),
          n2=nrow(don.B), ref=FALSE)
comp.prop(p1=tt.C, p2=tt.A, n1=nrow(aux.C),
          n2=nrow(rec.A), ref=FALSE)


###################################################
### code chunk number 32: Statistical_Matching_with_StatMatch.Rnw:766-801
###################################################
# predicting prob of c.netI in A under the CI assumption
out <- comb.samples(svy.A=out.hz$cal.A, svy.B=out.hz$cal.B,
                    svy.C=NULL, y.lab="work", z.lab="c.netI",
                    form.x=~c.age:rb090-1, micro=TRUE)
head(out$Z.A)

# predicting prob of c.netI in A under the CI assumption
out <- comb.samples(svy.A=out.hz$cal.A, svy.B=out.hz$cal.B,
                    svy.C=NULL, y.lab="work", z.lab="c.netI",
                    form.x=~c.age:rb090-1, micro=TRUE)
head(out$Z.A)
sum(out$Z.A<0) # negative est. prob.
sum(out$Z.A>1) # est. prob. >1

# compare marginal distributions of Z
t.zA <- colSums(out$Z.A*out.hz$weights.A)
t.zB <- xtabs(out.hz$weights.B~don.B$c.netI)
comp.prop(p1=t.zA, p2=t.zB, n1=nrow(rec.A), ref=TRUE)  

# predicting class of netIncome in A
# randomized prediction with prob proportional to estimated prob.
pred.zA <- apply(out$Z.A,1,sample,x=1:ncol(out$Z.A), size=1,replace=F)
rec.A$c.netI <- factor(pred.zA, levels=1:nlevels(don.B$c.netI), 
                   labels=as.character(levels(don.B$c.netI)), ordered=T)

# comparing marginal distributions of Z
t.zA <- xtabs(out.hz$weights.A~rec.A$c.netI)
comp.prop(p1=t.zA, p2=t.zB, n1=nrow(rec.A), ref=TRUE)  

# comparing joint distributions of X vs. Z
t.xzA <- xtabs(out.hz$weights.A~c.age+rb090+c.netI, data=rec.A)
t.xzB <- xtabs(out.hz$weights.B~c.age+rb090+c.netI, data=don.B)
out.comp <- comp.prop(p1=t.xzA, p2=t.xzB, n1=nrow(rec.A), ref=TRUE)  
out.comp$meas
out.comp$chi.sq


###################################################
### code chunk number 33: Statistical_Matching_with_StatMatch.Rnw:842-853
###################################################
#comparing joint distribution of the X_M variables in A and in B
t.xA <- xtabs(wwA~c.age+rb090, data=rec.A)
t.xB <- xtabs(wwB~c.age+rb090, data=don.B)
comp.prop(p1=t.xA, p2=t.xB, n1=nrow(rec.A), n2=nrow(don.B), ref=FALSE)
#
#computing tables needed by Frechet.bounds.cat
t.xy <- xtabs(wwA~c.age+rb090+work, data=rec.A)
t.xz <- xtabs(wwB~c.age+rb090+c.netI, data=don.B)
out.fb <- Frechet.bounds.cat(tab.x=t.xA, tab.xy=t.xy, tab.xz=t.xz, 
                             print.f="data.frame")
out.fb


###################################################
### code chunk number 34: Statistical_Matching_with_StatMatch.Rnw:860-866
###################################################
# continuous variables
don.B$log.netI <- log( ifelse(don.B$netIncome>0, don.B$netIncome, 0)+1 )
X.mtc <- c("age","rb090")
mix.3 <- mixed.mtc(data.rec=rec.A, data.don=don.B, match.vars=X.mtc,
                   y.rec="work", z.don="log.netI", 
                   method="MS")


