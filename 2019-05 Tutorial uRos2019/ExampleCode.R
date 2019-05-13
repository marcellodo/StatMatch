# Example of Statistical matching in R
# using StatMatch


# install.packages("StatMatch")
library(StatMatch)

data("samp.A")
str(samp.A)

data("samp.B")
str(samp.B)

# variables with same names in both samp.A and samp.B
intersect(colnames(samp.A), colnames(samp.B)) 

# variables only in samp.A (Y)
setdiff(colnames(samp.A), colnames(samp.B))
summary(samp.A$n.income)
summary(samp.A$c.neti)

# variables only in samp.B (Z)
setdiff(colnames(samp.B), colnames(samp.A))
summary(samp.B$labour5)


###########################
# check target population

# estimated population size
sum(samp.A$ww) # sum of survey weights
sum(samp.B$ww) # sum of survey weights

# distribution by regions
ttA <- xtabs(ww~area5, data=samp.A)
ttA
ttB <- xtabs(ww~area5, data=samp.B)
ttB
cbind(A=prop.table(ttA),
      B=prop.table(ttB))

# measure closeness between distributions
comp.prop(p1 = ttA, p2 = ttB, 
          n1 = nrow(samp.A), n2 = nrow(samp.B), ref = F)


# distribution by gender
ttA <- xtabs(ww~sex, data=samp.A)
ttB <- xtabs(ww~sex, data=samp.B)
comp.prop(p1 = ttA, p2 = ttB, 
          n1 = nrow(samp.A), n2 = nrow(samp.B), ref = F)


# distribution by age-classes
ttA <- xtabs(ww~c.age, data=samp.A)
ttB <- xtabs(ww~c.age, data=samp.B)
comp.prop(p1 = ttA, p2 = ttB, 
          n1 = nrow(samp.A), n2 = nrow(samp.B), ref = F)

# distribution by (region x gender x age-classes)
ttA <- xtabs(ww~sex+c.age+area5, data=samp.A)
ttB <- xtabs(ww~sex+c.age+area5, data=samp.B)
comp.prop(p1 = ttA, p2 = ttB, 
          n1 = nrow(samp.A), n2 = nrow(samp.B), ref = F)


# measure closeness between distributions
comp.prop(p1 = ttA, p2 = ttB, 
          n1 = nrow(samp.A), n2 = nrow(samp.B), ref = F)


###########################
# check marginal distr. of common variables
ttA <- xtabs(ww~hsize5, data=samp.A)
ttB <- xtabs(ww~hsize5, data=samp.B)
comp.prop(p1 = ttA, p2 = ttB, 
          n1 = nrow(samp.A), n2 = nrow(samp.B), ref = F)

ttA <- xtabs(ww~urb, data=samp.A)
ttB <- xtabs(ww~urb, data=samp.B)
comp.prop(p1 = ttA, p2 = ttB, 
          n1 = nrow(samp.A), n2 = nrow(samp.B), ref = F)


ttA <- xtabs(ww~edu7, data=samp.A)
ttB <- xtabs(ww~edu7, data=samp.B)
cc <- comp.prop(p1 = ttA, p2 = ttB, 
                n1 = nrow(samp.A), n2 = nrow(samp.B), 
                ref = F)
cc$meas

##############################################
# best predictors of n.income (Y, is continuous)

# install.packages("Hmisc")
library(Hmisc)
spearman2(n.income~area5+urb+hsize5+c.age+sex+marital+edu7, 
          data=samp.A)

# best predictors of labour5 (Z, is categorical)

pws <- pw.assoc(labour5~area5+urb+hsize5+c.age+sex+marital+edu7, 
                data=samp.B, out.df = TRUE)

pws[, c("norm.mi", "U", "AIC", "npar")]

# matching variables
x_y <- c("c.age", "sex", "edu7")
x_z <- c("c.age", "sex", "marital","edu7")
intersect(x_y, x_z)
union(x_y, x_z)


################################################
## selecting matching variables by uncertainty

# rescale weights to sum up to n
wwA <- samp.A$ww / sum(samp.A$ww) * nrow(samp.A)
wwB <- samp.B$ww / sum(samp.B$ww) * nrow(samp.B)

#estimate joint ditribution of starting Xs
txA <- xtabs(wwA~area5+urb+hsize5+c.age+sex+marital+edu7,
             data=samp.A)

txB <- xtabs(wwB~area5+urb+hsize5+c.age+sex+marital+edu7,
             data=samp.B)
txx <- txA+txB

# joint X vs. Y
txyA <- xtabs(wwA~area5+urb+hsize5+c.age+sex+marital+edu7+c.neti,
             data=samp.A)

# joint X vs. Z
txzB <- xtabs(wwB~area5+urb+hsize5+c.age+sex+marital+edu7+labour5,
              data=samp.B)

# launch search
unc <- selMtc.by.unc(tab.x=txx, tab.xy=txyA,  tab.xz=txzB, 
                     corr.d=2)

unc$ini.ord
unc$av.df


#################################################################
# assessment of uncertainty for categorical variables

txA <- xtabs(wwA~area5+c.age+sex+edu7,
             data=samp.A)

txB <- xtabs(wwB~area5+c.age+sex+edu7,
             data=samp.B)
txx <- txA+txB

# joint X vs. Y
txyA <- xtabs(wwA~area5+c.age+sex+edu7+c.neti,
              data=samp.A)

# joint X vs. Z
txzB <- xtabs(wwB~area5+c.age+sex+edu7+labour5,
              data=samp.B)

# estimate frechet-bonferroni bounds for relative frequancies in table
# c.neti vs. labour5
fbw <- Frechet.bounds.cat(tab.x = txx, 
                          tab.xy = txyA, tab.xz = txzB, 
                          print.f = "data.frame", align.margins =TRUE)
head(fbw$bounds, 4)

#####################################################
## Random Hot-deck within fixed classes

# check for empty classes in donor
dcA <- xtabs(~c.age+sex+edu7+area5, data=samp.A)
dcB <- xtabs(~c.age+sex+edu7+area5, data=samp.B)

tst <- dcA>0 & dcB==0
sum(tst)

# discard area5
dcA <- xtabs(~c.age+sex+edu7, data=samp.A)
dcB <- xtabs(~c.age+sex+edu7, data=samp.B)

tst <- dcA>0 & dcB==0
sum(tst)

# discard edu7

dcA <- xtabs(~c.age+sex, data=samp.A)
dcB <- xtabs(~c.age+sex, data=samp.B)

tst <- dcA>0 & dcB==0
sum(tst)

# run random hot deck
out.rnd <- RANDwNND.hotdeck(data.rec = samp.A, data.don = samp.B,
                       don.class = c("c.age", "sex")) 

head(out.rnd$mtc.ids, 4)

# create synthetic data set, samp.A is the recipient
fillA.rnd.1 <- create.fused(data.rec = samp.A, data.don = samp.B,
                          mtc.ids = out.rnd$mtc.ids, dup.x =T,
                          match.vars =  c("c.age", "sex"), 
                          z.vars = "labour5")

head(fillA.rnd.1, 2)


#####################################################
## Random Hot-deck within NON-fixed classes
# random selection of one of k=5 closest donors in term of age
# having same gender and education level
out.rnd2 <- RANDwNND.hotdeck(data.rec = samp.A, data.don = samp.B,
                             don.class = c("edu7", "sex"), 
                             match.vars = "age", cut.don = "exact", 
                             k = 5)

head(out.rnd2$sum.dist, 4)

# create synthetic data set, samp.A is the recipient
fillA.rnd.2 <- create.fused(data.rec = samp.A, data.don = samp.B,
                            mtc.ids = out.rnd2$mtc.ids, dup.x =T,
                            match.vars =  c("c.age", "sex"), 
                            z.vars = "labour5")

#####################################################
## Nearest neighbour distance hot-deck 
# within classes formed on sex
# distance calculated on age
# unconstrained: a donor can be used more than once

out.nnd1 <- NND.hotdeck(data.rec = samp.A, data.don = samp.B,
                        don.class = "sex", 
                        match.vars = "age")

summary(out.nnd1$dist.rd)


# create synthetic data set, samp.A is the recipient
fillA.nnd.1 <- create.fused(data.rec = samp.A, data.don = samp.B,
                            mtc.ids = out.nnd1$mtc.ids, dup.x =T,
                            match.vars =  c("age", "sex"), 
                            z.vars = "labour5")
head(fillA.nnd.1, 3)

# constrained: a donor can be used just once (k=1)

out.nnd2 <- NND.hotdeck(data.rec = samp.A, data.don = samp.B,
                        don.class = "sex", 
                        match.vars = "age", 
                        constrained = T, k = 1,
                        constr.alg = "hungarian")

sum(out.nnd2$dist.rd)
sum(out.nnd1$dist.rd)


# create synthetic data set, samp.A is the recipient
fillA.nnd.2 <- create.fused(data.rec = samp.A, data.don = samp.B,
                            mtc.ids = out.nnd2$mtc.ids, dup.x =T,
                            match.vars =  c("age", "sex"), 
                            z.vars = "labour5")
head(fillA.nnd.2, 3)

######################################################################
# compare marginal distribution of imputed income vs. reference one
# unweighted
t.imp <- table(fillA.rnd.1$labour5)
t.don <- table(samp.B$labour5)

ck1 <- comp.prop(p1 = t.imp, p2 = t.don, 
                 n1 = nrow(fillA.rnd.1),
                 n2 = nrow(samp.B), 
                 ref = T)
ck1$meas
          
# weighted
t.imp <- xtabs(ww~labour5, data = fillA.rnd.1)
t.don <- xtabs(ww~labour5, data = samp.B)

ck1.w <- comp.prop(p1 = t.imp, p2 = t.don, 
                   n1 = nrow(fillA.rnd.1),
                   n2 = nrow(samp.B), 
                   ref = T)
ck1.w$meas

# check joint distribution of imputed income with matching variables
# unweighted
t.imp <- xtabs(~labour5+c.age+sex, data = fillA.rnd.1)
t.don <- xtabs(~labour5+c.age+sex, data = samp.B)

cc <- comp.prop(p1 = t.imp, p2 = t.don, 
                n1 = nrow(fillA.rnd.1),
                n2 = nrow(samp.B), 
                ref = T)
cc$meas

# weighted
t.imp <- xtabs(ww~labour5+c.age+sex, data = fillA.rnd.1)
t.don <- xtabs(ww~labour5+c.age+sex, data = samp.B)

cc <- comp.prop(p1 = t.imp, p2 = t.don, 
                n1 = nrow(fillA.rnd.1),
                n2 = nrow(samp.B), 
                ref = T)
cc$meas

#######################################
# comparison of methods

# marginal of imputed (weighted)
t.don <- xtabs(ww~labour5, data = samp.B)
t.imp.rnd1 <- xtabs(ww~labour5, data = fillA.rnd.1)
t.imp.rnd2 <- xtabs(ww~labour5, data = fillA.rnd.2)
t.imp.nndu <- xtabs(ww~labour5, data = fillA.nnd.1)
t.imp.nndc <- xtabs(ww~labour5, data = fillA.nnd.2)

a <- comp.prop(p1 = t.imp.rnd1, p2 = t.don, 
               n1 = nrow(fillA.rnd.1),
               n2 = nrow(samp.B), 
               ref = T)
b <- comp.prop(p1 = t.imp.rnd2, p2 = t.don, 
               n1 = nrow(fillA.rnd.1),
               n2 = nrow(samp.B), 
               ref = T)

c <- comp.prop(p1 = t.imp.nndu, p2 = t.don, 
               n1 = nrow(fillA.rnd.1),
               n2 = nrow(samp.B), 
               ref = T)

d <- comp.prop(p1 = t.imp.nndc, p2 = t.don, 
               n1 = nrow(fillA.rnd.1),
               n2 = nrow(samp.B), 
               ref = T)
rbind(rnd.1=a$meas, rnd.2=b$meas,
      nnd.unc=c$meas, nnd.c=d$meas)

# Joint of imputed vs. matching (weighted)
t.don <- xtabs(ww~labour5+c.age+sex, data = samp.B)
t.imp.rnd1 <- xtabs(ww~labour5+c.age+sex, data = fillA.rnd.1)
t.imp.rnd2 <- xtabs(ww~labour5+c.age+sex, data = fillA.rnd.2)
t.imp.nndu <- xtabs(ww~labour5+c.age+sex, data = fillA.nnd.1)
t.imp.nndc <- xtabs(ww~labour5+c.age+sex, data = fillA.nnd.2)
a <- comp.prop(p1 = t.imp.rnd1, p2 = t.don, 
               n1 = nrow(fillA.rnd.1),
               n2 = nrow(samp.B), 
               ref = T)
b <- comp.prop(p1 = t.imp.rnd2, p2 = t.don, 
               n1 = nrow(fillA.rnd.1),
               n2 = nrow(samp.B), 
               ref = T)

c <- comp.prop(p1 = t.imp.nndu, p2 = t.don, 
               n1 = nrow(fillA.rnd.1),
               n2 = nrow(samp.B), 
               ref = T)

d <- comp.prop(p1 = t.imp.nndc, p2 = t.don, 
               n1 = nrow(fillA.rnd.1),
               n2 = nrow(samp.B), 
               ref = T)
rbind(rnd.1=a$meas, rnd.2=b$meas,
      nnd.unc=c$meas, nnd.c=d$meas)


#######################################à
# estimate table Y vs. Z in filled A
# based on the chosen method: NND unconstrained

t.yz <- xtabs(ww~c.neti+labour5, data = fillA.nnd.1)

round(addmargins( prop.table(t.yz))*100,2)

# estimate association
assoc <- pw.assoc(c.neti~labour5, data = fillA.nnd.1, 
                  out.df = T, weights = "ww")
assoc$V
assoc$norm.mi



