library(devtools)
library(Biobase)
#############Q1##############
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
no_tf = edata
tf1 = log2(edata + 1)
tf2 = tf1 - rowMeans(tf1)
svd1=svd(no_tf)
svd2=svd(tf1)
svd3=svd(tf2)

par(pch = 19)
par(mfrow=c(1,3))

# NO Transformations
plot(svd1$d^2 / sum(svd1$d^2),
     main = "%Var w NO Transform",
     ylim = c(0, 1),
     ylab = "Percent Variance Explained",
     col = 2)
# Log2 Transformation
plot(svd2$d^2 / sum(svd2$d^2),
     main = "%Var w Log2 Transform",
     ylim = c(0, 1),
     ylab = "Percent Variance Explained",
     col = 2)

# Log2 Transformation and Subtract rowMeans
plot(svd3$d^2 / sum(svd3$d^2),
     main = "%Var w Log2 & Subtract rowMeans",
     ylim = c(0, 1),
     ylab = "Percent Variance Explained",
     col = 2)



###q2####
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
edata = log2(edata +1)
edata = edata - rowMeans(edata)
set.seed(333)
kmeans1 = kmeans(t(edata),centers=2)
svd1 = svd(edata)
q2 = cor(kmeans1$cluster, svd1$v[,1])
boxplot(svd1$v[,1] ~ kmeans1$cluster, border=c(1,2))
points(svd1$v[,1] ~ jitter(kmeans1$cluster),col=kmeans1$cluster)


#####Q3#########
library(broom)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)
gene1 = as.matrix(edata[1,])
mod1 = model.matrix(~ as.factor(pdata_bm$num.tech.reps))
lm1 = lm(gene1 ~ mod1)
plot(pdata_bm$num.tech.reps, gene1, col=2)


#####Q4#######
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)
gene1 = as.matrix(edata[1,])
lm1 = lm(gene1 ~ pdata_bm$age + pdata_bm$gender)


####Q5#####
library(limma)
library(edgeR)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
edata = log2(edata+1)
mod1 = model.matrix(~ pdata$population)
fit1= lm.fit(mod1, t(edata))

par(mfrow = c(1, 2))
hist(fit1$coefficients[1, ], breaks = 10, col = 2, ylim = c(0, 50000), xlab = "Intercept")
hist(fit1$coefficients[2, ], breaks = 10, col = 2, ylim = c(0, 50000), xlab = "PopulationYRI")
abline(v = 0, lwd = 2, col = 1)


#######Q6#######
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
edata = log2(edata+1)
mod1 = model.matrix(~ pdata$population)
fit1= lm.fit(mod1, t(edata))


#####Q7######
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

p_1 = pdata_bm[!is.na(pdata_bm$age),]
e_1 = edata[,!is.na(pdata_bm$age)]
mod1 = model.matrix(~ p_1$age)
fit1= lm.fit(mod1, t(e_1))
q7=fit1$coefficients[,1000][2]
plot(p_1$age,e_1[1000,],col=2)
abline(fit1$coefficients[,1000],col=1,lwd=3)



#####Q8######
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

p_1 = pdata_bm[!is.na(pdata_bm$age),]
e_1 = edata[,!is.na(pdata_bm$age)]
mod1 = model.matrix(~ p_1$age + as.factor(p_1$tissue.type))
fit1= lm.fit(mod1, t(e_1))
fit1$coefficients[,1]

##Q9###
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
table(pdata$population,pdata$study)


####Q10#######
library(sva)
library(snpStats)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)
edata = log2(edata+1)
edata=edata[!rowMeans(edata)<1,]
p_1 = pdata_bm[!is.na(pdata_bm$age),]
e_1 = edata[,!is.na(pdata_bm$age)]
set.seed(33353)

mod_age = model.matrix(~age,data=p_1)
mod0 = model.matrix(~1, data=p_1)
sva1 = sva(e_1,mod_age,mod0,n.sv=2)