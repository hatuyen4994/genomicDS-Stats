####Q1
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

sub.10 <- snps.10[,use]
snp3 = as.numeric(snpdata[,3])
snp3[snp3==0] = NA
glm3 = glm(status ~ snp3,family="binomial")
tidy(glm3)

lm3 = lm(status ~ snp3)
tidy(lm3)

###Q3####
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

sub.10 <- snps.10[,use]
snp10 = as.numeric(snpdata[,10])
snp10[snp10==0] = NA
#additive model
glm10 = glm(status ~ snp10,family="binomial")
tidy(glm10)
table(status[!is.na(snp10)],glm10$fitted.values)


###recessive model
snp10_dom = (snp10 == 2)
glm10_dom = glm(status ~ snp10_dom,family="binomial")
tidy(glm10_dom)
table(status[!is.na(snp10_dom)],glm10_dom$fitted.values)


###Q4####
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

results = rep(NA, dim(snpdata)[2])
for (i in 1:ncol(snpdata)){
  snp_i = as.numeric(snpdata[,i])
  snp_i[snp_i ==0] = NA
  glm_i = glm(status~snp_i, family = "binomial")
  results[i] = tidy(glm_i)$statistic[2]
}
mean(results)
min(results)
max(results)

###Q5####
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc
glm_all = snp.rhs.tests(status ~ 1,snp.data=sub.10)

results = rep(NA, dim(snpdata)[2])
for (i in 1:ncol(snpdata)){
  snp_i = as.numeric(snpdata[,i])
  snp_i[snp_i ==0] = NA
  glm_i = glm(status~snp_i, family = "binomial")
  results[i] = tidy(glm_i)$statistic[2]
}
cor(glm_all@chisq,results^2)


###Q6####
library(devtools)
library(Biobase)
library(limma)
library(edge)
library(genefilter)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

edata = log2(as.matrix(edata) + 1)

# perform rowttests
tstats_obj = rowttests(edata, as.factor(pdata$population))
tidy(tstats_obj)

fstats_obj = rowFtests(edata, as.factor(pdata$population))
tidy(fstats_obj)

par(mfrow=c(1,2))
hist(tstats_obj$statistic, col=2)
hist(fstats_obj$statistic, col=2)

###Q7###
library(DESeq2)
library(limma)
library(edge)
library(genefilter)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
edata = edata[rowMeans(edata) > 100,]
fdata = fData(mp)

# using DESeq2 test the differences between the studies
de = DESeqDataSetFromMatrix(edata, pdata, ~study)
glm_de = DESeq(de)
result_de = results(glm_de)

# using limma test the differences
edata = log2(as.matrix(edata) + 1)
mod = model.matrix(~ as.factor(pdata$study))
fit_limma = lmFit(edata, mod)
ebayes_limma = eBayes(fit_limma) 
top = topTable(ebayes_limma,number=dim(edata)[1], sort.by="none")

# correlation in the statistics between two analyses
cor(result_de$stat, top$t)
y = cbind(result_de$stat, top$t)
limma::plotMA(y)


###Q8#####
# DESeq analysis
fp_bh = p.adjust(result_de$pvalue, method="BH")
sum(fp_bh < 0.05)
# limma analysis
fp_bh = p.adjust(top$P.Value, method="BH")
sum(fp_bh < 0.05)
