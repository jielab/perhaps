library(data.table); library(dplyr); library(tidyverse); library(magrittr); library(ggpubr);  
setwd("C:/Users/黄捷/Desktop")
ukbdir="D:/projects/001UKB/"

## e1: C-T;  e2: T-T;  e3: T-C;  e4: C-C

### for those samples with WES ###
covars <- read.table(paste0(ukbdir,"pheno/raw/ukb.cov"), header=T, as.is=T) %>%
	select(IID, ethnicity, age, sex)
garray <- read.table(paste0(ukbdir,"pheno/raw/ukb_sqc_v2.txt"), header=T, as.is=T) %>%
	select(IID, genotyping.array) %>%
	rename(garray=genotyping.array)
typed = read.table(gzfile(paste0(ukbdir,"typed/typed.raw.gz"),"r"), header=T, as.is=T) %>% 
	select(IID, apoe.rs429358_T, apoe.rs7412_C)
phased <- read.table(paste0(ukbdir,"hap/apoe.txt"), header=T, as.is=T) %>%
	rename(apoe.phased=apoe)
wes <- read.table("D:/projects/001students/001Perhaps/results/apoe.wes.txt", header=T, as.is=T) %>%
	mutate (
	hap2new=ifelse(hap2=="---", hap1, hap2),  
	hap1r=ifelse(hap1=="C-T", "e1", ifelse(hap1=="T-T", "e2", ifelse(hap1=="T-C", "e3", ifelse(hap1=="C-C", "e4", NA)))),
	hap2r=ifelse(hap2new=="C-T", "e1", ifelse(hap2new=="T-T", "e2", ifelse(hap2new=="T-C", "e3", ifelse(hap2new=="C-C", "e4", NA)))),
	apoe=ifelse(hap1r<hap2r, paste0(hap1r, hap2r), paste0(hap2r, hap1r))
	)
dat0 <- Reduce(function(x,y) merge(x,y,by="IID"), list(covars, garray, typed, phased, wes)) 
dat <- dat0 %>% 
	filter (!is.na(apoe)) %>%
	mutate( concordant= ifelse(apoe==apoe.phased, "Yes","No"),
	cnt_min = ifelse(cnt2==0, cnt1, pmin(cnt1,cnt2)),
	cnt_maj = ifelse(cnt2==0, cnt1, pmax(cnt1,cnt2)),
	cnt_tot = cnt1 + cnt2 )
subset(dat, grepl("e1",apoe) | grepl("e1",apoe.phased))
ftable(dat$garray, dat$apoe.phased, dat$apoe, dnn=c("garray","phased","wes"))
aggregate(cnt_min ~ concordant + apoe, data=dat, mean)
dat$apoe <- factor(dat$apoe, levels=c("e1e2","e1e3","e1e4","e2e2","e2e3","e2e4","e3e3","e3e4","e4e4"), label=c("*1/*2","*1/*3","*1/*4","*2/*2","*2/*3","*2*/*4","*3/*3","*3/*4","*4/*4"))
png("figure2.png", w=5000, h=5000, res=600)
p <- ggboxplot(dat, x="apoe", y="cnt_min", color="concordant", notch=F, ylim=c(0,50), xlab="APOE derived from WES paired reads", ylab="Paired reads of the minor haplotype", outlier.shape=NA, font.label=list(size=100, face="bold"), size=1.5, palette="jco", add="none") # outlier.shape=NA,
p + font("xlab", size=16) + font("ylab", size=16) + font("xy.text", size=16, face="bold") +
	stat_compare_means(aes(label=..p.format.., group=concordant), method="wilcox.test", label.y=50)
dev.off()


### for all samples with imputed APOE ###
dat0 <- readRDS(paste0(ukbdir,"Rdata/ukb.phe.rds") %>% 
rename(apoe=hap.apoe, HDL=bc_HDL, LDL=bc_LDL, TC=bc_TC, TG=bc_TG)
dat <- subset(dat0, !is.na(apoe) & !grepl("e1", apoe) & grepl("White|Asian|Black",race))
dat <- subset(dat, select=c("IID","race","age","sex","apoe","HDL","LDL","TC","TG"))
table(dat$apoe, dat$race, useNA="always", dnn=c("apoe", "race"))
aggregate(HDL ~ apoe + race, data=dat, mean)
pdf("lipids.pdf", w=10, h=10)
ggline(dat, x="apoe", y=c("HDL","LDL","TC","TG"), xlab="APOE haplotype", ylab="", combine=T, add="mean_se", color="race", palette="jco", scales="free_y") +rotate_x_text(45)
dev.off()


### 
myhist <- hist(dat$age, breaks=seq(35,75,5))
avg <- by(dat$e4, dat$age_cat, function(x) mean(x,na.rm=T)); avg[is.na(avg)] <- 0
sdev <- by(dat$e4, dat$age_cat, sd); sdev[is.na(sdev)] <- 0
par(new=T)
plot(myhist$mids, avg, ylim=range(c(avg-2*sdev, avg+2*sdev)), pch=16, axes=F, xlab=NA, ylab=NA, cex=1.2)
#arrows(myhist$mids, avg-sdev, myhist$mids, avg+sdev, length=0.05, angle=90, code=3)
axis(side=4); mtext(side=4, line=3, 'measured')


### list by frequency ###
ggplot(dat, aes(fill=age_cat, y=e4, x=race)) + geom_bar(position="dodge", stat="identity")
bp <- boxplot(e4 ~ race*age_cat, data=dat, xlab="", ylab="", main="", las=2, col=rainbow(6), font=2)
aggregate(bc_LDL ~ race*apoe, data=dat, FUN=function(x) {round(c(length(x), mean(x), sd(x)), 2)} )
tbl <- table(m$ad01, m$apoe, m$smoke01, dnn=c("AD", "apoe","smoke")); tbl
round(prop.table(tbl, 2), 3)
aggregate(bc_LDL ~ e4*sex, data=dat, FUN=function(x) { round( c(length(x), mean(x), sd(x), quantile(x,probs=c(0,0.5,1))), 2) } )
bp <- boxplot(bc_LDL ~ e4*sex, data=dat, main="", las=2, col=rainbow(3), font=2); bp$stats
summary(lmod <- lm(bc_LDL ~ apoe, data=dat))
anova(lmod)
