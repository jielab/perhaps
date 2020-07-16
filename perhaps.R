library(data.table); library(dplyr)
library(tidyverse); library(magrittr) 
library(ggpubr)

ukbdir="D:/projects/001UKB/"

### for those samples with WES ###

covars <- read.table(paste0(ukbdir,"pheno/raw/ukb.cov"), header=T, as.is=T) %>%
	select(IID, ethnicity, age, sex)

garray <- read.table(paste0(ukbdir,"pheno/raw/ukb_sqc_v2.txt"), header=T, as.is=T) %>%
	select(IID, genotyping.array) %>%
	rename(garray=genotyping.array)

typed = read.table(gzfile(paste0(ukbdir,"typed/typed.raw.gz),"r"), header=T, as.is=T) %>% 
	select(IID, apoe.rs429358_T, apoe.rs7412_C)

phased <- read.table(paste0(ukbdir,"hap/apoe.hap.txt"), header=T, as.is=T) %>%
	rename(apoe.phased=apoe)

wes <- read.table(paste0(ukbdir,"wes/apoe.wes.txt"), header=T, as.is=T) 

dat <- Reduce(function(x,y) merge(x,y,by="IID"), list(covars, garray, typed, phased, wes)) %>% 
	mutate( same= ifelse(apoe==apoe.phased, "Y","N"),
	cnt_min = ifelse(hap2=="---", cnt1, pmin(cnt1,cnt2)),
	cnt_maj = ifelse(hap2=="---", cnt1, pmax(cnt1,cnt2)),
	cnt_tot = cnt1 + cnt2 )

ftable(dat$garray, dat$apoe.phased, dat$apoe, dnn=c("garray","phased","wes"))

dat1 <- subset(dat, grepl("e1", apoe)) 
aggregate(cnt_min ~ same + apoe, data=dat, mean)

pdf("fig2.pdf", w=8, h=8)
p <- ggboxplot(dat, x="apoe", y="cnt_min", color="same", notch=F, ylim=c(0,50), xlab="APOE derived from WES paired reads", ylab="Paired reads of the minor haplotype", outlier.shape=NA, font.label=list(size=100, face="bold"), size=1.5, palette="jco", add="none") # outlier.shape=NA,
p + stat_compare_means(aes(label=..p.format.., group=same), method="wilcox.test", label.y=50)
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
