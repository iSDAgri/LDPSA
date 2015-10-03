# PSD changes in water or calgon following ultra-sonification as ...
# potential indicators of soil aggregate stability/erodibility
# M. Walsh, October 2014

# Set your local working directory here, e.g.
dat_dir <- "/Users/markuswalsh/Documents/LDSF/LDPSA"
setwd(dat_dir)

# Load packages
require(downloader)
require(arm)

# Either run LDPSA_starter.R first or download the data -------------------
download("https://www.dropbox.com/s/8mxyz1jeacczhme/LDPSA_comp.csv.zip?dl=0", "LDPSA_comp.csv.zip", mode="wb")
unzip("LDPSA_comp.csv.zip", overwrite=T)
ldps.comp <- read.table("LDPSA_comp.csv", header=T, sep=",")

# Lab covariate data
# download("https://www.dropbox.com/s/gcga9uyt2b8cv9h/Lab_cov.csv.zip?dl=0", "Lab_cov.csv.zip", mode="wb")
# unzip("Lab_cov.csv.zip", overwrite=T)
# labcov <- read.table("Lab_cov.csv", header=T, sep=",")

# Data setup --------------------------------------------------------------
# Calculate ultra-sonification PSD differences between samples dispersed in water
w1 <- subset(ldps.comp, TRT=="w1", select=c(Site, GID, SSN, Depth, V1, V2))
colnames(w1) <- c("Site", "GID", "SSN", "Depth", "w1v1", "w1v2")
w4 <- subset(ldps.comp, TRT=="w4", select=c(SSN, V1, V2))
colnames(w4) <- c("SSN", "w4v1", "w4v2")
dwater <- merge(w1, w4, by="SSN")
dwater$dwV1 <- dwater$w4v1-dwater$w1v1
dwater$dwV2 <- dwater$w4v2-dwater$w1v2
vars <- c("dwV1","dwV2")
w <- dwater[vars]
Sw <- cov(w)
dwater$dw <- mahalanobis(w, c(quantile(dwater$dwV1, prob=0.975), quantile(dwater$dwV2, prob=0.975)), Sw)
dwater <- merge(dwater, labcov, by="SSN")

# Calculate ultra-sonification PSD differences between samples dispersed in calgon
c1 <- subset(ldps.comp, TRT=="c1", select=c(Site, GID, SSN, Depth, V1, V2))
colnames(c1) <- c("Site", "GID", "SSN", "Depth", "c1v1", "c1v2")
c4 <- subset(ldps.comp, TRT=="c4", select=c(SSN, V1, V2))
colnames(c4) <- c("SSN", "c4v1", "c4v2")
dcalgon <- merge(c1, c4, by="SSN")
dcalgon$dcV1 <- dcalgon$c4v1-dcalgon$c1v1
dcalgon$dcV2 <- dcalgon$c4v2-dcalgon$c1v2
vars <- c("dcV1","dcV2")
c <- dcalgon[vars]
Sc <- cov(c)
dcalgon$dc <- mahalanobis(c, c(quantile(dcalgon$dcV1, prob=0.975), quantile(dcalgon$dcV2, prob=0.975)), Sc)
dcalgon <- merge(dcalgon, labcov, by="SSN")

# Example REML models -----------------------------------------------------
# Main effects PSD models for ultra-sonification changes in water 
dwV1.lmer <- lmer(dwV1~Depth+log(SOC)+EC+pH+log(ECEC)+(1|Site)+(1|GID:Site), data=dwater)
summary(dwV1.lmer)
dwV2.lmer <- lmer(dwV2~Depth+log(SOC)+EC+pH+log(ECEC)+(1|Site)+(1|GID:Site), data=dwater)
summary(dwV2.lmer)
dw.lmer <- lmer(dw~Depth+log(SOC)+EC+pH+log(ECEC)+(1|Site)+(1|GID:Site), data=dwater)
summary(dw.lmer)

# Site-level random effects plot for samples dispersed in water
dwV1.ranef <- ranef(dwV1.lmer)
dwV2.ranef <- ranef(dwV2.lmer)
plot(dwV1.ranef$Site[,1], dwV2.ranef$Site[,1], type="n", xlim=c(-2,2), ylim=c(-1,1), xlab="delta(ilr[Sand|Silt,Clay])", ylab="delta(ilr[Silt|Clay])")
text(dwV1.ranef$Site[,1], dwV2.ranef$Site[,1], labels = row.names(dwV1.ranef$Site), cex=0.7)

# Main effects PSD models for ultra-sonification changes in calgon
dcV1.lmer <- lmer(dcV1~Depth+log(SOC)+EC+pH+log(ECEC)+(1|Site)+(1|GID:Site), data=dcalgon)
summary(dcV1.lmer)
dcV2.lmer <- lmer(dcV2~Depth+log(SOC)+EC+pH+log(ECEC)+(1|Site)+(1|GID:Site), data=dcalgon)
summary(dwV2.lmer)
dc.lmer <- lmer(dc~Depth+log(SOC)+EC+pH+log(ECEC)+(1|Site)+(1|GID:Site), data=dcalgon)
summary(dc.lmer)

# Site-level random effects plot for samples dispersed in calgon
dcV1.ranef <- ranef(dcV1.lmer)
dcV2.ranef <- ranef(dcV2.lmer)
plot(dcV1.ranef$Site[,1], dcV2.ranef$Site[,1], type="n", xlim=c(-2.5,2.5), ylim=c(-1.5,1.5), xlab="delta(ilr[Sand|Silt,Clay])", ylab="delta(ilr[Silt|Clay])")
text(dcV1.ranef$Site[,1], dcV2.ranef$Site[,1], labels = row.names(dcV1.ranef$Site), cex=0.7)
