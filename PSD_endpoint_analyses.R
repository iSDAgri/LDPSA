# PSD changes in water or calgon following ultra-sonification endpoints as ...
# potential indicators of soil aggregate stability/erodibility
# M. Walsh, October 2014

# Set your local working directory here, e.g.
dat_dir <- "/Users/markuswalsh/Documents/LDSF/LDPSA"
setwd(dat_dir)

# Load packages
require(downloader)
require(arm)

# Either run LDPSA_starter.R first or download the data -------------------

# LDPSA data
# download("https://www.dropbox.com/s/8mxyz1jeacczhme/LDPSA_comp.csv.zip?dl=0", "LDPSA_comp.csv.zip", mode="wb")
# unzip("LDPSA_comp.csv.zip", overwrite=T)
# ldps.comp <- read.table("LDPSA_comp.csv". header=T, sep=",")

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
attach(dwater)
dwater$dwV1 <- w4v1-w1v1
dwater$dwV2 <- w4v2-w1v2
dwater$dw <- w4v1-w1v1+w4v2-w1v2
detach()
dwater <- merge(dwater, labcov, by="SSN")

# Calculate ultra-sonification PSD differences between samples dispersed in calgon
c1 <- subset(ldps.comp, TRT=="c1", select=c(Site, GID, SSN, Depth, V1, V2))
colnames(c1) <- c("Site", "GID", "SSN", "Depth", "c1v1", "c1v2")
c4 <- subset(ldps.comp, TRT=="c4", select=c(SSN, V1, V2))
colnames(c4) <- c("SSN", "c4v1", "c4v2")
dcalgon <- merge(c1, c4, by="SSN")
attach(dcalgon)
dcalgon$dcV1 <- c4v1-c1v1
dcalgon$dcV2 <- c4v2-c1v2
dcalgon$dc <- c4v1-c1v1+c4v2-c1v2
detach()
dcalgon <- merge(dcalgon, labcov, by="SSN")

# Example REML models -----------------------------------------------------

# Main effects PSD models for ultra-sonification changes in water 
dwV1.lmer <- lmer(dwV1~I(Depth/100)+log(SOC)+EC+pH+ECEC+(1|Site)+(1|GID:Site), data=dwater)
summary(dwV1.lmer)

dwV2.lmer <- lmer(dwV2~I(Depth/100)+log(SOC)+EC+pH+ECEC+(1|Site)+(1|GID:Site), data=dwater)
summary(dwV2.lmer)

dw.lmer <- lmer(dw~I(Depth/100)+log(SOC)+EC+pH+ECEC+(1|Site)+(1|GID:Site), data=dwater)
summary(dw.lmer)
