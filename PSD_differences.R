# run the LDPSA_starter.R script first
# M. Walsh, October 2014

# Load lab covariates -----------------------------------------------------

download("https://www.dropbox.com/s/gcga9uyt2b8cv9h/Lab_cov.csv.zip?dl=0", "Lab_cov.csv.zip", mode="wb")
unzip("Lab_cov.csv.zip", overwrite=T)
labcov <- read.table("Lab_cov.csv", header=T, sep=",")

# REML analyses ---------------------------------------------------

# Ultrasonic treatment differences of samples dispersed in water
water <- subset(ldps.comp, Disp=="water", select=c(Site, GID, SSN, Ultra, Depth, V1, V2))
water <- merge(water, labcov, by="SSN")

wV1.lmer <- lmer(V1~I(Depth/100)+Ultra*log(SOC)+(1|Site)+(1|GID:Site), data=water)
summary(wV1.lmer)
wV2.lmer <- lmer(V2~I(Depth/100)+Ultra*log(SOC)+(1|Site)+(1|GID:Site), data=water)
summary(wV2.lmer)

# Ultrasonic treatment differences of samples dispersed in calgon
calgon <- subset(ldps.comp, Disp=="calgon", select=c(Site, GID, SSN, Ultra, Depth, V1, V2))
calgon <- merge(calgon, labcov, by="SSN")

cV1.lmer <- lmer(V1~I(Depth/100)+Ultra*log(SOC)+(1|Site)+(1|GID:Site), data=calgon)
summary(cV1.lmer)
cV2.lmer <- lmer(V2~I(Depth/100)+Ultra*log(SOC)+(1|Site)+(1|GID:Site), data=calgon)
summary(cV2.lmer)

