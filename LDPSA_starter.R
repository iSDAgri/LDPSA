#' AfSIS laser diffraction particle size analyses (LDPSA) starter
#' M. Walsh, September 2014

# Load packages
require(downloader)
require(proj4)
require(raster)
require(compositions)
require(arm)
require(soiltexture)

#+ Data download -----------------------------------------------------------
# Set your local working directory here, e.g.
dir.create("LDPSA_data", showWarnings=F)
dat_dir <- "./LDPSA_data"

# LDPSA_60 lab data
download("https://www.dropbox.com/s/je066cib8zt9nh6/LDPSA_60.zip?dl=0", "./LDPSA_data/LDPSA_60.zip", mode="wb")
unzip("./LDPSA_data/LDPSA_60.zip", exdir="./LDPSA_data", overwrite=T)
profile <- read.table(paste(dat_dir, "/Profiles.csv", sep=""), header=T, sep=",")
profile <- na.omit(profile)
samples <- read.table(paste(dat_dir, "/Samples.csv", sep=""), header=T, sep=",")
samples <- na.omit(samples)

# Lab covariate data
download("https://www.dropbox.com/s/gcga9uyt2b8cv9h/Lab_cov.csv.zip?dl=0", "./LDPSA_data/Lab_cov.csv.zip", mode="wb")
unzip("./LDPSA_data/Lab_cov.csv.zip", exdir="./LDPSA_data", overwrite=T)
labcov <- read.table(paste(dat_dir, "/Lab_cov.csv", sep=""), header=T, sep=",")

#+ Generate coordinate reference and GID's ----------------------------------
# Project profile coords to Africa LAEA from LonLat
profile.laea <- as.data.frame(project(cbind(profile$Lon, profile$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(profile.laea) <- c("x","y")
profile <- cbind(profile, profile.laea)

# Generate AfSIS grid cell ID's (GID)
res.pixel <- 1000
xgid <- ceiling(abs(profile$x)/res.pixel)
ygid <- ceiling(abs(profile$y)/res.pixel)
gidx <- ifelse(profile$x<0, paste("W", xgid, sep=""), paste("E", xgid, sep=""))
gidy <- ifelse(profile$y<0, paste("S", ygid, sep=""), paste("N", ygid, sep=""))
GID <- paste(gidx, gidy, sep="-")
profile.gid <- cbind(profile, GID)
ldpsdat <- merge(profile.gid, samples, by="PID")

# Texture triangle plot examples ------------------------------------------
# Water dispersed samples
fractions <- c("Clay","Silt","Sand")
TRTw1  <- subset(ldpsdat, TRT=="w1", select=c(PID, Depth, SSN, Clay, Silt, Sand)) 
TRTw1.comp <- as.data.frame(acomp(TRTw1[fractions], total=100))
TT.plot(tri.data=TRTw1.comp, class.sys="HYPRES.TT", cex=0.6, cex.lab=1, cex.axis=0.8, main="Dispersed in water", css.names=fractions)

# Water dispersed samples, ultra-sonified for 4 min
TRTw4 <- subset(ldpsdat, TRT=="w4", select=c(PID, Depth, SSN, Clay, Silt, Sand))
TRTw4.comp <- as.data.frame(acomp(TRTw4[fractions], total=100))
TT.plot(tri.data=TRTw4.comp, class.sys="HYPRES.TT", cex=0.6, cex.lab=1, cex.axis=0.8, main="Dispersed in water + 4 min ultra-sonification", css.names=fractions)

# Integrated log ratio (ilr) transformation -------------------------------
# Define binary partion
cdata <- acomp(ldpsdat[fractions], total=100)
bpart <- t(matrix(c(-1,-1, 1,
                    -1, 1, 0), ncol=3, nrow=2, byrow=T))
CoDaDendrogram(X=acomp(cdata), signary=bpart)					
idata <- as.data.frame(ilr(cdata, V=bpart))
ldps.comp <- cbind(ldpsdat, idata)
# write.csv(ldps.comp, "LDPSA_comp.csv", row.name = F)

# Example REML analyses ---------------------------------------------------
# Main effects model ilr[Sand|Silt,Clay] = V1
V1.lmer <- lmer(V1~Disp*Ultra+I(Depth/100)+(1|Site)+(1|GID:Site), data=ldps.comp)
summary(V1.lmer)

# Extract and plot Site-level random effects and standard errors
V1.ranef <- ranef(V1.lmer)
V1.se <- se.coef(V1.lmer)
coefplot(V1.ranef$Site[,1], V1.se$Site[,1], varnames=rownames(V1.ranef$Site), xlim=c(-3,3), CI=2, cex.var=0.6, cex.pts=0.9, main="")

# Main effects model ilr[Silt|Clay] = V2
V2.lmer <- lmer(V2~Disp*Ultra+I(Depth/100)+(1|Site)+(1|GID:Site), data=ldps.comp)
summary(V2.lmer)
V2.ranef <- ranef(V2.lmer)
V2.se <- se.coef(V2.lmer)
coefplot(V2.ranef$Site[,1], V2.se$Site[,1], varnames=rownames(V2.ranef$Site), xlim=c(-1,1), CI=2, cex.var=0.6, cex.pts=0.9, main="")

# Ultrasonic treatment differences of samples dispersed in water with covariates
water <- subset(ldps.comp, Disp=="water", select=c(Site, GID, SSN, Ultra, Depth, V1, V2))
water <- merge(water, labcov, by="SSN")

wV1.lmer <- lmer(V1~Ultra+I(Depth/100)+log(SOC)+(1|Site)+(1|GID:Site), data=water)
summary(wV1.lmer)
wV1.ranef <- ranef(wV1.lmer)
wV1.se <- se.coef(wV1.lmer)
coefplot(wV1.ranef$Site[,1], wV1.se$Site[,1], varnames=rownames(wV1.ranef$Site), xlim=c(-2.5,2.5), CI=2, cex.var=0.6, cex.pts=0.9, main="")

wV2.lmer <- lmer(V2~Ultra+I(Depth/100)+log(SOC)+(1|Site)+(1|GID:Site), data=water)
summary(wV2.lmer)
wV2.ranef <- ranef(wV2.lmer)
wV2.se <- se.coef(wV2.lmer)
coefplot(wV2.ranef$Site[,1], wV2.se$Site[,1], varnames=rownames(wV2.ranef$Site), xlim=c(-1,1), CI=2, cex.var=0.6, cex.pts=0.9, main="")

# Site-level random effects plot for samples dispersed in water
plot(wV1.ranef$Site[,1], wV2.ranef$Site[,1], type="n", xlim=c(-3,3), ylim=c(-1.5,1.5), xlab="ilr[Sand|Silt,Clay]", ylab="ilr[Silt|Clay]")
text(wV1.ranef$Site[,1], wV2.ranef$Site[,1], labels = row.names(wV1.ranef$Site), cex=0.7)

# Ultrasonic treatment differences of samples dispersed in calgon with covariates
calgon <- subset(ldps.comp, Disp=="calgon", select=c(Site, GID, SSN, Ultra, Depth, V1, V2))
calgon <- merge(calgon, labcov, by="SSN")

cV1.lmer <- lmer(V1~Ultra+I(Depth/100)+log(SOC)+(1|Site)+(1|GID:Site), data=calgon)
summary(cV1.lmer)
cV1.ranef <- ranef(cV1.lmer)
cV1.se <- se.coef(cV2.lmer)
coefplot(cV1.ranef$Site[,1], cV1.se$Site[,1], varnames=rownames(cV1.ranef$Site), xlim=c(-2,2), CI=2, cex.var=0.6, cex.pts=0.9, main="")

cV2.lmer <- lmer(V2~Ultra+I(Depth/100)+log(SOC)+(1|Site)+(1|GID:Site), data=calgon)
summary(cV2.lmer)
cV2.ranef <- ranef(cV2.lmer)
cV2.se <- se.coef(cV2.lmer)
coefplot(cV2.ranef$Site[,1], cV2.se$Site[,1], varnames=rownames(cV2.ranef$Site), xlim=c(-1.5,1.5), CI=2, cex.var=0.6, cex.pts=0.9, main="")

# Site-level random effects plot for samples dispersed in calgon
plot(cV1.ranef$Site[,1], cV2.ranef$Site[,1], type="n", xlim=c(-3,3), ylim=c(-1.5,1.5), xlab="ilr[Sand|Silt,Clay]", ylab="ilr[Silt|Clay]")
text(cV1.ranef$Site[,1], cV2.ranef$Site[,1], labels = row.names(cV1.ranef$Site), cex=0.7)
