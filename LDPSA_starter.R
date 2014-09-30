# AfSIS laser diffraction particle size analyses starter (LDPSA)
# M. Walsh, September 2014

# Set your local working directory here, e.g.
dat_dir <- "/Users/markuswalsh/Documents/LDSF/LDPSA"
setwd(dat_dir)

# Load packages
require(downloader)
require(proj4)
require(raster)
require(compositions)
require(arm)
require(soiltexture)

# Load data and na.omit samples w missing values ---------------------------

# LDPSA_60 lab data
download("https://www.dropbox.com/s/je066cib8zt9nh6/LDPSA_60.zip?dl=0", "LDPSA_60.zip", mode="wb")
unzip("LDPSA_60.zip", overwrite=T)
profile <- read.table("Profiles.csv", header=T, sep=",")
profile <- na.omit(profile)
samples <- read.table("Samples.csv", header=T, sep=",")
samples <- na.omit(samples)

# Generate coordinate reference and GID's ----------------------------------

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

# Water dispersed topsoils
TRTw1  <- subset(ldpsdat, TRT=="w1" & Depth=="Topsoil", select=c(Site, Clay, Silt, Sand)) 
fractions <- c("Clay","Silt","Sand")
TRTw1.comp <- as.data.frame(acomp(TRTw1[fractions], total=100))
TT.plot(tri.data=TRTw1.comp, class.sys="HYPRES.TT", cex=0.6, cex.lab=1, cex.axis=0.8, main="Topsoils dispersed in water", css.names=fractions)

# Water dispersed topsoils ultra-sonified for 4 min
TRTw4 <- subset(ldpsdat, TRT=="w4" & Depth=="Topsoil", select=c(Site, Clay, Silt, Sand))
TRTw4.comp <- as.data.frame(acomp(TRTw4[fractions], total=100))
TT.plot(tri.data=TRTw4.comp, class.sys="HYPRES.TT", cex=0.6, cex.lab=1, cex.axis=0.8, main="Topsoils dispersed in water & ultra-sonified 4 min", css.names=fractions)

# Integrated log ratio (ilr) transformation -------------------------------

# Define binary partion
cdata <- acomp(ldpsdat[fractions])
bpart <- t(matrix(c(1,-1,-1,
                    0, 1,-1), ncol=3, nrow=2, byrow=T))
CoDaDendrogram(X=acomp(cdata), signary=bpart)					
idata <- as.data.frame(ilr(cdata, V=bpart))
ldpsdat <- cbind(ldpsdat, idata)

# Write cleaned data frame ------------------------------------------------

write.csv(ldpsdat, "LDPSA.dat.csv")

# Example REML analyses ---------------------------------------------------

# Main effects model ilr[Sand|Silt,Clay]
sand.lmer <- lmer(V1~Depth+TRT+(1|Site)+(1|GID), data=ldpsdat)
display(sand.lmer)
plot(V1~fitted(sand.lmer), ldpsdat)

# Alternatively substituting dispersal agent by ultrasonification time interaction for treatments
sand1.lmer <- lmer(V1~Depth+Disp*Ultra+(1|Site)+(1|GID), data=ldpsdat)
display(sand1.lmer)

# Extract and plot Site-level random effects and standard errors
sand.ranef <- ranef(sand.lmer)
sand.se <- se.coef(sand.lmer)
coefplot(sand.ranef$Site[,1], sand.se$Site[,1], varnames=rownames(sand.ranef$Site), xlim=c(-3,3), CI=2, cex.var=0.6, cex.pts=0.9, main="")

# Main effects model ilr[Silt|Clay]
sicl.lmer <- lmer(V2~Depth+TRT+(1|Site)+(1|GID), data=ldpsdat)
display(sicl.lmer)
sicl.ranef <- ranef(sicl.lmer)
sicl.se <- se.coef(sicl.lmer)
coefplot(sicl.ranef$Site[,1], sicl.se$Site[,1], varnames=rownames(sicl.ranef$Site), xlim=c(-1,1), CI=2, cex.var=0.6, cex.pts=0.9, main="")

