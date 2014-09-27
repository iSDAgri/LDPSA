# AfSIS laser diffraction particle size analyses (LDPSA)
# M. Walsh, September 2014

# Set your local working directory here, e.g.
dat_dir <- "/Users/markuswalsh/Documents/LDSF/LDPSA"
setwd(dat_dir)

# Load packages
library(downloader)
library(arm)
library(compositions)

# Load data
download("https://www.dropbox.com/s/je066cib8zt9nh6/LDPSA_60.zip?dl=0", "LDPSA_60.zip", mode="wb")
unzip("LDPSA_60.zip", overwrite=T)
profile <- read.table("Profiles.csv", header=T, sep=",")
samples <- read.table("Samples.csv", header=T, sep=",")
ldpsdat <- merge(profile, samples, by="PID")
ldpsdat <- na.omit(ldpsdat)

# Particle size compositions and isometric log ratio transforms -----------

fractions <- c("Sand","Silt","Clay")
plot(acomp(ldpsdat[fractions]), axes=T, main=" ")
cdata <- acomp(ldpsdat[fractions])

# Define binary partion
bpart <- t(matrix(c(1,-1,-1,
					          0, 1,-1), ncol=3, nrow=2, byrow=T))
CoDaDendrogram(X=acomp(cdata), signary=bpart)					
idata <- as.data.frame(ilr(cdata, V=bpart))
ldpsdat <- cbind(ldpsdat, idata)

 # Example REML analyses ---------------------------------------------------

# Main effects model [Sand|Silt,Clay] ilr
sand.lmer <- lmer(V1~Depth+TRT+(1|Site), data=ldpsdat)
display(sand.lmer)
plot(V1~fitted(sand.lmer), ldpsdat)

# Extract and plot Site-level random effects and standard errors
sand.ranef <- ranef(sand.lmer)
sand.se <- se.coef(sand.lmer)
coefplot(sand.ranef$Site[,1], sand.se$Site[,1], varnames=rownames(sand.ranef$Site), xlim=c(-3,3), CI=2, cex.var=0.6, cex.pts=0.9, main="")

# Main effects model [Silt|Clay] ilr
sicl.lmer <- lmer(V2~Depth+TRT+(1|Site), data=ldpsdat)
display(sicl.lmer)
sicl.ranef <- ranef(sicl.lmer)
sicl.se <- se.coef(sicl.lmer)
coefplot(sicl.ranef$Site[,1], sicl.se$Site[,1], varnames=rownames(sicl.ranef$Site), xlim=c(-2,2), CI=2, cex.var=0.6, cex.pts=0.9, main="")
