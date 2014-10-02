# run the LDPSA_starter.R first
# M. Walsh, October 2014

# Load lab covariates -----------------------------------------------------

download("https://www.dropbox.com/s/gcga9uyt2b8cv9h/Lab_cov.csv.zip?dl=0", "Lab_cov.csv.zip", mode="wb")
unzip("Lab_cov.csv.zip", overwrite=T)
labcov <- read.table("Lab_cov.csv", header=T, sep=",")

# Data setup --------------------------------------------------------------

# Ultrasonic treatment differences of samples dispersed in water
w1  <- subset(ldps.comp, TRT=="w1", select=c(Site, PID, SSN, Depth, V1, V2))
w4  <- subset(ldps.comp, TRT=="w4", select=c(SSN, V1, V2))
colnames(w1) <- c("Site", "PID", "SSN", "Depth", "w1v1", "w1v2")
colnames(w4) <- c("SSN", "w4v1", "w4v2")
water <- merge(w1, w4, by="SSN")
attach(water)
water$dv1 <- w4v1-w1v1
water$dv2 <- w4v2-w1v2
water$dw4 <- (w4v1-w1v1)+(w4v2-w1v2)
detach(water)
psdw <- merge(water, labcov, by="SSN")

# Ultrasonic treatment differences of samples dispersed in Calgon
c1  <- subset(ldps.comp, TRT=="c1", select=c(Site, PID, SSN, Depth, V1, V2))
c4  <- subset(ldps.comp, TRT=="c4", select=c(SSN, V1, V2))
colnames(c1) <- c("Site", "PID", "SSN", "Depth", "c1v1", "c1v2")
colnames(c4) <- c("SSN", "c4v1", "c4v2")
calgon <- merge(c1, c4, by="SSN")
attach(calgon)
calgon$dv1 <- c4v1-c1v1
calgon$dv2 <- c4v2-c1v2
calgon$dc4 <- (c4v1-c1v1)+(c4v2-c1v2)
detach(calgon)
psdc <- merge(calgon, labcov, by="SSN")





