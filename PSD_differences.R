# run the LDPSA_starter.R first
# M. Walsh, October 2014

# Data setup --------------------------------------------------------------

# Ultrasonic treatment differences between samples dispersed in water
w1  <- subset(ldps.comp, TRT=="w1", select=c(SSN, V1, V2))
w4  <- subset(ldps.comp, TRT=="w4", select=c(SSN, V1, V2))
colnames(w1) <- c("SSN", "w1v1", "w1v2")
colnames(w4) <- c("SSN", "w4v1", "w4v2")
wsa <- merge(w1, w4, by="SSN")
attach(wsa)
wsa$dwsa <- (w4v1-w1v1)+(w4v2-w1v2)
detach(wsa)

# Ultrasonic treatment differences between samples dispersed in Calgon
c1  <- subset(ldps.comp, TRT=="c1", select=c(SSN, V1, V2))
c4  <- subset(ldps.comp, TRT=="c4", select=c(SSN, V1, V2))
colnames(c1) <- c("SSN", "c1v1", "c1v2")
colnames(c4) <- c("SSN", "c4v1", "c4v2")
csa <- merge(c1, c4, by="SSN")
attach(csa)
csa$dcsa <- (c4v1-c1v1)+(c4v2-c1v2)
detach(csa)





