library(pracma)
library(grImport)
#library(grConvert)
library(plyr)
library(surv2sampleComp)
source("extractFromKMcurveOrig.R")
source("funktioner.R")
library(miscTools)
setwd("H:/Teori/Rituximab")


library(XML)

# figur 2C - svg elementer
# id red: g3897 
# id blue: g3893
# id xaxis: g3525
# id yaxis: g3529
# id 48 mdr: g3625
doc <- htmlParse("Extract from pdf/2C48.svg")
# Extract the coordinates, as strings
path <- getNodeSet(doc, "//svg/g/path")
g <- getNodeSet(doc, "//svg/g/g")
red <- getTransformCoordinates("g3897", g)
blue <- getTransformCoordinates("g3893", g)
xaxis <- getTransformCoordinates("g3525", g)
yaxis <- getTransformCoordinates("g3529", g)
mth48 <- getTransformCoordinates("g3625", g)
red_res <- transformCoordinates(svgStrToCoordinates(red[[2]]), red[[1]])
blue_res <- transformCoordinates(svgStrToCoordinates(blue[[2]]), blue[[1]])
xaxis_res <- transformCoordinates(svgStrToCoordinates(xaxis[[2]]), xaxis[[1]])
yaxis_res <- transformCoordinates(svgStrToCoordinates(yaxis[[2]]), yaxis[[1]])
mth48_res <- transformCoordinates(svgStrToCoordinates(mth48[[2]]), mth48[[1]])
# replace x-coordinate for 48 months in last point of curve
red_res[1,length(red_res[1,])] <- mth48_res[1,1]
blue_res[1,length(blue_res[1,])] <- mth48_res[1,1]
# "move" graph to (0,0) by subtracting coordinates for first point in curve
moveX <- red_res[1,1]
moveY <- red_res[2,1]
red_res <- moveTo0(red_res, moveX, moveY)
blue_res <- moveTo0(blue_res, moveX, moveY)
xaxis_res <- moveTo0(xaxis_res, moveX, moveY)
yaxis_res <- moveTo0(yaxis_res, moveX, moveY)
mth48_res <- moveTo0(mth48_res, moveX, moveY)
monthsGained <- calculateMonthsGained(blue_res, red_res, xaxis_res, yaxis_res, mth48_res)
convert_red <- convertToKMcoord(red_res,xaxis_res,mth48_res)
file <- "2CRituximab"
nrisk_red <- c(505,484,459,444,428,325,220,93,19)
createInputFiles(convert_red, nrisk_red, file)
convert_blue <- convertToKMcoord(blue_res,xaxis_res,mth48_res)
nrisk_blue <- c(513,492,460,425,393,302,188,75,20)
file <- "2CControl"
createInputFiles(convert_blue, nrisk_blue, file)
digisurvfile<- "2CControl_coord.txt" #Input survival times from graph reading
nriskfile<-"2CControl_int.txt" #Input reported number at risk
KMdatafile<-"2CControl_events_cens.txt" #Output file events and cens
KMdataIPDfile<-"2CControl_output_IPD.txt" #Output file for IPD
tot.events<-129 #tot.events = total no. of events reported. If not reported, then tot.events="NA"
arm.id<-1 #arm indicator
dataFromKMCurve(digisurvfile, nriskfile, KMdatafile, KMdataIPDfile, tot.events, arm.id)
digisurvfile<-"2CRituximab_coord.txt" #Input survival times from graph reading
nriskfile<-"2CRituximab_int.txt" #Input reported number at risk
KMdatafile<-"2CRituximab_events_cens.txt" #Output file events and cens
KMdataIPDfile<-"2CRituximab_output_IPD.txt" #Output file for IPD
tot.events<-80 #tot.events = total no. of events reported. If not reported, then tot.events="NA"
arm.id<-1 #arm indicatordataFromKMCurve(digisurvfile, nriskfile, KMdatafile, KMdataIPDfile, tot.events, arm.id)
dataFromKMCurve(digisurvfile, nriskfile, KMdatafile, KMdataIPDfile, tot.events, arm.id)



# bootstrap CI
res <- bootAreal(convert_red, convert_blue, 1000)
quantile(res[,3], c(.025, .975))
mean(res[,3])



# figur 2B - svg elementer
# id red: g3901 
# id blue: g3499
# id xaxis: g3371
# id yaxis: g3375
# id 48 mdr: g3455
doc <- htmlParse("Extract from pdf/Rituxima maintenence figure 2B.svg")
# Extract the coordinates, as strings
path <- getNodeSet(doc, "//svg/g/path")
g <- getNodeSet(doc, "//svg/g/g")
red <- getTransformCoordinates("g3901", g)
blue <- getTransformCoordinates("g3499", g)
xaxis <- getTransformCoordinates("g3371", g)
yaxis <- getTransformCoordinates("g3375", g)
mth48 <- getTransformCoordinates("g3455", g)
red_res <- transformCoordinates(svgStrToCoordinates(red[[2]]), red[[1]])
blue_res <- transformCoordinates(svgStrToCoordinates(blue[[2]]), blue[[1]])
xaxis_res <- transformCoordinates(svgStrToCoordinates(xaxis[[2]]), xaxis[[1]])
yaxis_res <- transformCoordinates(svgStrToCoordinates(yaxis[[2]]), yaxis[[1]])
mth48_res <- transformCoordinates(svgStrToCoordinates(mth48[[2]]), mth48[[1]])
# replace x-coordinate for 48 months in last point of curve
red_res[1,length(red_res[1,])] <- mth48_res[1,1]
blue_res[1,length(blue_res[1,])] <- mth48_res[1,1]
# "move" graph to (0,0) by subtracting coordinates for first point in curve
moveX <- red_res[1,1]
moveY <- red_res[2,1]
red_res <- moveTo0(red_res, moveX, moveY)
blue_res <- moveTo0(blue_res, moveX, moveY)
xaxis_res <- moveTo0(xaxis_res, moveX, moveY)
yaxis_res <- moveTo0(yaxis_res, moveX, moveY)
mth48_res <- moveTo0(mth48_res, moveX, moveY)
monthsGained <- calculateMonthsGained(blue_res, red_res, xaxis_res, yaxis_res, mth48_res)
convert_red <- convertToKMcoord(red_res,xaxis_res,mth48_res)
file <- "2BRituximab"
nrisk_red <- c(505,483,455,441,414,312,209,91,17)
createInputFiles(convert_red, nrisk_red, file)
convert_blue <- convertToKMcoord(blue_res,xaxis_res,mth48_res)
nrisk_blue <- c(513,487,452,417,380,286,170,71,18)
file <- "2BControl"
createInputFiles(convert_blue, nrisk_blue, file)
digisurvfile<-"2BControl_coord.txt" #Input survival times from graph reading
nriskfile<-"2BControl_int.txt" #Input reported number at risk
KMdatafile<-"2BControl_events_cens.txt" #Output file events and cens
KMdataIPDfile<-"2BControl_output_IPD.txt" #Output file for IPD
tot.events<-167 #tot.events = total no. of events reported. If not reported, then tot.events="NA"
arm.id<-1 #arm indicator
dataFromKMCurve(digisurvfile, nriskfile, KMdatafile, KMdataIPDfile, tot.events, arm.id)
digisurvfile<-"2BRituximab_coord.txt" #Input survival times from graph reading
nriskfile<-"2BRituximab_int.txt" #Input reported number at risk
KMdatafile<-"2BRituximab_events_cens.txt" #Output file events and cens
KMdataIPDfile<-"2BRituximab_output_IPD.txt" #Output file for IPD
tot.events<-102 #tot.events = total no. of events reported. If not reported, then tot.events="NA"
arm.id<-1 #arm indicator
dataFromKMCurve(digisurvfile, nriskfile, KMdatafile, KMdataIPDfile, tot.events, arm.id)

df <- data.frame(t(convert_blue))
df <- df[order(df[,1],df[,2]),]
df <- df[!duplicated(df$X1),]

# figur 2A - svg elementer
# id red: g3885 
# id blue: g3177
# id xaxis: g3191
# id yaxis: g3195
# id 48 mdr: g3275
doc <- htmlParse("Extract from pdf/Rituxima maintenence figure 2A.svg")
# Extract the coordinates, as strings
path <- getNodeSet(doc, "//svg/g/path")
g <- getNodeSet(doc, "//svg/g/g")
red <- getTransformCoordinates("g3885", g)
blue <- getTransformCoordinates("g3177", g)
xaxis <- getTransformCoordinates("g3191", g)
yaxis <- getTransformCoordinates("g3195", g)
mth48 <- getTransformCoordinates("g3275", g)
red_res <- transformCoordinates(svgStrToCoordinates(red[[2]]), red[[1]])
blue_res <- transformCoordinates(svgStrToCoordinates(blue[[2]]), blue[[1]])
xaxis_res <- transformCoordinates(svgStrToCoordinates(xaxis[[2]]), xaxis[[1]])
yaxis_res <- transformCoordinates(svgStrToCoordinates(yaxis[[2]]), yaxis[[1]])
mth48_res <- transformCoordinates(svgStrToCoordinates(mth48[[2]]), mth48[[1]])
# replace x-coordinate for 48 months in last point of curve
# for den blå kurve er det 2. sidste koordinat
red_res[1,length(red_res[1,])] <- mth48_res[1,1]
blue_res[1,122] <- mth48_res[1,1]
blue_res <- blue_res[,-(123:126)] 
# "move" graph to (0,0) by subtracting coordinates for first point in curve
moveX <- red_res[1,1]
moveY <- red_res[2,1]
red_res <- moveTo0(red_res, moveX, moveY)
blue_res <- moveTo0(blue_res, moveX, moveY)
xaxis_res <- moveTo0(xaxis_res, moveX, moveY)
yaxis_res <- moveTo0(yaxis_res, moveX, moveX)
mth48_res <- moveTo0(mth48_res, moveX, moveX)
monthsGained <- calculateMonthsGained(blue_res, red_res, xaxis_res, yaxis_res, mth48_res)
convert_red <- convertToKMcoord(red_res,xaxis_res,mth48_res)
file <- "2ARituximab"
nrisk_red <- c(505,472,445,423,404,307,207,84,17)
createInputFiles(convert_red, nrisk_red, file)
convert_blue <- convertToKMcoord(blue_res,xaxis_res,mth48_res)
# af en eller anden grund begynder kurven ikke i (1,0)...?
convert_blue <- cbind(c(0,1), convert_blue)
nrisk_blue <- c(513,469,415,367,334,247,161,70,16)
file <- "2AControl"
createInputFiles(convert_blue, nrisk_blue, file)
digisurvfile<-"2AControl_coord.txt" #Input survival times from graph reading
nriskfile<-"2AControl_int.txt" #Input reported number at risk
KMdatafile<-"2AControl_events_cens.txt" #Output file events and cens
KMdataIPDfile<-"2AControl_output_IPD.txt" #Output file for IPD
tot.events<-218 #tot.events = total no. of events reported. If not reported, then tot.events="NA"
arm.id<-1 #arm indicator
dataFromKMCurve(digisurvfile, nriskfile, KMdatafile, KMdataIPDfile, tot.events, arm.id)
digisurvfile<-"2ARituximab_coord.txt" #Input survival times from graph reading
nriskfile<-"2ARituximab_int.txt" #Input reported number at risk
KMdatafile<-"2ARituximab_events_cens.txt" #Output file events and cens
KMdataIPDfile<-"2ARituximab_output_IPD.txt" #Output file for IPD
tot.events<-130 #tot.events = total no. of events reported. If not reported, then tot.events="NA"
arm.id<-1 #arm indicator
dataFromKMCurve(digisurvfile, nriskfile, KMdatafile, KMdataIPDfile, tot.events, arm.id)

# bootstrap CI
res <- bootAreal(convert_red, convert_blue, 1000)
quantile(res[,3], c(.025, .975))
mean(res[,3])

######################################################################################
# beregning af konfidensintervaller ud fra rekonstruerede data
######################################################################################
D=pbc.sample()
surv2sample(D$time, D$status, D$group, npert=500, timepoints=c(2,4,6,8),
           quanprobs =c(0.25, 0.5), tau=8, procedure="KM")
# fig 2C
ritData = read.table("2CRituximab_output_IPD.txt")
ritData$Group <- 1
controlData = read.table("2CControl_output_IPD.txt")
controlData$Group <- 0
D <- rbind(ritData, controlData)
View(allData)
surv2sample(D$V1, D$V2, D$Group, npert=500, timepoints=c(12,24,36,48),
            quanprobs =c(0.25, 0.5), tau=48, procedure="KM")
# fig 2B
ritData = read.table("2BRituximab_output_IPD.txt")
ritData$Group <- 1
controlData = read.table("2BControl_output_IPD.txt")
controlData$Group <- 0
D <- rbind(ritData, controlData)
View(allData)
surv2sample(D$V1, D$V2, D$Group, npert=500, timepoints=c(12,24,36,48),
            quanprobs =c(0.25, 0.5), tau=48, procedure="KM")

# fig 2A
ritData = read.table("2ARituximab_output_IPD.txt")
ritData$Group <- 1
controlData = read.table("2AControl_output_IPD.txt")
controlData$Group <- 0
D <- rbind(ritData, controlData)
View(allData)
surv2sample(D$V1, D$V2, D$Group, npert=500, timepoints=c(12,24,36,48),
            quanprobs =c(0.25, 0.5), tau=48, procedure="KM")




# # plot af de transformerede koordinater
plot(red_res[1,], red_res[2,], col="white", xlim=c(300,500), ylim=c(600,800),pch=".")
lines(red_res[1,], red_res[2,], col="blue", xlim=c(300,500), ylim=c(600,800))
par(new=TRUE)
plot(blue_res[1,], blue_res[2,], col="white", xlim=c(300,500), ylim=c(600,800), pch=".")
lines(blue_res[1,], blue_res[2,], col="red", xlim=c(300,500), ylim=c(600,800))
par(new=TRUE)
plot(xaxis_res[1,], xaxis_res[2,], col="black", xlim=c(300,500), ylim=c(600,800))
lines(xaxis_res[1,], xaxis_res[2,], col="green", xlim=c(300,500), ylim=c(600,800))
par(new=TRUE)
plot(yaxis_res[1,], yaxis_res[2,], col="red", xlim=c(300,500), ylim=c(600,800))
lines(yaxis_res[1,], yaxis_res[2,], col="green", xlim=c(300,500), ylim=c(600,800))
par(new=TRUE)
plot(mth48_res[1,], mth48_res[2,], col="blue", xlim=c(300,500), ylim=c(600,800))
lines(mth48_res[1,], mth48_res[2,], col="green", xlim=c(300,500), ylim=c(600,800))

# plot af konverterede KM kurver
plot(convert_red[1,], convert_red[2,], col="white", xlim=c(0,50), ylim=c(0,1),pch=".")
lines(convert_red[1,], convert_red[2,], col="red", xlim=c(0,50), ylim=c(0,1))
par(new=TRUE)
plot(convert_blue[1,], convert_blue[2,], col="white", xlim=c(0,50), ylim=c(0,1), pch=".")
lines(convert_blue[1,], convert_blue[2,], col="blue", xlim=c(0,50), ylim=c(0,1))
par(new=TRUE)
plot(xaxis_res[1,], xaxis_res[2,], col="black", xlim=c(0,50), ylim=c(0,1))
lines(xaxis_res[1,], xaxis_res[2,], col="green", xlim=c(0,50), ylim=c(0,1))
par(new=TRUE)
plot(yaxis_res[1,], yaxis_res[2,], col="red", xlim=c(0,50), ylim=c(0,1))
lines(yaxis_res[1,], yaxis_res[2,], col="green", xlim=c(0,50), ylim=c(0,1))
