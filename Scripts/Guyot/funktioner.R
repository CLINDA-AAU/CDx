##################################################################################
# Funktioner relateret til på basis af svg fil med to KM-kurver at beregne arealet
# mellem de to kurver
##################################################################################
calculateMonthsGained <- function(blue, red, xaxis, yaxis, mth48){
  areaControl <- trapz(blue[1,], blue[2,])
  areaRituximab <- trapz(red[1,], red[2,])
  areaDiff <- areaControl - areaRituximab
  rectangleX <- xaxis[2,2] 
  rectangleY <- mth48[1,1]
  areaTotal <- rectangleX * rectangleY
  monthsGained <- (areaDiff / areaTotal) * 48
  return(monthsGained)
}

moveTo0 <- function(object, moveX, moveY){
  object[1,] <- object[1,] - moveX
  object[2,] <- object[2,] - moveY
  return(object)
}

getTransformCoordinates <- function(id, g){
  for(i in 1:length(g)){
    g_id <- xmlGetAttr(g[[i]],"id")
    if(g_id == id){
      transform <- xmlGetAttr(g[[i]],"transform")
      coordinates <- xmlGetAttr(g[[i]][[1]],"d")
    }
  }
  return(c(transform, coordinates) )
} 

transformCoordinates <- function(coordinates, transformation){
  p <- substr(transformation,8,nchar(transformation)-1)
  p <- lapply( strsplit(p, " "), function(u) 
    matrix(as.numeric(unlist(strsplit(u, ","))),ncol=3, byrow=FALSE) )
  p <- data.frame(p)
  p_mat <- matrix(, nrow=nrow(p) , ncol =3)
  p_mat[,1] <- p$X1
  p_mat[,2] <- p$X2
  p_mat[,3] <- p$X3
  p_mat <- rbind(p_mat, c(0,0,1))
  res <- p_mat %*% t(coordinates)
  return(res)
}

svgStrToCoordinates <- function(p){
  if (substr(p,1,1) == "m" ){
    rel <- TRUE
  }else if (substr(p,1,1) == "M" ){
    rel <- FALSE
  }
  p <- substr(p,3,nchar(p))
  p <- lapply( strsplit(p, " "), function(u) 
    matrix(as.numeric(unlist(strsplit(u, ","))),ncol=2,byrow=TRUE) )
  
  p <- data.frame(p)
  p_mat <- matrix(, nrow=nrow(p), ncol =3)
  nrow(p)
  p_mat[,1] <- p$X1
  p_mat[,2] <- p$X2
  p_mat[,3] <- 1
  p_mat
  # absolutte koordinater i stedet for relative
  if (rel) { 
    for (i in 1:length(p_mat[,1])){
      if (i>1){
        p_mat[i,1] <- p_mat[i-1,1] +  p_mat[i,1] 
        p_mat[i,2] <- p_mat[i-1,2] +  p_mat[i,2] 
      }
    }
  }
  return(p_mat)
}


#################################################################################
# Funktion relateret til rekonstruktion af KM data (Guyots metode)
#################################################################################
# konverter koordinater til sandsynligheder og måneder
convertToKMcoord <- function(curve, xaxis, mth48){
  converted <- matrix(0, ncol = length(curve[1,]), nrow=2)
  converted[2,] <- 1 - (curve[2,]/xaxis[2,2])
  converted[1,] <- 48 * curve[1,]/mth48[1,1]
  return(converted)
}

# create input files for algorithm (Guyot)
createInputFiles <- function(coord, nriskList,file){
  df <- data.frame(t(coord))
  df <- df[order(df[,1],df[,2]),]
  df <- df[!duplicated(df$X1),]
  df$Interval <- floor(df$X1 / 6) +1
  if (nrow(df[df$Interval == 0,])!=0){
    df[df$Interval == 0,]$Interval <- 1
  }
  df$k <-seq.int(nrow(df))
  i <-  0
  trisk <- 0
  lower <-  0
  upper <-  0
  nrisk <-  0
  df2 <- data.frame(i, trisk, lower, upper, nrisk)
  for(i in levels(factor(df$Interval))){
    subset <- df[df$Interval == i,]
    lower <- subset[which.min(subset$k),]$k
    upper <- subset[which.max(subset$k),]$k
    newrow = c(i, as.numeric(i)*6-6, lower, upper, nriskList[as.numeric(i)])
    df2 = rbind(df2,newrow)
  }
#   if (max(df2$i) != (max(df$Interval)+1)){
#     df2 = rbind(df2, c(max(df$Interval)+1,48,max(df$k),max(df$k), nriskList[max(df$Interval)+1]))    
#   }
  df2 <- df2[df2$i > 0,]
  df2[, c(1:5)] <- sapply(df2[, c(1:5)], as.numeric)
  #df2$nrisk <- nrisk
  fileName <- paste(file, "_int.txt", sep="")
  write.table(df2, fileName, sep="\t", append = FALSE, row.names = FALSE)
  df <- rename(df, c("X1"="TK", "X2"="Sk"))
  df <- df[c(4,1,2)]
  df[, c(1:3)] <- sapply(df[, c(1:3)], as.numeric)
  fileName <- paste(file, "_coord.txt", sep="")
  write.table(df, fileName, sep="\t", append = FALSE, row.names = FALSE)
}

#################################################################################
# Funktioner til bootstrap af konfidensintervaller
#################################################################################
# funktion, der danner matrix med x-koordinater og spring
spring <- function(KMcoord){
  N <- length(KMcoord[1,])
  n <- (length(KMcoord[1,])-2)/2
  k <- 1
  i <- 1
  spring = matrix(0,nrow=n,ncol=2)
  for(i in 1:(N-1)){
    if(KMcoord[1,i]==KMcoord[1,i+1]){
      spring[k,1] <-  KMcoord[1,i]
      spring[k,2] <- KMcoord[2,i] - KMcoord[2,i+1]
      k <- k + 1
    }
  }
  return(spring)
}

# funktion, der ud fra spring matrix tager en sample med tilbagelægning og beregner areal, som returneres
sampleAreal <- function(spring){
  N <- length(spring[,1])*2+2
  n <- length(spring[,1])
  s <- sample(spring[,1], n, replace = TRUE, prob = spring[,2])
  xCoord <- spring[,1]
  KMcoord = matrix(0, N, ncol=2)
  s <- cbind(xCoord, s)
  seq <- seq(2,N-1,by=2)
  KMcoord[1,] <- c(0,1)
  for(i in seq){
    KMcoord[i,1] <- s[i/2,1]
    KMcoord[i,2] <- KMcoord[i-1,2]
    KMcoord[i+1,1] <- KMcoord[i,1]
    KMcoord[i+1,2] <- KMcoord[i,2] - s[i/2,2]  
  }
  KMcoord[N,1] <- 48
  KMcoord[N,2] <- KMcoord[N-1,2]
  area <- trapz(KMcoord[,1], KMcoord[,2])
  return(area)
}

# funktion, der beregner forskellen i arealer mellem red og blue nRep gange for at kunne bootstrappe CI 
# for forskellen i areal
bootAreal <- function(red, blue, nRep){
  springRed <- spring(red)
  springBlue <- spring(blue) 
  res <- matrix(0, nrow = nRep, ncol = 3)
  for (i in 1:nRep){
    res[i,1] <- sampleAreal(springBlue)
  }
  for (i in 1:1000){
    res[i,2] <- sampleAreal(springRed)
  }
  res[,3] = res[,2]-res[,1]
  return(res)
}
