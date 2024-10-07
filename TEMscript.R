library(fields)
library(prodlim)
library(dplyr)

# Function that associates to each point "Pixel", the minimal distance to the points in the matrix Points1 (LD) 
# and the coordinates of that point in Points1 and the same for the points in the matrix Points2 (PM), 
# returning a data.frame with the x and y coordinate of the pixel point, the minimal distance to Points1, 
# the minimal distance to Points2, x and y coordinate of the closest point in Points1,
# x and y coordinate of the closest point in Points2
ClosestPixel <- function(Pixel, Points1, Points2) {
  Points1 <- as.matrix(Points1)
  Points2 <- as.matrix(Points2)
  colnames(Points1) <- colnames(Pixel)
  colnames(Points2) <- colnames(Pixel)
  Elements1 <- rbind(Pixel, Points1)
  Elements2 <- rbind(Pixel, Points2)
  DistElements1 <- rdist(Elements1)
  DistElements2 <- rdist(Elements2)
  MinDist1 <- min(DistElements1[1,c(2:dim(DistElements1)[1])])
  MinDist2 <- min(DistElements2[1,c(2:dim(DistElements2)[1])])
  Point1 <- which(DistElements1[1, c(2:dim(DistElements1)[1])] == MinDist1)[1]
  Point2 <- which(DistElements2[1, c(2:dim(DistElements2)[1])] == MinDist2)[1]
  return(cbind(Pixel, MinDist1, MinDist2, t(Points1[Point1,]), t(Points2[Point2,]), Point1, Point2))
} 

# Function to assign the closest grid point to every caveolae
Cav2Grid <- function(x) {
  il <- c(floor((x[1] - floor(min(PixelCoord[,1])))/GridSize)*GridSize + floor(min(PixelCoord[,1])), floor((x[2] - floor(min(PixelCoord[,2])))/GridSize)*GridSize + floor(min(PixelCoord[,2])))
  ir <- c(ceiling((x[1] - floor(min(PixelCoord[,1])))/GridSize)*GridSize + floor(min(PixelCoord[,1])), floor((x[2] - floor(min(PixelCoord[,2])))/GridSize)*GridSize + floor(min(PixelCoord[,2])))
  sl <- c(floor((x[1] - floor(min(PixelCoord[,1])))/GridSize)*GridSize + floor(min(PixelCoord[,1])), ceiling((x[2] - floor(min(PixelCoord[,2])))/GridSize)*GridSize + floor(min(PixelCoord[,2])))
  sr <- c(ceiling((x[1] - floor(min(PixelCoord[,1])))/GridSize)*GridSize + floor(min(PixelCoord[,1])), ceiling((x[2] - floor(min(PixelCoord[,2])))/GridSize)*GridSize + floor(min(PixelCoord[,2])))
  Elements <- rbind(x[c(1,2)], il, ir, sl, sr)
  DistElements <- rdist(Elements)
  DistGrid <- min(DistElements[1,c(2:dim(DistElements)[1])])
  Point <- which(DistElements[1, c(2:dim(DistElements)[1])] == DistGrid)[1] + 1
  print(paste0('Caveolae ', i,' of a total of ', dim(CavData)[1]))
  return(cbind(x[c(1,2)], Elements[Point,], x[3], DistGrid))
}

options("scipen"=100, "digits"=9)
dir <- "/data_lab_MAP/vjimenezj/MABOY"
# dir <- Sys.getenv("MABOY_PATH")
setwd(dir)

# Data load
AllCaveolaeData <- read.csv(file = paste0(dir, "/TEMCAVEOLAE/TEMCAVEOLAE.txt"), header=TRUE, sep=" ")

# Read arguments for file loading
args = commandArgs(trailingOnly=TRUE)

# Parameters for the summary of caveolae number and Pixel information
GridSize <- 20  ##Related to maximum distance from caveolae to grid center

# Pixfile <- args[1]
# LDfile <- args[2]
# PMfile <- args[3]
# PartIDs <- args[4]

Pixfile <- "Fused.CFF2028.2863.PARTE2.tif_CFF2028.2863.PARTE2.tif.1._PART_2_PIXELCOORDINATES_REANALYSIS.txt"
LDfile <- "Fused.CFF2028.2863.PARTE2.tif_CFF2028.2863.PARTE2.tif.1._PART_2_ALLDISTANCESLDPM_REANALYSIS.txt"
PMfile <- "Fused.CFF2028.2863.PARTE2.tif_CFF2028.2863.PARTE2.tif.1._PART_2_ALLDISTANCESPMLD_REANALYSIS.txt"
PartIDs <- "CFF2028.2863.PARTE2.tif.1.__2"

print(args)
print(Sys.time())

# Reading the data of the points in the plama membrane
PMdata <- as.matrix(read.csv(file = paste0(dir,"/",PMfile), header=TRUE, sep="\t"))
# Select the x and y coordinates
PMcoordinates <- PMdata[,c(2,3)]

# Computation of the distances among the PM points
DistancesPM <- rdist(PMcoordinates)
# Selection of the range of points which are closer than twice the distance to the LD plus 50 pixels: EffRangePM
RangePM <- rep(PMdata[,1]*2+100, nrow(PMdata))
EffRangePM <- DistancesPM < RangePM

# Reading the data of the points in the plama membrane
LDdata <- as.matrix(read.csv(file = paste0(dir,"/",LDfile), header=TRUE, sep="\t"))
# Select the x and y coordinates
LDcoordinates <- LDdata[,c(2,3)]

# Computation of the distances among the LD points
DistancesLD <- rdist(LDcoordinates)
# Selection of the range of points which are closer than twice the distance to the LD plus 50 pixels: EffRangePM
RangeLD <- rep(LDdata[,1]*2+100, nrow(LDdata))
EffRangeLD <- DistancesLD < RangeLD

# Reading the pixels inside the two membranes
PixelCoord <- read.csv(file = paste0(dir,"/",Pixfile), header=TRUE, sep="\t")
PixelCoord <- PixelCoord[,c('C1','C2')]

# Creation of a grid to summarize observations
GridX <- seq(from = floor(min(PixelCoord[,1])), to = ceiling(max(PixelCoord[,1])), by = GridSize)
GridY <- seq(from = floor(min(PixelCoord[,2])), to = ceiling(max(PixelCoord[,2])), by = GridSize)
grid <- expand.grid(x = GridX, y = GridY)
grid <- data.frame(grid$x,grid$y)
colnames(PixelCoord) <- colnames(grid)

# Selection of grid elements inside the two membranes
AllCoord <- rbind(PixelCoord, grid)
#sum(duplicated(AllCoord))
#dim(PixelCoord)
InsideGrid <- AllCoord[(duplicated(AllCoord)),]

# Creation of a matrix with the distance from each pixel to the PM
PixelMatrix <- data.frame(matrix (nrow = dim(InsideGrid)[1], ncol = 10))

# # Fast method
# PixelMatrix3[1,] <- ClosestPixel(InsideGrid[1,],PMcoordinates, LDcoordinates)
# for (i in 2:dim(PixelMatrix3)[1]) {
#   if (rdist(InsideGrid[i,], PixelMatrix3[i-1,c(1,2)])<50) {
#     ReducedRangePM <- PMcoordinates[EffRangePM[PixelMatrix3[i-1,9],],]
#     ReducedRangeLD <- LDcoordinates[EffRangeLD[PixelMatrix3[i-1,10],],]
#     PixelMatrix3[i,] <- ClosestPixel(InsideGrid[i,],ReducedRangePM, ReducedRangeLD)
#     h<-0
#     while ( sum(EffRangePM[PixelMatrix3[i-1,9],c(1:(PixelMatrix3[i,9]+h))]) < PixelMatrix3[i,9]) {
#       h <- h+1
#     }
#     PixelMatrix3[i,9] <- h + PixelMatrix3[i,9]
#     h<-0
#     while ( sum(EffRangeLD[PixelMatrix3[i-1,10],c(1:(PixelMatrix3[i,10]+h))]) < PixelMatrix3[i,10]) {
#       h <- h+1
#     }
#     PixelMatrix3[i,10] <- h + PixelMatrix3[i,10]
#   }  else {
#     PixelMatrix3[i,] <- ClosestPixel(InsideGrid[i,],PMcoordinates, LDcoordinates)
#   }
#   print(paste0('Pixel ', i-1,' of a total of ', dim(PixelMatrix3)[1]-1))
# }
# colnames(PixelMatrix3) <- c("GridX", "GridY", "PMdist", "LDdist", "PMx", "PMy", "LDx", "LDy", "PosPM", "PosLD")

# Slow method
PixelMatrix[1,] <- ClosestPixel(InsideGrid[1,],PMcoordinates, LDcoordinates)
for (i in 2:dim(PixelMatrix)[1]) {
  PixelMatrix[i,] <- ClosestPixel(InsideGrid[i,],PMcoordinates, LDcoordinates)
  print(paste0('Pixel ', i-1,' of a total of ', dim(PixelMatrix)[1]-1))
}
colnames(PixelMatrix) <- c("GridX", "GridY", "PMdist", "LDdist", "PMx", "PMy", "LDx", "LDy", "PosPM", "PosLD")

# Selecting just data from caveolae of the cell under analysis
CavData <- select(filter(AllCaveolaeData, PartID == PartIDs), c('X', 'Y', 'Type'))

# Assignment of the caveolae to the closest grid point
Caveodata <- data.frame(matrix (nrow = dim(CavData)[1], ncol = 6))
for (i in 1:dim(CavData)[1]) {
  Caveodata[i,] <- Cav2Grid(CavData[i,])
}
colnames(Caveodata) <- c("CavX", "CavY", "GridX", "GridY", 'Type', "DistGrid")

# Binding all data in two final datasets
PixelMatrix[,"Totaldist"] <- PixelMatrix[,"PMdist"] + PixelMatrix[,"LDdist"]
CavFinal <- merge(Caveodata, PixelMatrix, by = c("GridX", "GridY"))

PixelMatrix[,"PartID"] <- PartIDs
CavFinal[,"PartID"] <- PartIDs

saveRDS(PixelMatrix, file = paste0(dir,'/NewAnalysisSlow/',PartIDs,'_PixelFile.rds'))
saveRDS(CavFinal, file = paste0(dir,'/NewAnalysisSlow/',PartIDs,'_CavFile.rds'))

print(Sys.time())
