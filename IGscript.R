library(fields)
library(prodlim)
library(dplyr)

# Function that associates to each point "Pixel", the minimal distance to 
# the points in the matrix Points1 (LD) and the coordinates of that point
# in Points1 and the same for the points in the matrix Points2 (PM), 
# returning a data.frame with the x and y coordinate of the pixel point,
# the minimal distance to Points1, # the minimal distance to Points2, 
# x and y coordinate of the closest point in Points1,
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
  return(cbind(Pixel,
               MinDist1,
               MinDist2,
               t(Points1[Point1,]),
               t(Points2[Point2,]), 
               Point1,
               Point2))
} 

# Function to assign the closest grid point to every caveolae
Cav2Grid <- function(x) {
  il <- c(floor((x[1] - floor(min(PixelCoord[,1])))/GridSize)*GridSize + floor(min(PixelCoord[,1])), 
          floor((x[2] - floor(min(PixelCoord[,2])))/GridSize)*GridSize + floor(min(PixelCoord[,2])))
  ir <- c(ceiling((x[1] - floor(min(PixelCoord[,1])))/GridSize)*GridSize + floor(min(PixelCoord[,1])),
          floor((x[2] - floor(min(PixelCoord[,2])))/GridSize)*GridSize + floor(min(PixelCoord[,2])))
  sl <- c(floor((x[1] - floor(min(PixelCoord[,1])))/GridSize)*GridSize + floor(min(PixelCoord[,1])), 
          ceiling((x[2] - floor(min(PixelCoord[,2])))/GridSize)*GridSize + floor(min(PixelCoord[,2])))
  sr <- c(ceiling((x[1] - floor(min(PixelCoord[,1])))/GridSize)*GridSize + floor(min(PixelCoord[,1])), 
          ceiling((x[2] - floor(min(PixelCoord[,2])))/GridSize)*GridSize + floor(min(PixelCoord[,2])))
  Elements <- rbind(x[c(1,2)], il, ir, sl, sr)
  DistElements <- rdist(Elements)
  DistGrid <- min(DistElements[1,c(2:dim(DistElements)[1])])
  Point <- which(DistElements[1, c(2:dim(DistElements)[1])] == DistGrid)[1] + 1
  print(paste0('Caveolae ', i,' of a total of ', dim(CavData)[1]))
  return(cbind(x[c(1,2)], Elements[Point,], x[3], DistGrid))
}

options(scipen=999, digits=9)
dir <- "/data_lab_MAP/vjimenezj/MABOY"
# dir <- Sys.getenv("MABOY_PATH")
setwd(dir)

# Data load
AllCaveolaeData <- read.csv(file = paste0(dir, "/DataIG/IGTEMCAVEOLAE/igtemcaveolae.txt"), header=TRUE, sep=" ")

# Read arguments for file loading
args = commandArgs(trailingOnly=TRUE)

# Parameters for the summary of caveolae number and Pixel information
GridSize <- 20  ##Related to maximum distance from caveolae to grid center

Pixfile <- args[1]
LDfile <- args[2]
PMfile <- args[3]
PartIDs <- args[4]

# Pixfile <- "I01_02082019_HCO3321_F01_P05.tif_C13_PIXELCOORDINATES_REANALYSIS.txt"
# LDfile <- "I01_02082019_HCO3321_F01_P05.tif_C13_ALLDISTANCESLDPM_REANALYSIS.txt"
# PMfile <- "I01_02082019_HCO3321_F01_P05.tif_C13_ALLDISTANCESPMLD_REANALYSIS.txt"
# PartIDs <- "I01_02082019_HCO3321_F01_P05.tif_C13"

print(args)
print(Sys.time())

# Reading the data of the points in the plama membrane
PMdata <- as.matrix(read.csv(file = paste0(dir,"/DataIG/",PMfile), header=TRUE, sep="\t"))
# Select the x and y coordinates
PMcoordinates <- PMdata[,c(2,3)]

# Computation of the distances among the PM points
DistancesPM <- rdist(PMcoordinates)
# Selection of the range of points which are closer than twice the distance to the LD plus 50 pixels: EffRangePM
RangePM <- rep(PMdata[,1]*2+100, nrow(PMdata))
EffRangePM <- DistancesPM < RangePM

# Reading the data of the points in the plama membrane
LDdata <- as.matrix(read.csv(file = paste0(dir,"/DataIG/",LDfile), header=TRUE, sep="\t"))
# Select the x and y coordinates
LDcoordinates <- LDdata[,c(2,3)]

# Computation of the distances among the LD points
DistancesLD <- rdist(LDcoordinates)
# Selection of the range of points which are closer than twice the distance to the LD plus 50 pixels: EffRangePM
RangeLD <- rep(LDdata[,1]*2+100, nrow(LDdata))
EffRangeLD <- DistancesLD < RangeLD

# Reading the pixels inside the two membranes
PixelCoord <- read.csv(file = paste0(dir,"/DataIG/",Pixfile), header=TRUE, sep="\t")
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

# Final Pixel dataset
PixelMatrix[,"Totaldist"] <- PixelMatrix[,"PMdist"] + PixelMatrix[,"LDdist"]
PixelMatrix[,"PartID"] <- PartIDs

# Final dataset of caveolae inside the cell
CavInsideOfCell <- merge(Caveodata, PixelMatrix, by = c("GridX", "GridY"))
CavInsideOfCell[,"PartID"] <- PartIDs
CavInsideOfCell[,"Inside"] <- 1

# Final dataset of caveolae outside the cell
CavWithoutPixel <- anti_join(Caveodata, CavInsideOfCell)

if (dim(Caveodata)[1] - dim(CavInsideOfCell)[1] - dim(CavWithoutPixel)[1]) {
  print("Error seleccionando caveolas dentro y fuera")
} else {print("Caveolas procesadas correctamente")}

CavOutOfCell <- data.frame(matrix (nrow = dim(CavWithoutPixel)[1], ncol = 10))
CavOutOfCell[1,] <- ClosestPixel(CavWithoutPixel[1, c(1,2)],PMcoordinates, LDcoordinates)
for (i in 2:dim(CavOutOfCell)[1]) {
  CavOutOfCell[i,] <- ClosestPixel(CavWithoutPixel[i, c(1,2)],PMcoordinates, LDcoordinates)
  print(paste0('Pixel ', i-1,' of a total of ', dim(CavOutOfCell)[1]-1))
}
colnames(CavOutOfCell) <- c("CavX", "CavY", "PMdist", "LDdist", "PMx", "PMy", "LDx", "LDy", "PosPM", "PosLD")
CavOutOfCell$Type <- CavWithoutPixel$Type
CavOutOfCell$Inside <- 0
CavOutOfCellUnchanged <- CavOutOfCell

for (i in 1:dim(CavOutOfCell)[1]) {
  if (CavOutOfCell[i, "Type"] == 0) {
    CavOutOfCell[i, "GridX"] <- NA
    CavOutOfCell[i, "GridY"] <- NA
    CavOutOfCell[i, "DistGrid"] <- NA
    CavOutOfCell[i, "Totaldist"] <- NA
    CavOutOfCell[i, "Inside"] <- NA
    
  } else {
    if (CavOutOfCell[i, "PMdist"] > CavOutOfCell[i, "LDdist"]) {
      CavOutOfCell[i, "GridX"] <- CavOutOfCell[i, "LDx"]
      CavOutOfCell[i, "GridY"] <- CavOutOfCell[i, "LDy"]
      CavOutOfCell[i, "DistGrid"] <- CavOutOfCell[i, "LDdist"]
      CavOutOfCell[i, "LDdist"] <- 0
      CavOutOfCell[i, "PMdist"] <- LDdata[CavOutOfCell[i, "PosLD"], 1]
      CavOutOfCell[i, "PMx"] <- LDdata[CavOutOfCell[i, "PosLD"], "closerPMX"]
      CavOutOfCell[i, "PMy"] <- LDdata[CavOutOfCell[i, "PosLD"], "closerPMY"]
      CavOutOfCell[i, "Totaldist"] <- LDdata[CavOutOfCell[i, "PosLD"], "results"]
    } else {
      CavOutOfCell[i, "GridX"] <- CavOutOfCell[i, "PMx"]
      CavOutOfCell[i, "GridY"] <- CavOutOfCell[i, "PMy"]
      CavOutOfCell[i, "DistGrid"] <- CavOutOfCell[i, "PMdist"]
      CavOutOfCell[i, "PMdist"] <- 0
      CavOutOfCell[i, "LDdist"] <- PMdata[CavOutOfCell[i, "PosPM"], 1]
      CavOutOfCell[i, "LDx"] <- PMdata[CavOutOfCell[i, "PosPM"], "closerLDX"]
      CavOutOfCell[i, "LDy"] <- PMdata[CavOutOfCell[i, "PosPM"], "closerLDY"]
      CavOutOfCell[i, "Totaldist"] <- PMdata[CavOutOfCell[i, "PosPM"], "results"]
      
    }
    if (CavOutOfCell[i, "DistGrid"] > 100) {CavOutOfCell[i, "Inside"] <- NA}
  }
}
CavOutOfCell[,"PartID"] <- PartIDs


# Creating the final caveolae dataset by binding caveolae out and inside the cell
CavFinal <- rbind(CavInsideOfCell, CavOutOfCell)


saveRDS(PixelMatrix, file = paste0(dir,'/ResultsIG/',PartIDs,'_PixelFile.rds'))
saveRDS(CavFinal, file = paste0(dir,'/ResultsIG/',PartIDs,'_CavFile.rds'))

print(Sys.time())


