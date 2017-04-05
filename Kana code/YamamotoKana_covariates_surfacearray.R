# YamamotoKana_covariates_surfacearray.R
# 
# go through the surface O2 and moisture data for Kana's decomp project
#
# O2 - redox - GHG project
# CS O'Connell, UCB, Silver Lab
#
# see also: Combine-array-sensor-data.R, which points to where the O2 csv file lives


########################################################################
# BRING IN DATA

library(ggplot2)
library(plyr)
library(dplyr)
library(lubridate)

# summarySE using plyr
source("~/Documents/GITHUB/RPersonalFunctionsChristine/summarySE.r")

# where does the O2 csv live
pathdata = "~/Documents/GITHUB/cso044code_HotSpotsHotMoments/HotSpotsHotMomentsAnalysis/HotSpotsHotMoments-Data-Raw/Sensors/SurfaceDataArchive/"
# where to save figures
pathsavefigs = "~/Desktop/RESEARCH PROJECTS/cso043_ArrayDecomp/Kana Code/"

# bring in O2 data
fulldaily <- as.data.frame(read.csv(paste(pathdata,"fulldaily",".csv",sep=""), stringsAsFactors=FALSE))

########################################################################
# QUALITY CONTROL DATA

# make factors where needed
fulldaily$TransectID <- as.factor(fulldaily$TransectID)
fulldaily$TopoLocation <- as.factor(fulldaily$TopoLocation)

# get rid of any O2 values above 0.23
fulldaily$avgO2pct[fulldaily$avgO2pct>22.5] <- NA

# get only the dates when the litterbag experiment was happening
fulldaily$Date3 <- ymd(fulldaily$Date2)
fulldailysubset <- fulldaily[fulldaily$Date3 > ymd(20160601),]


########################################################################
# GRAPH THINGS

# graph to see if this is normal looking
png(file = paste(pathsavefigs, "decomplitterbag_O2.png", sep=""),width=10,height=10,units="in",res=400)
ggplot(fulldailysubset,aes(x=as.Date(Date2),y=avgO2pct,color=TopoLocation)) + 
  geom_point() +
  labs(x="Date",y="O2 concentration")
dev.off()

# graph moisture, same reason
png(file = paste(pathsavefigs, "decomplitterbag_moisture.png", sep=""),width=10,height=10,units="in",res=400)
ggplot(fulldailysubset,aes(x=as.Date(Date2),y=avgVWC,color=TopoLocation)) + 
  geom_point() +
  labs(x="Date",y="Vol Water Content")
dev.off()



