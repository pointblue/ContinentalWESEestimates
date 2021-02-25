# TODO: Add comment
# 
# Author: lsalas
###############################################################################


## This is a quick script to change the definition of fast ice predictability
## It is now defined as it being there in early October the past 5 years.

load("c:/users/lsalas/downloads/studyarea_points_wNearLandPenguins_old.RData")
libs <- c("rgdal", "proj4","rgeos","maptools","raster","stringr","plyr","dplyr","xml2","httr","data.table","SDraw")
lapply(libs, require, character.only = TRUE)
pathToGit<-"c:/users/lsalas/git/fasticecovars/data/"
source("c:/users/lsalas/git/fasticecovars/scripts/fastIceCovars_utils.R")
spdf<-studyarea_pointswLandPenguins[,"pointid"]
iceyear<-2011
primaryproj<-CRS(projection(studyarea_pointswLandPenguins))
icemonth="Nov";iceyear=2011
filename<-getNICfilename(getmonth=icemonth,getyear=iceyear)	#   ********************* User specifies date + point buffer and the function does the rest
nicdtdf<-data.frame(NICdate=1:NROW(filename),FileName=filename);print(nicdtdf)

print("******************************************************************************************************************")
print("***** STOP: review the NIC dates printed above; choose one date by entering the line number below  ***************")
print("******************************************************************************************************************")

fn<-1 #  ******************************************** User specifies the desired available NIC date

#nf<-readline(paste0("Which NICdate to use? (1 to ",NROW(filename),"): "))
#Downloading one file using the date selected, and unzipping.
if(fn<1 | fn>NROW(filename)){fn<-1}

nicsavedir<-"c:/temp/"
savename<-nicDownload(filename[fn],nicsavedir=nicsavedir)	

myareas<-getFastIce(fileloc=savename,dataproj=primaryproj,nicsavedir=nicsavedir)

sppdata<-getFastIceSPPdata(iceyear=iceyear,iceareas=myareas$fast,spdf=spdf,keyfield="pointid",nicsavedir=nicsavedir)

names(sppdata)<-gsub("recId","pointid",names(sppdata))
head(sppdata)


sppdata<-sppdata[,c("pointid","PredictabilityOct5Years")]

studyarea_pointswLandPenguins<-merge(studyarea_pointswLandPenguins,sppdata,by="pointid",all.x=T)
studyarea_pointswLandPenguins$PredictabilityOct5Years<-ifelse(is.na(studyarea_pointswLandPenguins$PredictabilityOct5Years),0,studyarea_pointswLandPenguins$PredictabilityOct5Years)
save(studyarea_pointswLandPenguins,file="c:/users/lsalas/downloads/studyarea_points_wNearLandPenguins.RData")

load("c:/users/lsalas/downloads/continentalWESE_old.RData")
stapdf<-as.data.frame(studyarea_pointswLandPenguins)
stapdf<-stapdf[,c("pointid","PredictabilityOct5Years")]
wesedf<-merge(wesedf,stapdf,by.x="gridCellId",by.y="pointid",all.x=TRUE)
save(wesedf,file="c:/users/lsalas/downloads/continentalWESE.RData")

