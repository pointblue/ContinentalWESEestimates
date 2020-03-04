# TODO: Add comment
# 
# Author: lsalas
###############################################################################


libs<-c("ggplot2","plyr")
lapply(libs, require, character.only = TRUE)
pathToLocalGit<-"C:/Users/lsalas/git/ContinentalWESEestimates/data/"

######################  FUNCTIONS WE'LL NEED
source("C:/Users/lsalas/git/ContinentalWESEestimates/scripts/countSealsFromTags_functions.R")

## Load the data - from script countSealsFromTags.R
load(file=paste0(pathToLocalGit,"estimatesByMap_unadjusted.RData"))

## Adjusting the estimates for differences with ground counts is done with call to this function:
adjRates<-predictDetRates(pathToGit=pathToLocalGit,dat=df,keyFieldName="regionMapId",islandWeight=0.98)
## Here we are saying that 90% of the colonies out there are island colonies
## Two methods: assuming all colonies are line those in Erebus Bay, weighing the average by the proportion of island/mainland colonies, or...
## simply using a model that corrects for island or mainland, using the Erebus data as refrence

## Now we calculate the numbers
countdf<-merge(df,adjRates,by="regionMapId")

## Since detectionRate = count-in-map/count-on-the-ground, and we want count-on-the-ground, then we divide the count-in-map/detectionRate
countdf$mdlColEstimate<-round(countdf$estimateH/countdf$wgtPredColRate)
countdf$mdlColUpper<-round(countdf$upperH/countdf$wgtPredColRate)
countdf$mdlColLower<-round(countdf$lowerH/countdf$wgtPredColRate)

countdf$mdlIslEstimate<-round(countdf$estimateH/countdf$wgtPredIslRate)
countdf$mdlIslUpper<-round(countdf$upperH/countdf$wgtPredColRate)
countdf$mdlIslLower<-round(countdf$lowerH/countdf$wgtPredColRate)


estByRegionCol<-as.data.frame(countdf[,c("region","mdlColLower","mdlColEstimate","mdlColUpper")] %>% group_by(region) %>% dplyr::summarize(Lower=round(sum(mdlColLower)),Estimate=round(sum(mdlColEstimate)),Upper=round(sum(mdlColUpper))))
estByRegionCol<-rbind(estByRegionCol,data.frame(region="Total",Lower=round(sum(countdf$mdlColLower)),Estimate=round(sum(countdf$mdlColEstimate)),Upper=round(sum(countdf$mdlColUpper))))
print(estByRegionCol)

save(countdf,file=paste0(pathTogit,"data/correctedCounts.RData"))
## Use 95% and Colony model
## Prepare data for regression model
## Setup the notebook

