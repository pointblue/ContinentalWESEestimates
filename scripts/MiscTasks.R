# TODO: Edit code to estimate corrected abundances
# 
# Author: lsalas
###############################################################################


###############################################################################
## Summary of the data
getGeneralSummary<-function(){
	################################### DATA
	
	mid<-21758509	## This is Michelle LaRue's tagger Id
	
	## load the data, prepare for analyses
	load("//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/ContinentalEst2/compiledData.R")
	
	print("How many maps inspected?")
	print(NROW(unique(views$mapId)))
	print("By how many taggers?")
	print(NROW(unique(views$taggerId)))
	print("How many tags?")
	print(NROW(unique(tags$tagId)))
	
	crd<-merge(crowd,tags[,c("tagId","mapViewId")],all.x=T)
	crd<-merge(crd,maps[,c("mapId","overlayId")],by="mapId",all.x=T)
	ww<-subset(overlays, overlayId %in% unique(crd$overlayId))
	ww$acDate<-format(ww$acquisition_date,"%Y-%m-%d %H:%M:%S")
	aa<-unique(ww$acDate)
	print("How many acquired satellite images overall? (Here estimated as unique date-times)")
	print(NROW(aa))
	
	print("How many maps had tags?")
	tags<-unique(tags)
	tags<-merge(tags,unique(views[,c("mapViewId","mapId")]),by="mapViewId",all.x=T)
	print(NROW(unique(tags$mapId)))
	
	print("How many maps inspected by ML?")
	print(NROW(unique(subset(views,taggerId==21758509)$mapId)))
	print("In how many did she place tags?")
	print(NROW(unique(subset(tags,taggerId==21758509)$mapId)))
	print("How many tags did she place?")
	print(NROW(unique(subset(tags,taggerId==21758509)$tagId)))
	
	print("How many maps had features?")
	crt<-merge(crowd[,c("tagId","score","agremnt","sensor")],tags[,c("tagId","mapViewId")],by="tagId",all.x=T)
	crt<-subset(crt,!is.na(mapViewId))
	crtm<-merge(crt,views[,c("mapViewId","mapId")],by="mapViewId",all.x=T)
	crtm<-subset(crtm,!is.na(mapId))
	crtm$numFeatures<-1
	print(NROW(unique(crtm$mapId)))
	print("How many tags in these?")
	print(nrow(crtm))
	
	print("How many of these maps did ML inspect?")
	mlMaps<-unique(subset(views,taggerId==21758509)$mapId)
	print(NROW(unique(subset(crtm,mapId %in% mlMaps)$mapId)))
	
	print("In how many did she place tags?")
	crtmMLmaps<-unique(subset(crtm,mapId %in% mlMaps)$mapId)
	print(NROW(unique(subset(tags,(taggerId==21758509) & (mapId %in% crtmMLmaps))$mapId)))
	print("How many tags?")
	print(NROW(unique(subset(tags,(taggerId==21758509) & (mapId %in% crtmMLmaps))$tagId)))
	
	print("How many taggers tagged features?")
	print(NROW(unique(subset(tags,mapId %in% unique(crtm$mapId))$taggerId))-1)
	print("with how many taggers did she overlap?")
	crtmTaggers<-unique(subset(tags,mapId %in% unique(crtm$mapId)))
	crtmTaggers<-subset(crtmTaggers,taggerId!=21758509)
	crtmTaggers<-subset(crtmTaggers,mapId %in% crtmMLmaps)
	print(NROW(unique(crtmTaggers$taggerId)))
	
	
	
}

####  Get the general statistics
getGeneralSummary()


################################################################################
## Question 1:
# Compare filtered count of seals per map for maps inpsected by Michelle, after filtering by CR value
libs<-c("ggplot2","plyr")
lapply(libs, require, character.only = TRUE)

#find the maps that Michelle inspected, and determine how many seals she counted...
mid<-21758509

## load the data, prepare for analyses
load("//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/ContinentalEst2/compiledData.R")
mlMaps<-unique(subset(views,taggerId==mid,select=c("mapId","region")))

#Get the number of tags per map for maps inspected by Michelle
mlTags<-subset(tags,taggerId==mid,select=c("taggerId","mapViewId","region"))
mlCounts<-aggregate(taggerId~mapViewId+region,mlTags,NROW);names(mlCounts)<-c("mapViewId","region","mlCount")
mlCounts<-merge(mlCounts,views[,c("mapViewId","mapId","region")],by=c("mapViewId","region"),all.x=T)
mlMapCounts<-aggregate(mlCount~mapId+region,mlCounts,sum,na.rm=T)
#Add the 0's
mlCountdf<-merge(mlMaps,mlMapCounts,by=c("mapId","region"),all.x=T)
mlCountdf$mlCount<-ifelse(is.na(mlCountdf$mlCount),0,mlCountdf$mlCount)

#Iteratively filter number of features per map by CR and regress against the numbers Michelle counted
crData<-merge(crowd,tags[,c("tagId","mapViewId","region")],by=c("tagId","region"),all.x=T)
crData<-subset(crData,mapId %in% mlCountdf$mapId)

crfvals<-seq(0.7,0.95,0.05)
dat<-data.frame();mdls<-list();mdlres<-data.frame();i<-0
for(crf in crfvals){
	i<-i+1
	tdf<-subset(crData,score >= crf)
	ccdf<-aggregate(sensor~mapId+region,tdf,NROW);names(ccdf)<-c("mapId","region","crCount")
	ccdf$threshold<-paste0("Threshold=",crf)
	ccdf<-merge(ccdf,mlCountdf,by=c("mapId","region"),all.x=T)
	ccdf<-subset(ccdf,!is.na(mlCount))
	dat<-rbind(dat,ccdf)
	mdl<-lm(crCount~mlCount,ccdf)
	mdls[[i]]<-list(mdl=mdl,threshold=crf)
	smdl<-summary(mdl)
	mdlrt<-data.frame(threshold=crf,adjRsq=smdl$adj.r.squared,slope=smdl$coefficients[2,1],intercept=smdl$coefficients[1,1],numMaps=nrow(ccdf),mlCount=sum(ccdf$mlCount),crCount=sum(ccdf$crCount))
	mdlres<-rbind(mdlres,mdlrt)
}

## FIRST RESULT is the table mdlres
## And report that MLR counted 2,854 seals, so at cr threshold of 0.7 we already lose 499 seals.
## TO NOTE: increasing the threshold actually decreases the adj.R square, and decreases both the total number of seals detected and the total number of maps with seals. So, we miss a lot

## SECOND RESULT
p1<-ggplot(data=dat,aes(x=mlCount,y=crCount)) + geom_point(size=1.5) + geom_abline(slope=1,intercept=0,color="black",size=1.2) +
		facet_wrap(~threshold,ncol=3) + theme_bw() + labs(x="Expert count",y="Crowd estimate") +
		theme(axis.text=element_text(size=12),axis.title=element_text(size=14),strip.text=element_text(size=14))
## TO NOTE: even at threshold=0.95 there are maps that may have 15 or more seals and the crowd only counted <2. No visible reduction in over-counting or undercounting.

################################################################################
## Question 2:
# Compare bias-adjusted count of seals per map for maps inpsected by Michelle
load(file="//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/Erebus/shrunkCounts_vsMLR.RData") ## from calcSealNumbers_v4.R

crtvals<-seq(0.7,0.95,0.05)
mdlsh<-list();mdlresh<-data.frame();i<-0
for(crt in crtvals){
	i<-i+1
	tdf<-subset(pdfs_ml,featuresShared==20 & round(crankThreshold,2)==round(crt,2))
	mdl<-lm(estSeals~countML,tdf)
	mdlsh[[i]]<-list(mdl=mdl,threshold=crf)
	smdl<-summary(mdl)
	mdlrt<-data.frame(threshold=crt,adjRsq=smdl$adj.r.squared,slope=smdl$coefficients[2,1],intercept=smdl$coefficients[1,1],
			numMaps=nrow(tdf),mlCount=sum(tdf$countML),crEstimate=sum(tdf$estSeals),meanShrinkage=mean(tdf$meanShrinkage),numMapViews=sum(tdf$numTaggers))
	mdlresh<-rbind(mdlresh,mdlrt)
}

## THIRD RESULT is the table mdlresh
## TO NOTE: increasing the threshold does increase the adj.R square, but as expected decreases both the total number of maps with seals to estimate. 
## Here we don't care about the total number of maps with estimates, because we will use the shrinkage estimate to adjust all maps's counts, including those not searched by ML, or by taggers who overlapped with ML
## BUT we care about representativity, so we want a good number of maps, and a good number of views of these maps (i.e., counting events) from which to obtain the shrinkage.
## At threshold 0.8, the number of features tagged by the surveyors is 4.7 more than the number of seals in the map, on average.

## FOURTH RESULT
pdfs_ml$thresholdLab<-paste("Crowdrank threshold =",pdfs_ml$crankThreshold)
p2<-ggplot(data=subset(pdfs_ml,featuresShared==10 & crankThreshold>0.65),aes(x=countML,y=estSeals)) + 
		geom_point() + geom_errorbar(aes(ymin=estSealslLower,ymax=estSealslUpper)) +
		geom_abline(slope=1,intercept=0,color="black",size=1.2) +
		facet_wrap(~thresholdLab,ncol=3,scales="free") +
		theme_bw() + labs(x="Expert count",y="Crowd estimate") +
		theme(axis.text=element_text(size=12),axis.title=element_text(size=14),strip.text=element_text(size=12))


### Printing out
write.csv(mdlres,file="//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/FinalFigures/Table1.csv")
write.csv(mdlresh,file="//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/FinalFigures/Table2.csv")
jpeg(filename = "//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/FinalFigures/Figure2.jpg",width=2200,height=1500,res=300,quality=100)
print(p1)
dev.off()
jpeg(filename = "//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/MethodsPaper/FinalFigures/Figure3.jpg",width=2500,height=1700,res=300,quality=100)
print(p2)
dev.off()

#####################
## Then switch to file Methods_paper_results2.R

