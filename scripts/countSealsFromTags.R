# TODO: Add comment
# 
# Author: lsalas
###############################################################################


## Load the features
## See in which maps there are features, and which of these were explored by Michelle
## ProbF is the probability of a feature in the collection of maps inspected by 2 or more people - universal
## ProbS is the probability of Michelle finding a seal in a map she inspected - universal
## ProbFS is ss/(ss+ns) - for the maps inspected by ML, is her tag within 3m of a feature? - depending on CR threshold
## ProbSF is ProbFS*ProbS/ProbF - universal
## Corr factor is ProbSF/ProbFS, or ProbS/ProbF, thus independent of threshold unless... 
## We calculate ProbFS by map, and thus ProbSF by map. This means that CorrFactor has expectation ProbS/ProbF. Let's see!

libs<-c("ggplot2","plyr","dplyr","fitdistrplus")
lapply(libs, require, character.only = TRUE)
pathToLocalGit<-"C:/Users/lsalas/git/ContinentalWESEestimates/data/"

######################  FUNCTIONS WE'LL NEED
source("C:/Users/lsalas/git/ContinentalWESEestimates/scripts/countSealsFromTags_functions.R")

################################### DATA
mid<-21758509	## This is Michelle LaRue's tagger Id

## load the data, prepare for analyses
load(paste0(pathToLocalGit,"compiledData.RData"))


##############################
#prepare to count tags 
tags<-unique(tags)
tgvutm<-merge(tags,views[,c("mapViewId","mapId","region")],by=c("mapViewId","region"),all.x=T)
tgvutm$numTags<-1
tgvutm$regionMapId<-paste0(tgvutm$region,tgvutm$mapId)

#prepare the feature data
crt<-merge(crowd[,c("tagId","score","agremnt","sensor","region")],tgvutm[,c("tagId","mapViewId","easting","northing","region","regionTagId","regionMapViewId","regionTaggerId")],by=c("tagId","region"),all.x=T)
crt<-subset(crt,!is.na(mapViewId))
crtm<-merge(crt,views[,c("mapViewId","mapId","region")],by=c("mapViewId","region"),all.x=T)
crtm<-subset(crtm,!is.na(mapId))
crtm$numFeatures<-1
crtm$regionMapId<-paste0(crtm$region,crtm$mapId)

#prepare the tagger crowdrank data
rankdata<-taggers[,c("taggerId","taggerScore","region","regionTaggerId")]

#prepare ML data
mlCounts<-aggregate(numTags~mapId+region+regionMapId,data=subset(tgvutm,taggerId==mid),FUN=sum)
names(mlCounts)<-c("mapId","region","regionMapId","mlcount")
mlMaps<-unique(subset(views,taggerId==mid))
mlMaps$regionMapId<-paste0(mlMaps$region,mlMaps$mapId)
tagsML<-subset(tgvutm,taggerId==mid)

#Subset the data to only those taggers that overlapped with ML - find vews, maps, etc.
taggersSel<-unique(subset(views,regionMapId %in% mlMaps$regionMapId,select=c("taggerId","region","regionMapId","regionMapViewId","regionTaggerId")))
taggersSel<-subset(taggersSel,taggerId!=mid)		#taggers who share maps with ML
tagsSel<-subset(tgvutm,regionTaggerId %in% taggersSel$regionTaggerId) 
viewsSel<-subset(views,regionTaggerId %in% c(unique(taggersSel$regionTaggerId),mid))  

## Estimate (aporx.) the probability of a feature being found on a map and the probability of a tag being placed on a map
totalMapsInspected<-NROW(unique(views$regionMapId))
totalMapsWtags<-NROW(unique(tgvutm$regionMapId))
totalMapsWfeat<-NROW(unique(crtm$regionMapId))
probFeatInMap<-totalMapsWfeat/totalMapsInspected
probTagInMap<-totalMapsWtags/totalMapsInspected
probTagAsFeat<-totalMapsWfeat/totalMapsWtags

##Adding more stringent tagger quality filters
subtaggers<-subset(taggers,taggerScore>0.5 & numApprovedTags > 0)
subtaggers<-subset(subtaggers,regionTaggerId %in% taggersSel$regionTaggerId)   


############ Approximating ProbF and ProbS
## ProbF as the probability of placing a tag on a map: num maps with tags/ total num maps inspected
estProbF<-NROW(unique(tgvutm$regionMapId))/NROW(unique(maps$regionMapId))
## ProbS as the probability of ML finding a seal on a map
estProbS<-NROW(unique(tagsML$regionMapId))/NROW(unique(mlMaps$regionMapId))

#if only a fraction of tags are seals, the approximate correction would be...
(estCorr<-estProbS/estProbF)
#That is, roughtly we count 1.3 seals for every tag placed
## THIS IS NOT what we expected: more tags than seals...


############ SHRINK
## ProbF (for each tagger separately) = number of  tags in maps inspected/total number of maps inspected 
## ProbF is the prob that any tag is found on a map.
## Version 1A: region-specific estimates
gspm<-getTaggerProbabilities_tagsOnlybyTag(rankdata=rankdata,crthr=0.5,taggersSel=subtaggers,tagsSel=tagsSel,tagsML=tagsML,viewsSel=viewsSel,ProbS=estProbS)

## Version 1B: general
gspmG<-getTaggerProbabilities_tagsOnlybyTag_general(rankdata=rankdata,crthr=0.5,taggersSel=subtaggers,tagsSel=tagsSel,tagsML=tagsML,viewsSel=viewsSel,ProbS=estProbS)

## Option 2 - use the ratio of seals to features per tagger to calculate Q
## Version 2A: region-specific estimates
gspmQ<-getTaggerProbabilities_fromQ(rankdata=rankdata,crthr=0.5,taggersSel=subtaggers,tagsSel=tagsSel,tagsML=tagsML,mapsML=mlMaps,viewsSel=viewsSel)

## Version 2B: general
gspmQG<-getTaggerProbabilities_fromQ_general(rankdata=rankdata,crthr=0.5,taggersSel=subtaggers,tagsSel=tagsSel,tagsML=tagsML,mapsML=mlMaps,viewsSel=viewsSel)


############################################################################
## Review the fits and estimates for the correction factors
## Option 1:
summary(gspm$corrFactor); nrow(gspm)
summary(gspmG$corrFactor); nrow(gspmG)
## Safe to say that the values that are region-specific are better and more numerous.

## Checking the distribution of values
ggplot(gspm,aes(x=corrFactor)) + geom_density()
## The above plot shows that it is safe to remove factors of value > 4. There are only 3 of these. Here removing the top 5%
gspm<-subset(gspm,corrFactor < 3.5)
## Review the distribution fit options: log-normal, gamma or Weibull
tstdist<-fitdist(log(gspm$corrFactor),"norm"); plot(tstdist)	#is log-normal best fit??
tstdist<-fitdist(gspm$corrFactor,"gamma"); plot(tstdist)	#I think gamma/weibull is best
tstdist<-fitdist(gspm$corrFactor,"weibull"); plot(tstdist)

cfdist<-fitdist(log(gspm$corrFactor),"norm")$estimate
# lower 95% - mean - upper 95%
paste(round(exp(cfdist[1] - (1.96*cfdist[2])),3),"-",round(exp(cfdist[1]),3),"-",round(exp(cfdist[1] + (1.96*cfdist[2])),3))
# versus
cfdist<-fitdist(gspm$corrFactor,"gamma")$estimate
paste(round(qgamma(0.025,cfdist[1],cfdist[2]),3),"-",round(cfdist[1]/cfdist[2],3),"-",round(qgamma(0.975,cfdist[1],cfdist[2]),3))
## USING the gamma, the number of seals is anywhere from 3% to 300% of the tag counts, mean being 87%

## Option 2:
summary(gspmQ$Qval); nrow(gspmQ)
summary(gspmQG$Qval); nrow(gspmQG)
## Both look good,but the non-general form is tighter

## Checking the distribution of values
ggplot(gspmQ,aes(x=Qval)) + geom_density()
#No need to remove data... and uses the largest number of taggers: 135
## Review the distribution fit options...
tstdist<-fitdist(log(gspmQ$Qval),"norm"); plot(tstdist)	
tstdist<-fitdist(gspmQ$Qval,"gamma"); plot(tstdist)	
tstdist<-fitdist(gspmQ$Qval,"weibull"); plot(tstdist)   #Weibull is certainly the best fit

#SO:
qdist<-fitdist(gspmQG$Qval,"weibull")$estimate
weimean<-qdist[2]*gamma(1+1/qdist[1])
paste(round(qweibull(0.025,qdist[1],qdist[2]),3),"-",round(weimean,3),"-",round(qweibull(0.975,qdist[1],qdist[2]),3))
## The q-values are more reasonable: the number of seals is anywhere from 7% to 157% of the tag counts, mean being 62% - more conservative!

#############################################################################
## Let's calculate seal numbers...
## Using the region-specific form:
#cdf,crthr=0.5,taggers,nSealsFilt=10,tgvutm,maps,overlays,corrMethod,dist="gamma"
countByQ<-getMapEstimates(cdf=gspmQ,crthr=0.5,taggers=subtaggers,nSealsFilt=8,tgvutm=tgvutm,maps=maps,overlays=overlays,corrMethod="taggerQval",dist="weibull")
estByRegionQ<-as.data.frame(countByQ %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)),estNumSeals=ceiling(sum(estNumSeals)),uclNumSeals=ceiling(sum(uclNumSeals))))
estByRegionQ<-rbind(estByRegionQ,data.frame(region="Total",lclNumSeals=round(sum(countByQ$lclNumSeals)),estNumSeals=round(sum(countByQ$estNumSeals)),uclNumSeals=round(sum(countByQ$uclNumSeals))))
print(estByRegionQ)

## If we were to use the general form:
countByQG<-getMapEstimates(cdf=gspmQG,crthr=0.5,taggers=subtaggers,nSealsFilt=8,tgvutm=tgvutm,maps=maps,overlays=overlays,corrMethod="taggerQval",dist="weibull")
estByRegionQG<-as.data.frame(countByQG %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)),estNumSeals=ceiling(sum(estNumSeals)),uclNumSeals=ceiling(sum(uclNumSeals))))
estByRegionQG<-rbind(estByRegionQG,data.frame(region="Total",lclNumSeals=round(sum(countByQG$lclNumSeals)),estNumSeals=round(sum(countByQG$estNumSeals)),uclNumSeals=round(sum(countByQG$uclNumSeals))))
print(estByRegionQG)

## And if we used the region-specific ss values:
countBySS<-getMapEstimates(cdf=gspm,crthr=0.5,taggers=subtaggers,tgvutm=tgvutm,maps=maps,overlays=overlays,corrMethod="taggerCorrFactor",dist="gamma")
estByRegionSS<-as.data.frame(countBySS %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)),estNumSeals=ceiling(sum(estNumSeals)),uclNumSeals=ceiling(sum(uclNumSeals))))
estByRegionSS<-rbind(estByRegionSS,data.frame(region="Total",lclNumSeals=round(sum(countBySS$lclNumSeals)),estNumSeals=round(sum(countBySS$estNumSeals)),uclNumSeals=round(sum(countBySS$uclNumSeals))))
print(estByRegionSS)

## However, these numbers are probably too high - according to this graph,  ~38% too high
countBySS<-merge(countBySS,mlCounts[,c("regionMapId","mlcount")],by="regionMapId",all.x=T)
ggplot(subset(countBySS,!is.na(mlcount)),aes(x=mlcount,y=estNumSeals)) + geom_point() + geom_abline(slope=1,intercept=0, color="blue") + stat_smooth(method="lm",se=F, formula=y~x-1,color="red")
## The sample of surveyors is biased!

## Testing countByQ
countByQ<-merge(countByQ,mlCounts[,c("regionMapId","mlcount")],by="regionMapId",all.x=T)
ggplot(subset(countByQ,!is.na(mlcount)),aes(x=mlcount,y=estNumSeals)) + geom_point() + geom_abline(slope=1,intercept=0, color="blue") + stat_smooth(method="lm",se=F, formula=y~x-1,color="red")
## Seems to underestimate a little

## The general method is best
countByQG<-merge(countByQG,mlCounts[,c("regionMapId","mlcount")],by="regionMapId",all.x=T)
ggplot(subset(countByQG,!is.na(mlcount)),aes(x=mlcount,y=estNumSeals)) + geom_point() + geom_abline(slope=1,intercept=0, color="blue") + stat_smooth(method="lm",se=F, formula=y~x-1,color="red")

## We have some very high tag counts for low true seal counts. Inspecting these...
w<-tags; w$tagCount<-1
w<-merge(w,views[,c("regionMapViewId","regionMapId")],by="regionMapViewId",all.x=T)
## Incredibly, there are still problems in the data like this one:
#  mapId region mapViewId taggerId year regionMapViewId regionMapId regionTaggerId
# 357856    QMA    218627   138101 2011       QMA218627   QMA357856      QMA138101
# 357856    QMA    218621   138101 2011       QMA218621   QMA357856      QMA138101
## NOTE: same taggerId, same mapId, different viewId???

ww<-aggregate(tagCount~regionMapId+regionTaggerId+region+year,w,sum)
ggplot(ww,aes(x=tagCount)) + geom_density()
## So, rarely there are >50 seals in a map. Let's say we try to filter out the tag counts > 100?
## Can we remove these and still have other counts for these maps?
sum(ww$tagCount>100)  #10 cases in 9 maps
NROW(unique(subset(ww,tagCount>100)$regionMapId))/nrow(ww)  #9 maps out of 18,031, or 0.1%

## Can we filter out these crazy counts?
tcdf<-ldply(.data=unique(subset(ww,tagCount>100)$regionMapId), .fun=function(z,ww,views){
			ovis<-NROW(unique(subset(ww,regionMapId==z)$regionTaggerId))
			mval<-subset(ww,tagCount>100 & regionMapId==z)
			mval<-mval[,c("regionTaggerId","region","regionMapId","tagCount")]
			mval$totalViews<-ovis
			return(mval)
		},ww=ww,views=views)

tcdf<-tcdf[order(tcdf$regionMapId),]
print(tcdf)
## Several map counts could not be filtered, but some could.
## Let's look at the counts from maps for each of the above taggers vs other people's counts...
sum(tcdf$totalViews)-2 #one of the maps has high counts from both taggers inspecting it: QMA357856, and both very close counts!
ww$taggerId<-substr(ww$regionTaggerId,4,nchar(ww$regionTaggerId))
tcdf$taggerId<-substr(tcdf$regionTaggerId,4,nchar(tcdf$regionTaggerId))
wwt<-subset(ww,taggerId %in% tcdf$taggerId) 	#all the counts from the high taggers
wwm<-subset(ww,regionMapId %in% wwt$regionMapId)	#all the map counts from maps inspected by the high taggers
wwm$highTagger<-ifelse(wwm$taggerId %in% tcdf$taggerId,"Yes","No")

#let's try to visualize how much hogher these taggers count...
ndf<-aggregate(tagCount~regionMapId+highTagger,subset(wwm,highTagger=="No"),max)
ydf<-subset(wwm,highTagger=="Yes" & regionMapId %in% ndf$regionMapId,select=c("regionMapId","highTagger","tagCount"))
pdf<-rbind(ndf,ydf)
ggplot(pdf,aes(x=tagCount)) + geom_density() + facet_wrap(~highTagger,scales="free")
## The Yes taggers (high taggers) have much larger numbers!

## For these hightaggers, use the mean from the qvals/corrfactors for those in the appropriate gspm and correct their counts individually.


## So, we get numbers again and evaluate again - using the general estimates for Qval:
indivQdf<-subset(gspmQ,regionTaggerId %in% tcdf$regionTaggerId, select=c("regionTaggerId","Qval"))
unique(subset(tcdf,!regionTaggerId %in% indivQdf$regionTaggerId)$regionTaggerId)   ## Do we have value for all? No.
indivQdf<-rbind(indivQdf,data.frame(regionTaggerId=unique(subset(tcdf,!regionTaggerId %in% indivQdf$regionTaggerId)$regionTaggerId),Qval=mean(indivQdf$Qval)))   ## For the missing we use the average of the highTaggers for which we have Qvals/corrFactors
countByQ<-getMapEstimates(cdf=gspmQ,crthr=0.5,taggers=subtaggers,nSealsFilt=8,tgvutm=tgvutm,maps=maps,overlays=overlays,corrMethod="taggerQval",dist="weibull",regional=TRUE,indivCFdf=indivQdf)
estByRegionQ<-as.data.frame(countByQ %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)),estNumSeals=ceiling(sum(estNumSeals)),uclNumSeals=ceiling(sum(uclNumSeals))))
estByRegionQ<-rbind(estByRegionQ,data.frame(region="Total",lclNumSeals=round(sum(countByQ$lclNumSeals)),estNumSeals=round(sum(countByQ$estNumSeals)),uclNumSeals=round(sum(countByQ$uclNumSeals))))
print(estByRegionQ)

## If we were to use the general form:
indivQdf<-subset(gspmQG,taggerId %in% tcdf$taggerId, select=c("taggerId","Qval"))
unique(subset(tcdf,!taggerId %in% indivQdf$taggerId)$taggerId)   
indivQdf<-rbind(indivQdf,data.frame(taggerId=unique(subset(tcdf,!taggerId %in% indivQdf$taggerId)$taggerId),Qval=mean(indivQdf$Qval)))   
countByQG<-getMapEstimates(cdf=gspmQG,crthr=0.5,taggers=subtaggers,nSealsFilt=8,tgvutm=tgvutm,maps=maps,overlays=overlays,corrMethod="taggerQval",dist="weibull",regional=FALSE,indivCFdf=indivQdf)
estByRegionQG<-as.data.frame(countByQG %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)),estNumSeals=ceiling(sum(estNumSeals)),uclNumSeals=ceiling(sum(uclNumSeals))))
estByRegionQG<-rbind(estByRegionQG,data.frame(region="Total",lclNumSeals=round(sum(countByQG$lclNumSeals)),estNumSeals=round(sum(countByQG$estNumSeals)),uclNumSeals=round(sum(countByQG$uclNumSeals))))
print(estByRegionQG)

## And if we used the region-specific ss values:
indivCFdf<-subset(gspm,regionTaggerId %in% tcdf$regionTaggerId, select=c("regionTaggerId","corrFactor"))
unique(subset(tcdf,!regionTaggerId %in% indivCFdf$regionTaggerId)$regionTaggerId)   
indivCFdf<-rbind(indivCFdf,data.frame(regionTaggerId=unique(subset(tcdf,!regionTaggerId %in% indivCFdf$regionTaggerId)$regionTaggerId),corrFactor=mean(indivCFdf$corrFactor)))   
countBySS<-getMapEstimates(cdf=gspm,crthr=0.5,taggers=subtaggers,tgvutm=tgvutm,maps=maps,overlays=overlays,corrMethod="taggerCorrFactor",dist="gamma",regional=TRUE,indivCFdf=indivCFdf)
estByRegionSS<-as.data.frame(countBySS %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)),estNumSeals=ceiling(sum(estNumSeals)),uclNumSeals=ceiling(sum(uclNumSeals))))
estByRegionSS<-rbind(estByRegionSS,data.frame(region="Total",lclNumSeals=round(sum(countBySS$lclNumSeals)),estNumSeals=round(sum(countBySS$estNumSeals)),uclNumSeals=round(sum(countBySS$uclNumSeals))))
print(estByRegionSS)

## However, these numbers are probably too high - according to this graph,  ~38% too high
countBySS<-merge(countBySS,mlCounts[,c("regionMapId","mlcount")],by="regionMapId",all.x=T)
ggplot(subset(countBySS,!is.na(mlcount)),aes(x=mlcount,y=estNumSeals)) + geom_point() + geom_abline(slope=1,intercept=0, color="blue") + stat_smooth(method="lm",se=F, formula=y~x-1,color="red")
summary(countBySS$estNumSeals)
## Underestimates by 30%

## Testing countByQ
countByQ<-merge(countByQ,mlCounts[,c("regionMapId","mlcount")],by="regionMapId",all.x=T)
ggplot(subset(countByQ,!is.na(mlcount)),aes(x=mlcount,y=estNumSeals)) + geom_point() + geom_abline(slope=1,intercept=0, color="blue") + stat_smooth(method="lm",se=F, formula=y~x-1,color="red")
summary(countByQ$estNumSeals)
## Seems to underestimate by 25%

## The general method is still best - does not over-estimate, values are more precise (though less accurate than countBySS)
countByQG<-merge(countByQG,mlCounts[,c("regionMapId","mlcount")],by="regionMapId",all.x=T)
ggplot(subset(countByQG,!is.na(mlcount)),aes(x=mlcount,y=estNumSeals)) + geom_point() + geom_abline(slope=1,intercept=0, color="blue") + stat_smooth(method="lm",se=F, formula=y~x-1,color="red")
summary(countByQG$estNumSeals)
## Seems to underestimate by 20%

## one of two things are happening:
# 1)	The general sample used to generate the correction factors is really not representative and overall results in the underestimate we see; or…
# 2)	The set of MLR maps is a “special case” among all maps where the correction results in an underestimate.
# Since the correction factor is developed for all taggers that overlap with MLR, the underestimate should not happen if “the total is the sum of the parts”. 
# We correct each tagger separately, and expect that the collective should be correct (the graph is indeed against MLR counts). The fact that we see the underestimate 
# indicates that the collective set of values is not “the sum of the parts.” The sample of taggers results in a biased-low mean estimate of the correction factor. 
# In other words, we have evidence for case #1 above. We have no evidence for case #2.

## So...
inflQ<-1.315
estByRegionQ<-as.data.frame(countByQ %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)*inflQ),estNumSeals=ceiling(sum(estNumSeals)*inflQ),uclNumSeals=ceiling(sum(uclNumSeals)*inflQ)))
estByRegionQ<-rbind(estByRegionQ,data.frame(region="Total",lclNumSeals=round(sum(countByQ$lclNumSeals)*inflQ),estNumSeals=round(sum(countByQ$estNumSeals)*inflQ),uclNumSeals=round(sum(countByQ$uclNumSeals)*inflQ)))
print(estByRegionQ)
ggplot(subset(countByQ,!is.na(mlcount)),aes(x=mlcount,y=estNumSeals*inflQ)) + geom_point() + geom_abline(slope=1,intercept=0, color="blue") + stat_smooth(method="lm",se=F, formula=y~x-1,color="red")
summary(countByQ$estNumSeals*inflQ)

## Or...
inflQG<-1.245
estByRegionQG<-as.data.frame(countByQG %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)*inflQG),estNumSeals=ceiling(sum(estNumSeals)*inflQG),uclNumSeals=ceiling(sum(uclNumSeals)*inflQG)))
estByRegionQG<-rbind(estByRegionQG,data.frame(region="Total",lclNumSeals=round(sum(countByQG$lclNumSeals)*inflQG),estNumSeals=round(sum(countByQG$estNumSeals)*inflQG),uclNumSeals=round(sum(countByQG$uclNumSeals)*inflQG)))
print(estByRegionQG)
ggplot(subset(countByQG,!is.na(mlcount)),aes(x=mlcount,y=estNumSeals*inflQG)) + geom_point() + geom_abline(slope=1,intercept=0, color="blue") + stat_smooth(method="lm",se=F, formula=y~x-1,color="red")
summary(countByQG$estNumSeals*inflQG)

## Adjust countByQ and countByQG
countByQ$estNumSeals<-round(countByQ$estNumSeals*inflQ)
countByQ$uclNumSeals<-round(countByQ$uclNumSeals*inflQ)
countByQ$lclNumSeals<-round(countByQ$lclNumSeals*inflQ)

countByQG$estNumSeals<-round(countByQG$estNumSeals*inflQG)
countByQG$uclNumSeals<-round(countByQG$uclNumSeals*inflQG)
countByQG$lclNumSeals<-round(countByQG$lclNumSeals*inflQG)

## Need to add year and hour of the day to the results (Ross Sea time - NZDT), and map lat/lon
## PDT to NZDT is +19 hrs
## PST to NZDT is +21 hrs
nrow(countByQ)
countByQ<-merge(countByQ,unique(overlays[,c("overlayId","acquisition_date")]),by="overlayId",all.x=TRUE)
nrow(countByQ)
countByQ$originHour<-as.integer(format(countByQ$acquisition_date,"%H"))
countByQ$originTZ<-format(countByQ$acquisition_date,"%Z")
countByQ$RShour<-ifelse(countByQ$originTZ=="PST",(countByQ$originHour+21)%%24,
		ifelse(countByQ$originTZ=="PDT",countByQ$originHour+19%%24,countByQ$originTZ))
countByQ<-merge(countByQ,maps[,c("regionMapId","mapcoords.x1","mapcoords.x2")],by="regionMapId",all.x=TRUE)
nrow(countByQ)


save(estByRegionQ, estByRegionQG, countByQ, countByQG, file=paste0(pathToLocalGit,"estimatesByMap.RData"))

