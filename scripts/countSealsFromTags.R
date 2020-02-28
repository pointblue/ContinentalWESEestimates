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
## The above plot shows that it is safe to remove factors of value > 3. There are only 3 of these. 
gspm<-subset(gspm,corrFactor < 3)
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
qdist<-fitdist(gspmQ$Qval,"weibull")$estimate
weimean<-qdist[2]*gamma(1+1/qdist[1])
paste(round(qweibull(0.025,qdist[1],qdist[2]),3),"-",round(weimean,3),"-",round(qweibull(0.975,qdist[1],qdist[2]),3))
## The q-values are more reasonable: the number of seals is anywhere from 7% to 157% of the tag counts, mean being 62% - more conservative!
# versus
qdist<-fitdist(gspmQG$Qval,"weibull")$estimate
weimean<-qdist[2]*gamma(1+1/qdist[1])
paste(round(qweibull(0.025,qdist[1],qdist[2]),3),"-",round(weimean,3),"-",round(qweibull(0.975,qdist[1],qdist[2]),3))
## The general approach is more conservative, and better than using ss

#############################################################################
## For map estimates, setting the CI to 90% per Nadav's request
cival=90
## Let's calculate seal numbers...
## Using the region-specific form:
#cdf,crthr=0.5,taggers,nSealsFilt=10,tgvutm,maps,overlays,corrMethod,dist="gamma"
countByQ<-getMapEstimates(cdf=gspmQ,crthr=0.5,taggers=subtaggers,nSealsFilt=8,tgvutm=tgvutm,maps=maps,overlays=overlays,corrMethod="taggerQval",dist="weibull",cival=cival)
estByRegionQ<-as.data.frame(countByQ %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)),estNumSeals=ceiling(sum(estNumSeals)),uclNumSeals=ceiling(sum(uclNumSeals))))
estByRegionQ<-rbind(estByRegionQ,data.frame(region="Total",lclNumSeals=round(sum(countByQ$lclNumSeals)),estNumSeals=round(sum(countByQ$estNumSeals)),uclNumSeals=round(sum(countByQ$uclNumSeals))))
print(estByRegionQ)

## If we were to use the general form:
countByQG<-getMapEstimates(cdf=gspmQG,crthr=0.5,taggers=subtaggers,nSealsFilt=8,tgvutm=tgvutm,maps=maps,overlays=overlays,corrMethod="taggerQval",dist="weibull",cival=cival)
estByRegionQG<-as.data.frame(countByQG %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)),estNumSeals=ceiling(sum(estNumSeals)),uclNumSeals=ceiling(sum(uclNumSeals))))
estByRegionQG<-rbind(estByRegionQG,data.frame(region="Total",lclNumSeals=round(sum(countByQG$lclNumSeals)),estNumSeals=round(sum(countByQG$estNumSeals)),uclNumSeals=round(sum(countByQG$uclNumSeals))))
print(estByRegionQG)

## And if we used the region-specific ss values:
countBySS<-getMapEstimates(cdf=gspm,crthr=0.5,taggers=subtaggers,tgvutm=tgvutm,maps=maps,overlays=overlays,corrMethod="taggerCorrFactor",dist="gamma",cival=cival)
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
## So, rarely there are >50 seals in a map. Let's say we try to filter out the tag counts > 100?
## Can we remove these and still have other counts for these maps?
sum(ww$tagCount>100)  #44 cases in 40 maps
NROW(unique(subset(ww,tagCount>100)$regionMapId))/nrow(ww)  #40 maps out of 18,031, or 0.1%

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

#let's try to visualize how much higher these taggers count...
ndf<-aggregate(tagCount~regionMapId+highTagger,subset(wwm,highTagger=="No"),max)
ydf<-subset(wwm,highTagger=="Yes" & regionMapId %in% ndf$regionMapId,select=c("regionMapId","highTagger","tagCount"))
pdf<-rbind(ndf,ydf)
ggplot(pdf,aes(x=tagCount)) + geom_histogram(binwidth=20) + facet_wrap(~highTagger,scales="free")
## The Yes taggers (high taggers) have much larger numbers!
## But these are not unusual when their counts are low
ggplot(subset(pdf,highTagger=="Yes" & tagCount<100),aes(x=tagCount)) + geom_histogram(binwidth=20)

## So, we want to remove their high counts for a map if the map has been inspected by at least 3 other taggers
## We then use their own Qval for when they cannot be removed

## Removing the removable high counts
tcdff<-ldply(.data=unique(tcdf$regionMapId), .fun=function(m,tcdf){
			ttdf<-subset(tcdf,regionMapId==m)
			nr<-nrow(ttdf)
			nh<-unique(ttdf$totalViews)
			if(nh<(nr+2)){
				rdf<-ttdf
			}else{
				rdf<-data.frame(regionTaggerId=NA, region=NA, regionMapId=NA, tagCount=NA, totalViews=NA)
			}
			
		},tcdf=tcdf)
tcdff<-na.omit(tcdff)
tcdfr<-subset(tcdf,!regionMapId %in% tcdff$regionMapId)
tcdfr$taggerMap<-paste0(tcdfr$regionTaggerId,":",tcdfr$regionMapId)
## Confirming:
orv<-nrow(tgvutm)
tgvutm<-subset(tgvutm,!paste0(regionTaggerId,":",regionMapId) %in% tcdfr$taggerMap)
nrow(tgvutm)==orv-sum(tcdfr$tagCount)	## TRUE


## So, we get numbers again and evaluate again - using the general estimates for Qval plus the individual estimates for the high count taggers:
indivQdf<-subset(gspmQ,regionTaggerId %in% tcdff$regionTaggerId, select=c("regionTaggerId","Qval"))
unique(subset(tcdff,!regionTaggerId %in% indivQdf$regionTaggerId)$regionTaggerId)   ## Do we have value for all? No.
indivQdf<-rbind(indivQdf,data.frame(regionTaggerId=unique(subset(tcdff,!regionTaggerId %in% indivQdf$regionTaggerId)$regionTaggerId),Qval=mean(indivQdf$Qval)))   ## For the missing we use the average of the highTaggers for which we have Qvals/corrFactors
countByQ<-getMapEstimates(cdf=gspmQ,crthr=0.5,taggers=subtaggers,nSealsFilt=8,tgvutm=tgvutm,maps=maps,overlays=overlays,corrMethod="taggerQval",dist="weibull",regional=TRUE,indivCFdf=indivQdf,cival=cival)
estByRegionQ<-as.data.frame(countByQ %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)),estNumSeals=ceiling(sum(estNumSeals)),uclNumSeals=ceiling(sum(uclNumSeals))))
estByRegionQ<-rbind(estByRegionQ,data.frame(region="Total",lclNumSeals=round(sum(countByQ$lclNumSeals)),estNumSeals=round(sum(countByQ$estNumSeals)),uclNumSeals=round(sum(countByQ$uclNumSeals))))
print(estByRegionQ)

## If we were to use the general form:
indivQdf<-subset(gspmQG,taggerId %in% tcdff$taggerId, select=c("taggerId","Qval"))
unique(subset(tcdff,!taggerId %in% indivQdf$taggerId)$taggerId)   
indivQdf<-rbind(indivQdf,data.frame(taggerId=unique(subset(tcdff,!taggerId %in% indivQdf$taggerId)$taggerId),Qval=mean(indivQdf$Qval)))   
countByQG<-getMapEstimates(cdf=gspmQG,crthr=0.5,taggers=subtaggers,nSealsFilt=8,tgvutm=tgvutm,maps=maps,overlays=overlays,corrMethod="taggerQval",dist="weibull",regional=FALSE,indivCFdf=indivQdf,cival=cival)
estByRegionQG<-as.data.frame(countByQG %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)),estNumSeals=ceiling(sum(estNumSeals)),uclNumSeals=ceiling(sum(uclNumSeals))))
estByRegionQG<-rbind(estByRegionQG,data.frame(region="Total",lclNumSeals=round(sum(countByQG$lclNumSeals)),estNumSeals=round(sum(countByQG$estNumSeals)),uclNumSeals=round(sum(countByQG$uclNumSeals))))
print(estByRegionQG)

## And if we used the region-specific ss values:
indivCFdf<-subset(gspm,regionTaggerId %in% tcdff$regionTaggerId, select=c("regionTaggerId","corrFactor"))
unique(subset(tcdff,!regionTaggerId %in% indivCFdf$regionTaggerId)$regionTaggerId)   
indivCFdf<-rbind(indivCFdf,data.frame(regionTaggerId=unique(subset(tcdff,!regionTaggerId %in% indivCFdf$regionTaggerId)$regionTaggerId),corrFactor=mean(indivCFdf$corrFactor)))   
countBySS<-getMapEstimates(cdf=gspm,crthr=0.5,taggers=subtaggers,tgvutm=tgvutm,maps=maps,overlays=overlays,corrMethod="taggerCorrFactor",dist="gamma",regional=TRUE,indivCFdf=indivCFdf,cival=cival)
estByRegionSS<-as.data.frame(countBySS %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)),estNumSeals=ceiling(sum(estNumSeals)),uclNumSeals=ceiling(sum(uclNumSeals))))
estByRegionSS<-rbind(estByRegionSS,data.frame(region="Total",lclNumSeals=round(sum(countBySS$lclNumSeals)),estNumSeals=round(sum(countBySS$estNumSeals)),uclNumSeals=round(sum(countBySS$uclNumSeals))))
print(estByRegionSS)

## However, these numbers are probably too low - according to this graph
countBySS<-merge(countBySS,mlCounts[,c("regionMapId","mlcount")],by="regionMapId",all.x=T)
ggplot(subset(countBySS,!is.na(mlcount)),aes(x=mlcount,y=estNumSeals)) + geom_point() + geom_abline(slope=1,intercept=0, color="blue") + stat_smooth(method="lm",se=F, formula=y~x-1,color="red")
summary(countBySS$estNumSeals)
## Underestimates by 34%

## Testing countByQ
countByQ<-merge(countByQ,mlCounts[,c("regionMapId","mlcount")],by="regionMapId",all.x=T)
ggplot(subset(countByQ,!is.na(mlcount)),aes(x=mlcount,y=estNumSeals)) + geom_point() + geom_abline(slope=1,intercept=0, color="blue") + stat_smooth(method="lm",se=F, formula=y~x-1,color="red")
summary(countByQ$estNumSeals)
## Seems to underestimate by 30%

## The general method is still best - does not over-estimate, values are more precise (though less accurate than countBySS)
countByQG<-merge(countByQG,mlCounts[,c("regionMapId","mlcount")],by="regionMapId",all.x=T)
ggplot(subset(countByQG,!is.na(mlcount)),aes(x=mlcount,y=estNumSeals)) + geom_point() + geom_abline(slope=1,intercept=0, color="blue") + stat_smooth(method="lm",se=F, formula=y~x-1,color="red")
summary(countByQG$estNumSeals)
## Seems to underestimate by 34%

## one of two things are happening:
# 1)	The general sample used to generate the correction factors is really not representative and overall results in the underestimate we see; or…
# 2)	The set of MLR maps is a “special case” among all maps where the correction results in an underestimate.
# Since the correction factor is developed for all taggers that overlap with MLR, the underestimate should not happen if “the total is the sum of the parts”. 
# We correct each tagger separately, and expect that the collective should be correct (the graph is indeed against MLR counts). The fact that we see the underestimate 
# indicates that the collective set of values is not “the sum of the parts.” The sample of taggers results in a biased-low mean estimate of the correction factor. 
# In other words, we have evidence for case #1 above. We have no evidence for case #2.

## So...
inflQ<-1.385
estByRegionQ<-as.data.frame(countByQ %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)*inflQ),estNumSeals=ceiling(sum(estNumSeals)*inflQ),uclNumSeals=ceiling(sum(uclNumSeals)*inflQ)))
estByRegionQ<-rbind(estByRegionQ,data.frame(region="Total",lclNumSeals=round(sum(countByQ$lclNumSeals)*inflQ),estNumSeals=round(sum(countByQ$estNumSeals)*inflQ),uclNumSeals=round(sum(countByQ$uclNumSeals)*inflQ)))
print(estByRegionQ)
pq<-ggplot(subset(countByQ,!is.na(mlcount)),aes(x=mlcount,y=estNumSeals*inflQ)) + geom_point() + geom_abline(slope=1,intercept=0, color="blue") + stat_smooth(method="lm",se=F, formula=y~x-1,color="red")
summary(countByQ$estNumSeals*inflQ)

## Or...
inflQG<-1.48
estByRegionQG<-as.data.frame(countByQG %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)*inflQG),estNumSeals=ceiling(sum(estNumSeals)*inflQG),uclNumSeals=ceiling(sum(uclNumSeals)*inflQG)))
estByRegionQG<-rbind(estByRegionQG,data.frame(region="Total",lclNumSeals=round(sum(countByQG$lclNumSeals)*inflQG),estNumSeals=round(sum(countByQG$estNumSeals)*inflQG),uclNumSeals=round(sum(countByQG$uclNumSeals)*inflQG)))
print(estByRegionQG)
pqg<-ggplot(subset(countByQG,!is.na(mlcount)),aes(x=mlcount,y=estNumSeals*inflQG)) + geom_point() + geom_abline(slope=1,intercept=0, color="blue") + stat_smooth(method="lm",se=F, formula=y~x-1,color="red")
summary(countByQG$estNumSeals*inflQG)

## Or..
inflSS<-1.54
estByRegionSS<-as.data.frame(countBySS %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)*inflSS),estNumSeals=ceiling(sum(estNumSeals)*inflSS),uclNumSeals=ceiling(sum(uclNumSeals)*inflSS)))
estByRegionSS<-rbind(estByRegionSS,data.frame(region="Total",lclNumSeals=round(sum(countBySS$lclNumSeals)*inflSS),estNumSeals=round(sum(countBySS$estNumSeals)*inflSS),uclNumSeals=round(sum(countBySS$uclNumSeals)*inflSS)))
print(estByRegionSS)
pss<-ggplot(subset(countBySS,!is.na(mlcount)),aes(x=mlcount,y=estNumSeals*inflSS)) + geom_point() + geom_abline(slope=1,intercept=0, color="blue") + stat_smooth(method="lm",se=F, formula=y~x-1,color="red")
summary(countBySS$estNumSeals*inflSS)

## countByQ is the top choice!
## Some estimates are still a bit high:
summary(countByQ$estNumSeals*inflQ)


## Adjust countByQ and countByQG - we will use countByQ
countByQ$estNumSeals<-round(countByQ$estNumSeals*inflQ)
countByQ$uclNumSeals<-round(countByQ$uclNumSeals*inflQ)
countByQ$lclNumSeals<-round(countByQ$lclNumSeals*inflQ)

## Need to add year and hour of the day to the results (Ross Sea time - NZDT), and map lat/lon
nrow(countByQ)
countByQ<-merge(countByQ,unique(overlays[,c("overlayId","acquisition_date")]),by="overlayId",all.x=TRUE)
nrow(countByQ)
countByQ$originHour<-format(countByQ$acquisition_date,"%H")
countByQ$originTZ<-format(countByQ$acquisition_date,"%Z")
countByQ$RSdate<-as.POSIXlt(countByQ$acquisition_date,tz="Pacific/Auckland")	#convert to RS time zone
countByQ$RShour<-format(countByQ$RSdate,"%H")
countByQ<-merge(countByQ,maps[,c("regionMapId","mapcoords.x1","mapcoords.x2")],by="regionMapId",all.x=TRUE)
pdf<-countByQ
pdf$abBin<-ifelse(pdf$estNumSeals<20,"Under 20",ifelse(pdf$estNumSeals>19 & pdf$estNumSeals<51,"20 to 50","Over 50"))
pdf$abBin2<-ifelse(pdf$estNumSeals<40,"Under 40","Over 40")
pdf$wvq<-ifelse(pdf$satId=="QB02","QB","WV")
ggplot(pdf,aes(x=estNumSeals)) + geom_histogram(aes(fill=satId),position="dodge",binwidth=4) + facet_wrap(~abBin2,scales="free") + labs(x="Seals per map",y="Number of maps")

# It does not look like it's a satellite thing the high numbers...

# We can see that under 40 the distributions among satIds are the same for WV01, for WV02, and for QB (only 9 maps from GE01)
# The problem is the counts in 30 maps above 40 seals (out of 17,453!)
aggregate(regionMapId~satId,pdf,NROW)

# SO: maybe we don't need to worry about these maps with > 40 seals, or we can shrink them to the mean of the hour:
p1<-ggplot(countByQ,aes(x=RShour,y=estNumSeals))+geom_boxplot() + geom_hline(yintercept=40)
# Note the outliers
# The value for WV01 on the 21st hour is NOT credible, from only 1 tagger
g<-subset(tags,regionMapViewId %in% unique(subset(views,regionMapId=="RSS1255032")$regionMapViewId))
tt<-aggregate(tagId~regionTaggerId,g,NROW);names(tt)<-c("regionTaggerId","tagCount");print(tt)

# The value for WV1 for the 4th hour is also not credible - only 1 tagger
g<-subset(tags,regionMapViewId %in% unique(subset(views,regionMapId=="RSS938501")$regionMapViewId))
tt<-aggregate(tagId~regionTaggerId,g,NROW);names(tt)<-c("regionTaggerId","tagCount");print(tt)


# More generally, we can't trust any value with mean tags > 100, and we can't trust any value with tags > 40 from a single tagger
outs<-subset(countByQ,estNumSeals>40)
outvals<-ldply(.data=outs$regionMapId,.fun=function(mm,outs,tags,views){
			g<-subset(tags,regionMapViewId %in% unique(subset(views,regionMapId==mm)$regionMapViewId))
			tou<-subset(outs,regionMapId==mm)
			ns<-tou$estNumSeals; nt<-tou$numTaggers; hr<-tou$RShour
			adf<-aggregate(tagId~regionTaggerId,g,NROW)
			tdf<-data.frame(Hour=hr,regionMapId=mm,estNumSeals=ns,numTaggers=nt,minTags=min(adf$tagId),maxTags=max(adf$tagId),meanTags=mean(adf$tagId))
			return(tdf)
		},outs=outs,tags=tags,views=views)
outvals$Hour<-as.character(outvals$Hour)
outvals<-subset(outvals,(meanTags>100) | (meanTags>40 & numTaggers==1))

# Now we correct the counts for these outlier map counts by using the 99% value of a negative binomial distribution applied to the set of values for the hour
counts<-countByQ
for(mm in outvals$regionMapId){
	hv<-subset(outvals,regionMapId==mm)$Hour
	outh<-subset(counts, RShour==hv)	
	dv<-fitdist(outh$estNumSeals,"nbinom")$estimate
	nmv<-qnbinom(0.99,mu=dv[2],size=dv[1])
	counts$estNumSeals<-ifelse(counts$regionMapId==mm,nmv,counts$estNumSeals)
}

p2<-ggplot(counts,aes(x=RShour,y=estNumSeals))+geom_boxplot(aes(color=satId)) + geom_hline(yintercept=40)
p3<-ggplot(counts,aes(x=RShour,y=estNumSeals))+geom_boxplot() + geom_hline(yintercept=40)

## And estimate again...
# NOTE: we inflated the counts already, don't need to do it again - this is to adjust for the 30 over-estimated maps
estByRegionQ<-as.data.frame(counts[,c("region","lclNumSeals","estNumSeals","uclNumSeals")] %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)),estNumSeals=ceiling(sum(estNumSeals)),uclNumSeals=ceiling(sum(uclNumSeals))))
estByRegionQ<-rbind(estByRegionQ,data.frame(region="Total",lclNumSeals=round(sum(countByQ$lclNumSeals)),estNumSeals=round(sum(countByQ$estNumSeals)),uclNumSeals=round(sum(countByQ$uclNumSeals))))
print(estByRegionQ)
pq<-ggplot(subset(counts,!is.na(mlcount)),aes(x=mlcount,y=estNumSeals)) + geom_point() + geom_abline(slope=1,intercept=0, color="blue") + stat_smooth(method="lm",se=F, formula=y~x-1,color="red")
summary(counts$estNumSeals)

##HERE!!
## Careful! this takes time!
ntdf<-ldply(.data=unique(countByQ$regionMapId),.fun=function(m,views,tags){
			mv<-subset(views,regionMapId==m)
			numViews<-nrow(mv)
			mt<-subset(tags,regionMapViewId %in% unique(mv$regionMapViewId))
			numTaggers<-NROW(unique(mt$taggerId))
			totalTags<-nrow(mt)
			tdf<-data.frame(regionMapId=m,numViews=numViews,totalTags=totalTags,avgTags=(totalTags/numTaggers))
		},views=views,tags=tags)

counts<-merge(counts,ntdf,by="regionMapId",all.x=T)

#Correcting for hour effect by fitting the sinusoidal:
counts$scaledTotalTags<-scale(counts$totalTags); counts$scaledTotalTags<-as.numeric(counts$scaledTotalTags)
counts$scaledAvgTags<-scale(counts$avgTags); counts$scaledAvgTags<-as.numeric(counts$scaledAvgTags)

## Our base model:
mdlb<-lm(estNumSeals~region+satId+scaledAvgTags+numViews+numTaggers+region*scaledAvgTags+region:scaledTotalTags+satId*scaledAvgTags+year,data=counts) 	
#No year effects - we add the sinusoidal effect to this

hourdf<-data.frame(model="base",aicv=AIC(mdlb), bicv=BIC(mdlb), lmdlv=logLik(mdlb), dfv=nrow(counts)-mdlb$df.residual)

for(hh in 3:21){
	df<-counts
	df$numHour<-((as.integer(counts$RShour)-12) %% hh)/hh  #We test various moduli
	sindf<-unique(df[,c("RShour","numHour")])
	sindf$sinH<-sin(2*pi*sindf$numHour)
	df<-merge(df,sindf[,c("RShour","sinH")],by="RShour",all.x=T)
	
	mdlh<-lm(estNumSeals~region+satId+scaledAvgTags+numViews+numTaggers+region*scaledAvgTags+region:scaledTotalTags+satId*scaledAvgTags+sinH+I(sinH^2),data=df)
	hdf<-data.frame(model=paste0("mod",hh),aicv=AIC(mdlh), bicv=BIC(mdlh), lmdlv=logLik(mdlh), dfv=nrow(counts)-mdlh$df.residual)
	hourdf<-rbind(hourdf,hdf)
	
}
hourdf<-hourdf[order(hourdf$aicv),]

#looks like modulus 4 is best! Since the colony and island/mainland models were built with a 12hr modulus, we are here removing the hour effect
df<-counts
df$numHour<-((as.integer(counts$RShour)-12) %% 4)/4  #We use 4 as best modulus
sindf<-unique(df[,c("RShour","numHour")])
sindf$sinH<-sin(2*pi*sindf$numHour)
df<-merge(df,sindf[,c("RShour","sinH")],by="RShour",all.x=T)
mdlh<-lm(estNumSeals~region+satId+scaledAvgTags+numViews+numTaggers+region*scaledAvgTags+region:scaledTotalTags+satId*scaledAvgTags+sinH+I(sinH^2),data=df)

#For each hour's records, calculate the error (distance to predicted), then predict to a particular hour, and add the error to the prediction
df$predicted<-predict(mdlh)
df$resid<-mdlh$residuals

#########################################
##Find the best hour to predict to:
dfh<-unique(df[,c("RShour","sinH")])
dfh$region<-"RSS"; dfh$satId<-"WV02"; dfh$scaledAvgTags<-0; dfh$scaledTotalTags<-0; dfh$numViews<-9; dfh$numTaggers<-2
dfh$predicted<-predict(mdlh,newdata=dfh)
dfh[order(dfh$predicted),]
## RShour=20, sinH=0
#########################################

newdat<-df[,c("regionMapId","RShour","region","satId","scaledAvgTags","scaledTotalTags","numViews","numTaggers","resid","predicted","estNumSeals","lclNumSeals","uclNumSeals")]
newdat$propL<-newdat$lclNumSeals/newdat$estNumSeals
newdat$propU<-newdat$uclNumSeals/newdat$estNumSeals
newdat$sinH<-0
newdat$predictH<-predict(mdlh,newdata=newdat)
newdat$estimateH<-round(newdat$predictH+newdat$resid)
newdat$lowerH<-round(newdat$estimateH*newdat$propL)
newdat$upperH<-round(newdat$estimateH*newdat$propU)
#Print these!
estByRegionQ<-as.data.frame(newdat[,c("region","lowerH","estimateH","upperH")] %>% group_by(region) %>% dplyr::summarize(lower=ceiling(sum(lowerH)),estimate=ceiling(sum(estimateH)),upper=ceiling(sum(upperH))))
estByRegionQ<-rbind(estByRegionQ,data.frame(region="Total",lower=round(sum(newdat$lowerH)),estimate=round(sum(newdat$estimateH)),upper=round(sum(newdat$upperH))))
print(estByRegionQ)

df<-merge(df,newdat[,c("regionMapId","estimateH","lowerH","upperH")],by="regionMapId",all.x=T)
df<-df[,which(!names(df) %in% c("wgtEstNumSeals","wgtUclNumSeals","wgtLclNumSeals"))]


#Finally, calculte the estimates using the island model, then predict using the colony model assuming all are island, then all mainland, and then do a 90%/10% average
## Can use dff or df - both correcting models include sinH and predict the detection rate: estimate/actual count.  
## So, once we have the predicted detection rate, we divide by the estimate to obtain the inflated count
## We predict assuming each map is from each of the colonies, and then average, weighing the mainland colonies at 10%, all others at 90%
## The island model is: mdlIsl<-lm(detRate~scaledNumTags+sinH+I(sinH^2)+Island+acYear,data=numSealsF)
## We predict 100% as islands, 100% as mainland, then add up with 90% weight to island colonies.
## All this is done in the adjustCounts_withDetectionRate.R



save(estByRegionQ, estByRegionQG, countByQ, countByQG, counts, df, ntdf, file=paste0(pathToLocalGit,"estimatesByMap_unadjusted.RData"))

