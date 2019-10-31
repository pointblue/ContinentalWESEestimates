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

libs<-c("ggplot2","raster","rgdal","plyr","sp","rgeos","fitdistrplus")
lapply(libs, require, character.only = TRUE)

pathToLocalGit<-"C:/Users/lsalas/git/ContinentalWESEestimates/data/ContinentalEst2/"

######################  FUNCTIONS WE'LL NEED
## Convert to UTM
# data is the table with the geo data
# latfield is the string naming the latitude field
# lonfield is the string naming the longitude field
convertToUTM<-function(data,latfield,lonfield){
	
	coordinates(data)<-c(lonfield,latfield)
	proj4string(data) <- CRS("+proj=longlat +datum=WGS84")
	utmdata <- spTransform(data, CRS("+proj=utm +zone=58 +south ellps=WGS84"))
	utmdata<-as.data.frame(utmdata)
	names(utmdata)<-gsub(lonfield,"easting",names(utmdata))
	names(utmdata)<-gsub(latfield,"northing",names(utmdata))
	
	return(utmdata)
}

## Get the number of map by rank of tagger
# x is the filter rank value
# df is the table that lists the maps, taggers that inspected these, and their rank value
getNumMapsByRank<-function(x,df){
	mlMaps<-unique(subset(df,taggerId==21758509,select="mapId")$mapId)
	dfml<-subset(df,mapId %in% mlMaps)
	#keep only those where ML is not the sole inspector
	mlplus<-sapply(X=dfml$mapId,FUN=function(X,dfx){
				gg<-subset(dfx,mapId==X);
				mapchk<-ifelse(NROW(unique(gg$taggerId))>1,X,NA);
				return(mapchk)
			},dfx=dfml)
	dfml<-subset(dfml,mapId %in% na.omit(mlplus))
	
	dff<-subset(df,taggerScore>=x)
	numaps<-NROW(unique(dff$mapId))
	dffml<-subset(dfml,taggerScore>=x)
	numaps<-NROW(unique(dff$mapId));numlmaps<-NROW(unique(dffml$mapId))
	nutaggers<-NROW(unique(dff$taggerId));numltaggers<-NROW(unique(dffml$taggerId))
	resdf<-data.frame(threshold=x,Nmaps=numaps,NMLmaps=numlmaps,Ntaggers=nutaggers,NMLtaggers=numltaggers)
	return(resdf)
}

# x is the value to logit-transform to
toLogit<-function(x){
	cx<-ifelse(x==0,0.00001,
			ifelse(x==1,0.99999,x))
	lgx<-log(cx)-log(1-cx)
	return(lgx)
}

# x is the value to logit-transform from
fromLogit<-function(x){
	bt<-exp(x)/(1+exp(x))
	return(bt)
}

## Find the number of seals the tager correctly identified (compared to ML)
# taggerTags is the data.frame with all the tagger's tags in maps shared with ML, including the east/north of each tag
# mlTags is the data.frame with all of ML tags, including the east/north of each tag
getTaggerSS<-function(taggerTags,mlTags){
	#loop through each tag and see if it's within 3m of a ML tag
	mlTags$regionTagId<-paste0(mlTags$region,mlTags$tagId)
	ss<-0
	for(rr in 1:nrow(taggerTags)){
		mlTags$taggerEast<-taggerTags[rr,"easting"]
		mlTags$taggerNorth<-taggerTags[rr,"northing"]
		mlTags$dist<-sqrt(((mlTags$easting-mlTags$taggerEast)^2)+((mlTags$northing-mlTags$taggerNorth)^2))
		featTags<-subset(mlTags,dist <= 3)
		if(nrow(featTags)>0){
			#take the closest and remove from mlTags
			ss<-ss+1
			featTags<-featTags[order(featTags$dist),]
			topTag<-as.character(featTags[1,"regionTagId"])
			mlTags<-subset(mlTags,regionTagId!=topTag)
			if(nrow(mlTags)==0){break}
		}
	}
	return(ss)
}

## ALT10: ProbF = number of  tags in maps inspected/total number of maps inspected			viewsSel	region:taggertTags
## ProbF is the prob that any tag is found on a map.
## corrProbF should be the probability that any tagger will tag a map, which is the ratio of all tags/all maps
getTaggerProbabilities_tagsOnlybyTag<-function(rankdata,crthr,taggersSel,tagsSel,tagsML,crtm,viewsSel,ProbS,corrProbF=1){
	if(!is.na(crthr)){
		tsdf<-subset(rankdata,taggerScore>=crthr);	#only taggers matching or besting the threshold
		taggers<-subset(taggersSel,taggersSel$taggerId %in% tsdf$taggerId);	#those taggers that shared maps with Michelle and who also have a score higher or equal to threshold
		
	}else{
		taggers<-taggersSel
	}
	
	#calculate the probabilities for each tagger exceding the threshold crthr - need tags and feature geo
	taggerProbs<-ldply(.data=unique(taggers$taggerId),.fun=function(tt,viewsSel,crtm,tagsSel,tagsML){
				#calculating ProbF...(number of maps with tags among the maps inspected by this tagger)
				#need the number of maps inspected by this tagger
				taggerMaps<-unique(subset(viewsSel,taggerId==tt)$regionMapId)
				tagger_nMaps<-NROW(taggerMaps)
				#determine how many of these have tags
				taggerTags<-subset(tagsSel,taggerId==tt);totalTags<-nrow(taggerTags)
				taggerTags<-merge(taggerTags,viewsSel[,c("mapViewId","mapId","regionMapId","region")],by=c("mapViewId","region"),all.x=T)
				tagger_nTagMaps<-NROW(unique(taggerTags$regionMapId))
				mlTags<-subset(tagsML,regionMapId %in% taggerMaps)	#if this is 0, it's because ML did not find a seal in any of the maps this tagger inspected.
				totalSeals<-nrow(mlTags)
				ProbF<-totalTags/tagger_nMaps
				
				if(totalTags>0){	#features were found in maps inspected by the tagger
					#calculating ProbFS...
					#need the list of maps from this tagger also shared with ML - this is taggerMaps
					#need the list of tags in these maps
					
					if(nrow(mlTags)==0){
						ttdf<-data.frame(taggerId=tt,ss=0,nTags=totalTags,nSeals=totalSeals,nMaps=tagger_nMaps,ProbFS=0,ProbF=0)
					}else{
						ss<-getTaggerSS(taggerTags,mlTags)
						ProbFS<-ss/totalSeals
						
						ttdf<-data.frame(taggerId=tt,ss=ss,nTags=totalTags,nSeals=totalSeals,nMaps=tagger_nMaps,ProbFS=ProbFS,ProbF=ProbF)
					}
				}else{	#no features found in maps inspected by this tagger
					ttdf<-data.frame(taggerId=tt,ss=0,nTags=totalTags,nSeals=0,nMaps=tagger_nMaps,ProbFS=0,ProbF=0)
				}
				
				
				return(ttdf)},viewsSel=viewsSel,crtm=crtm,tagsSel=tagsSel,tagsML=tagsML)
	
	taggerProbs<-subset(taggerProbs,ProbFS>0)
	if(nrow(taggerProbs)>0){
		taggerProbs$ProbS<-ProbS
		taggerProbs$ProbSF<-taggerProbs$ProbS*taggerProbs$ProbFS/taggerProbs$ProbF	
		taggerProbs$corrFactor<-taggerProbs$ProbSF/taggerProbs$ProbFS
	}else{ #none of the computed probabilities permits the estimation of corrFactor
		taggerProbs<-NA
	}
	
	return(taggerProbs)
}


getMapEstimates_byTag<-function(gspm,tgvutm,maps,overlays,mlCounts,corrMethod){
	cutrows<-nrow(gspm)-ceiling(nrow(gspm)*0.05)
	gspm<-gspm[order(gspm$corrFactor),]
	gspm<-gspm[1:cutrows,]
	nTaggers<-nrow(gspm)
	gammdist<-fitdist(gspm$corrFactor, "gamma", method="mle")$estimate
	gammfit<-fitdist(gspm$corrFactor, "gamma", method="mle")
	meanCorrF<-qgamma(0.5,gammdist[1],gammdist[2])
	#meanCorrF<-gammdist[1]/gammdist[2]
	gdupper<-qgamma(0.975,gammdist[1],gammdist[2])
	gdlower<-qgamma(0.025,gammdist[1],gammdist[2])
	propUpper<-gdupper/meanCorrF;propLower<-gdlower/meanCorrF
	
	## Estimate number of seals per map
	numSeals<-aggregate(numTags~regionMapId,tgvutm,sum)
	numSeals$estNumSeals<-round(numSeals$numTags*meanCorrF)
	numSeals$uclNumSeals<-round(numSeals$estNumSeals*propUpper)
	numSeals$lclNumSeals<-round(numSeals$estNumSeals*propLower)
	
	maps$regionMapId<-paste0(maps$region,maps$mapId)
	mlCounts$regionMapId<-paste0(mlCounts$region,mlCounts$mapId)
	
	numSeals<-merge(numSeals,maps[,c("overlayId","regionMapId")],by="regionMapId",all.x=T)
	numSeals<-merge(numSeals,overlays[,c("overlayId","acquisition_date","satId","external_reference")],by="overlayId",all.x=T)
	numSeals$catalogId<-ifelse(grepl("amazonaws",numSeals$external_reference),substr(numSeals$external_reference,41,56),
			ifelse(grepl("ddnhehdam2vn0",numSeals$external_reference),substr(numSeals$external_reference,38,53),as.character(numSeals$external_reference)))
	#numSeals<-subset(numSeals,!is.na(acquisition_date))	#we lose 435 maps that do not have acDate, from 4500 maps we down to 4065
	numSeals$acDate<-format(numSeals$acquisition_date,"%Y%m%d")
	numSeals$acYear<-as.integer(format(numSeals$acquisition_date,"%Y"))
	numSeals$acHour<-as.integer(format(numSeals$acquisition_date,"%H"))
	numSeals$acquisition_date<-as.character(numSeals$acquisition_date)
	
	numEst<-merge(numSeals,mlCounts[,c("regionMapId","mlcount")],by="regionMapId",all.x=TRUE)
	numEst$corrMethod<-corrMethod
	
	retlist<-list(numEst=numEst,gammfit=gammfit,nTaggers=nTaggers)
	return(retlist)
}


################################### DATA
mid<-21758509	## This is Michelle LaRue's tagger Id

## load the data, prepare for analyses
load(paste0(pathToLocalGit,"compiledData.R"))

##############################
#need to filter everything to a collection of years... (all except the taggers table, right?)
overlays$year<-as.integer(format(overlays$acquisition_date,"%Y"))
aois<-merge(aois,overlays[,c("overlayId","region","year")], by=c("overlayId","region"), all.x=T)
maps<-merge(maps,overlays[,c("overlayId","region","year")], by=c("overlayId","region"), all.x=T)
views<-merge(views,maps[,c("mapId","region","year")],by=c("mapId","region"), all.x=T)
crowd<-merge(crowd,maps[,c("mapId","region","year")],by=c("mapId","region"), all.x=T)
tags<-merge(tags,views[,c("mapViewId","taggerId","region","year")],by=c("mapViewId","taggerId","region"), all.x=T)

overlays<-subset(overlays,year %in% c(2010,2011))
aois<-subset(aois,year %in% c(2010,2011))
maps<-subset(maps,year %in% c(2010,2011))
views<-subset(views,year %in% c(2010,2011))
crowd<-subset(crowd,year %in% c(2010,2011))
tags<-subset(tags,year %in% c(2010,2011))
##############################


tags<-unique(tags)
tgutm<-convertToUTM(tags,lonfield="tagcoords.x1",latfield="tagcoords.x2")

crt<-merge(crowd[,c("tagId","score","agremnt","sensor","region")],tgutm[,c("tagId","mapViewId","easting","northing","region")],by=c("tagId","region"),all.x=T)
crt<-subset(crt,!is.na(mapViewId))
crtm<-merge(crt,views[,c("mapViewId","mapId","region")],by=c("mapViewId","region"),all.x=T)
crtm<-subset(crtm,!is.na(mapId))
crtm$numFeatures<-1

tgvutm<-merge(tgutm,views[,c("mapViewId","mapId","region")],by=c("mapViewId","region"),all.x=T)
tgvutm$numTags<-1
tgvutm$regionMapId<-paste0(tgvutm$region,tgvutm$mapId)
mltdf<-subset(tgvutm,taggerId==mid)
mlCounts<-aggregate(numTags~mapId+region,data=mltdf,FUN=sum)
names(mlCounts)<-c("mapId","region","mlcount")

############################################################################
## This is probably a bad approximation...
## We want to know how many times a map was vewed:
aggviews<-aggregate(taggerId~mapId+region,data=views,NROW);names(aggviews)<-c("mapId","region","numViews")
GenProbF<-nrow(crtm)/sum(aggviews$numViews>1)	#the general probability that a feature will be found in a map = numFeatures/sum of all taggers visiting maps, for maps with at leat 2 taggers
#Calculating the overall Prob[S] and expectation for CorrFactor
mlMaps<-unique(subset(views,taggerId==mid))
crtm$regionMapId<-paste0(crtm$region,crtm$mapId)
mlMaps$regionMapId<-paste0(mlMaps$region,mlMaps$mapId)
#ProbS is the number of maps with seals (i'e', maps where MLR found a seal)/all maps visited by MLR
ProbS<-nrow(subset(crtm,regionMapId %in% mlMaps$regionMapId))/NROW(mlMaps)
expectCorrF<-ProbS/GenProbF
print(paste("Expectation for correction factor:",round(expectCorrF,2)))
############################################################################


################################################################################################################################################
###  HERE is where we split in how we estimate the numbers
###  Below we use estimates of correction factors per region individually
rankdata<-taggers[,c("taggerId","taggerScore","region")]

views$regionMapId<-paste0(views$region,views$mapId)
taggersSel<-unique(subset(views,regionMapId %in% mlMaps$regionMapId,select=c("taggerId","region","regionMapId")))
taggersSel<-subset(taggersSel,taggerId!=mid)		#taggers who share maps with ML

## What tags to use? For ProbFS...
tagsSel<-subset(tgutm,taggerId %in% taggersSel$taggerId)
tagsML<-subset(tgutm,taggerId==mid)
tagsML<-merge(tagsML,views[,c("mapId","mapViewId","region")],by=c("mapViewId","region"),all.x=T)
tagsML$regionMapId<-paste0(tagsML$region,tagsML$mapId)


## What maps can we use? For ProbF
viewsSel<-subset(views,taggerId %in% c(unique(taggersSel$taggerId),mid))

## Need the probability of a feature being found on a map and the probability of a tag being placed on a map
totalMapsInspected<-NROW(unique(views$regionMapId))
totalMapsWtags<-NROW(unique(tgvutm$regionMapId))
totalMapsWfeat<-NROW(unique(crtm$regionMapId))
probFeatInMap<-totalMapsWfeat/totalMapsInspected
probTagInMap<-totalMapsWtags/totalMapsInspected
probTagAsFeat<-totalMapsWfeat/totalMapsWtags


############ SHRINK
## ProbF = number of  tags in maps inspected/total number of maps inspected
## ProbF is the prob that any tag is found on a map.
## corrProbF should be the probability that any tagger will tag a map, which is the ratio of all tags/all maps 
gspm<-getTaggerProbabilities_tagsOnlybyTag(rankdata=rankdata,crthr=0.75,taggersSel=taggersSel,tagsSel=tagsSel,tagsML=tagsML,crtm=crtm,viewsSel=viewsSel,ProbS=ProbS)

uncorrCalc<-getMapEstimates_byTag(gspm=gspm,tgvutm=tgvutm,maps=maps,overlays=overlays,mlCounts=mlCounts,corrMethod="tagsOnlybyTag")
q<-uncorrCalc$numEst
gammfit<-uncorrCalc$gammfit;plot(gammfit)
print(uncorrCalc$nTaggers)

(sum(q$estNumSeals))
(sum(q$uclNumSeals))

## STOP HERE - review
##### NEED: to review the gamma fit
##### Also, try calculating correction factors by region.

#need to know all maps with tags, and the average number of tags per map
tagStats<-getTagStats(tags=tags,views=views,maps=maps,overlays=overlays)

## Save all results in a single data file
save(res_tagsOnlybyTag,	tagStats,compML_tag,
		file=paste0(pathToLocalGit,"continentalEstimates.RData"))


#show how #10 varies with threshold...
shrinkSeals<-data.frame()
for(crtval in c(0.7,0.75,0.8,0.85,0.9,0.95)){
	gspm<-getTaggerProbabilities_tagsOnlybyTag(rankdata=rankdata,crthr=crtval,taggersSel=taggersSel,tagsSel=tagsSel,tagsML=tagsML,crtm=crtm,viewsSel=viewsSel,ProbS=ProbS)
	res_tagsOnlybyTag<-getColonyEstimates_byTag(gspm=gspm,tgvutm=tgvutm,maps=maps,overlays=overlays,colinfo=colinfo,corrMethod="tagsOnlybyTag")
	tagsOnlybyTag<-res_tagsOnlybyTag$colEstimates
	numS<-res_tagsOnlybyTag$numSeals
	numS$threshold<-crtval
	numS<-subset(numS,mapId %in% mlMaps)
	numS<-merge(numS,mlCounts, by="mapId",all.x=T)
	numS$mlcount<-ifelse(is.na(numS$mlcount),0,numS$mlcount)
	numS<-numS[,c("estNumSeals","mlcount","threshold")]
	shrinkSeals<-rbind(shrinkSeals,numS)
}

save(res_featInMap,res_tagNearFeat,res_percTags,res_globalTagsInMap,res_tagsInMaps,res_taggerCR,res_numFeatInMap,res_trueFeatInMap,
		res_tagsOnlybyMap,res_tagsOnlybyTag,tagStats,compML_feat,compML_tag,shrinkSeals,
		file=paste0(pathToLocalGit,"/colonyEstimates.RData"))

