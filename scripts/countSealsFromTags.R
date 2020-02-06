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

libs<-c("ggplot2","raster","rgdal","plyr","dplyr","sp","rgeos","fitdistrplus")
lapply(libs, require, character.only = TRUE)

pathToLocalGit<-"C:/Users/lsalas/git/ContinentalWESEestimates/data/"

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

## Find the number of seals the tager correctly identified (compared to ML)
# taggerTags is the data.frame with all the tagger's tags in maps shared with ML, including the east/north of each tag
# mlTags is the data.frame with all of ML tags, including the east/north of each tag
getTaggerSS<-function(taggerTags,mlTags){
	#loop through each regiontag and see if it's within 3m of a ML regiontag
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

## Workhorse function: ProbF = number of  tags in maps inspected/total number of maps inspected			
## ProbF is the prob that any tag is found on a map.
## corrProbF should be the probability that any tagger will tag a map, which is the ratio of all tags/all maps
getTaggerProbabilities_tagsOnlybyTag<-function(rankdata,crthr,taggersSel,tagsSel,tagsML,crtm,viewsSel,ProbS,corrProbF=1){
	if(!is.na(crthr)){
		tsdf<-subset(rankdata,(region!="RSS" & taggerScore>=crthr) | (region=="RSS" & taggerScore>=0.75));	#only taggers matching or besting the threshold
		taggers<-subset(taggersSel,taggersSel$regionTaggerId %in% tsdf$regionTaggerId);	#those taggers that shared maps with Michelle and who also have a score higher or equal to threshold
	}else{
		taggers<-taggersSel
	}
	
	#calculate the probabilities for each tagger exceding the threshold crthr - need tags and feature geo
	taggerProbs<-ldply(.data=unique(taggers$regionTaggerId),.fun=function(tt,viewsSel,crtm,tagsSel,tagsML){
				#calculating ProbF...(number of maps with tags among the maps inspected by this tagger)
				#need the number of maps inspected by this tagger
				taggerMaps<-unique(subset(viewsSel,regionTaggerId==tt)$regionMapId)
				region<-substr(tt,1,3)
				tagger_nMaps<-NROW(taggerMaps)
				#determine how many of these have tags
				taggerTags<-subset(tagsSel,regionTaggerId==tt)
				taggerMapsWTags<-NROW(unique(taggerTags$regionMapId))
				totalTags<-nrow(taggerTags)
				
				#how many seals in the maps inspected by the tagger
				mlTags<-subset(tagsML,regionMapId %in% taggerMaps)	#if this is 0, it's because ML did not find a seal in any of the maps this tagger inspected.
				totalSeals<-nrow(mlTags)
				
				ProbF<-totalTags/tagger_nMaps 	#Use total number of features/total number of maps inspected??
				#ProbF<-taggerMapsWTags/tagger_nMaps	##DO NOT...use the number of maps with features/total number of maps inspected
				
				if(totalTags>0){	#features were found in maps inspected by the tagger
					#calculating ProbFS...
					#need the list of maps from this tagger also shared with ML - this is taggerMaps
					#need the list of tags in these maps
					
					if(nrow(mlTags)==0){	#but no seals...
						ttdf<-data.frame(regionTaggerId=tt,ss=0,nTags=totalTags,nSeals=totalSeals,nMaps=tagger_nMaps,ProbFS=0,ProbF=0,region=region)
					}else{		#and seals
						ss<-getTaggerSS(taggerTags,mlTags)
						ProbFS<-ss/totalSeals
						
						ttdf<-data.frame(regionTaggerId=tt,ss=ss,nTags=totalTags,nSeals=totalSeals,nMaps=tagger_nMaps,ProbFS=ProbFS,ProbF=ProbF,region=region)
					}
				}else{	#no features found in maps inspected by this tagger
					ttdf<-data.frame(regionTaggerId=tt,ss=0,nTags=totalTags,nSeals=0,nMaps=tagger_nMaps,ProbFS=0,ProbF=0,region=region)
				}
				
				
				return(ttdf)},viewsSel=viewsSel,crtm=crtm,tagsSel=tagsSel,tagsML=tagsML)
	
	taggerProbs<-subset(taggerProbs,ProbFS>0)	#use only data from taggers that had seals and features...
	if(nrow(taggerProbs)>0){
		taggerProbs$ProbS<-ProbS
		taggerProbs$ProbSF<-taggerProbs$ProbS*taggerProbs$ProbFS/taggerProbs$ProbF	
		taggerProbs$corrFactor<-taggerProbs$ProbSF/taggerProbs$ProbFS
	}else{ #none of the computed probabilities permits the estimation of corrFactor
		taggerProbs<-NA
	}
	
	return(taggerProbs)
}


## Workhorse function - same as above, but not region-specific: 
## ProbF = number of  tags in maps inspected/total number of maps inspected			
## ProbF is the prob that any tag is found on a map.
## corrProbF should be the probability that any tagger will tag a map, which is the ratio of all tags/all maps
getTaggerProbabilities_tagsOnlybyTag_general<-function(rankdata,crthr,taggersSel,tagsSel,tagsML,crtm,viewsSel,ProbS,corrProbF=1){
	if(!is.na(crthr)){
		tsdf<-subset(rankdata,(region!="RSS" & taggerScore>=crthr) | (region=="RSS" & taggerScore>=0.75));	#only taggers matching or besting the threshold
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
				taggerTags<-subset(tagsSel,taggerId==tt)
				taggerMapsWTags<-NROW(unique(taggerTags$regionMapId))
				totalTags<-nrow(taggerTags)
				
				#how many seals in the maps inspected by the tagger
				mlTags<-subset(tagsML,regionMapId %in% taggerMaps)	#if this is 0, it's because ML did not find a seal in any of the maps this tagger inspected.
				totalSeals<-nrow(mlTags)
				
				ProbF<-totalTags/tagger_nMaps 	#Use total number of features/total number of maps inspected??
				#ProbF<-taggerMapsWTags/tagger_nMaps	##DO NOT...use the number of maps with features/total number of maps inspected
				
				if(totalTags>0){	#features were found in maps inspected by the tagger
					#calculating ProbFS...
					#need the list of maps from this tagger also shared with ML - this is taggerMaps
					#need the list of tags in these maps
					
					if(nrow(mlTags)==0){	#but no seals...
						ttdf<-data.frame(taggerId=tt,ss=0,nTags=totalTags,nSeals=totalSeals,nMaps=tagger_nMaps,ProbFS=0,ProbF=0)
					}else{		#and seals
						ss<-getTaggerSS(taggerTags,mlTags)
						ProbFS<-ss/totalSeals
						
						ttdf<-data.frame(taggerId=tt,ss=ss,nTags=totalTags,nSeals=totalSeals,nMaps=tagger_nMaps,ProbFS=ProbFS,ProbF=ProbF)
					}
				}else{	#no features found in maps inspected by this tagger
					ttdf<-data.frame(taggerId=tt,ss=0,nTags=totalTags,nSeals=0,nMaps=tagger_nMaps,ProbFS=0,ProbF=0)
				}
				
				
				return(ttdf)},viewsSel=viewsSel,crtm=crtm,tagsSel=tagsSel,tagsML=tagsML)
	
	taggerProbs<-subset(taggerProbs,ProbFS>0)	#use only data from taggers that had seals and features...
	if(nrow(taggerProbs)>0){
		taggerProbs$ProbS<-ProbS
		taggerProbs$ProbSF<-taggerProbs$ProbS*taggerProbs$ProbFS/taggerProbs$ProbF	
		taggerProbs$corrFactor<-taggerProbs$ProbSF/taggerProbs$ProbFS
	}else{ #none of the computed probabilities permits the estimation of corrFactor
		taggerProbs<-NA
	}
	
	return(taggerProbs)
}


## Workhorse function: abbreviated, entirely confusion-matrix-based 
#   ProbF comes from the confuson matrix: totalFeatures/(totalFeatures + totalNotFeatures), where totalNotFeatures are the seals not tagged
#   ProbS also comes from the confusion matrix: totalSeals/(totalSeals + totalNonSeals), where totalNonSeals are the features that are not seals 
getTaggerProbabilities_fromQ<-function(rankdata,crthr,taggersSel,tagsSel,tagsML,mapsML,viewsSel){
	if(!is.na(crthr)){
		tsdf<-subset(rankdata,(region!="RSS" & taggerScore>=crthr) | (region=="RSS" & taggerScore>=0.75));	#only taggers matching or besting the threshold
		taggers<-subset(taggersSel,taggersSel$regionTaggerId %in% tsdf$regionTaggerId);	#those taggers that shared maps with Michelle and who also have a score higher or equal to threshold
		###REVIEW!!
	}else{
		taggers<-taggersSel
	}
	
	#calculate the quotient ProbS/ProbF for each tagger exceding the threshold crthr - need tags and feature geo
	taggerProbs<-ldply(.data=unique(taggers$regionTaggerId),.fun=function(tt,viewsSel,tagsSel,tagsML,mapsML){
				#calculating ProbF...(number of maps with tags among the maps inspected by this tagger)
				#need the number of maps inspected by this tagger
				taggerMaps<-unique(subset(viewsSel,regionTaggerId==tt)$regionMapId)
				region<-substr(tt,1,3)
				#determine how many of these overlap with MLRs and how many seals and features are we dealing with
				mlMaps<-subset(mapsML,regionMapId %in% taggerMaps); numOverMaps<-nrow(mlMaps)
				mlTags<-subset(tagsML,regionMapId %in% taggerMaps); TS<-nrow(mlTags)
				taggerTags<-subset(tagsSel,regionTaggerId==tt);TF<-nrow(taggerTags)
				# Need TS, TF and TT, where TS is total seals (or SS + SN), TF is total features (or SS + NS), and TT is total tags combined (or SS + SN + NS)
				# SS is seals tagged, NS is features that are not seals, and SN are seals not tagged
				# But since Qval = ProbS/ProbF, and ProbS = TS/TT and ProbF = TF/TT, then Qval=TS/TF
				if(TF>0){
					Qval=TS/TF
				}else{	#no features found in maps inspected by this tagger
					Qval<-NA
				}
				
				ttdf<-data.frame(regionTaggerId=tt,nTags=TF,nSeals=TS,nMaps=numOverMaps,Qval=Qval,region=region)
				return(ttdf)
			},viewsSel=viewsSel,tagsSel=tagsSel,tagsML=tagsML,mapsML=mapsML)
	
	taggerProbs<-subset(taggerProbs,!is.na(Qval) & nSeals>0)
	return(taggerProbs)
}


## Workhorse function: abbreviated, entirely confusion-matrix-based 
#   ProbF comes from the confuson matrix: totalFeatures/(totalFeatures + totalNotFeatures), where totalNotFeatures are the seals not tagged
#   ProbS also comes from the confusion matrix: totalSeals/(totalSeals + totalNonSeals), where totalNonSeals are the features that are not seals 
getTaggerProbabilities_fromQ_general<-function(rankdata,crthr,taggersSel,tagsSel,tagsML,mapsML,viewsSel){
	if(!is.na(crthr)){
		tsdf<-subset(rankdata,(region!="RSS" & taggerScore>=crthr) | (region=="RSS" & taggerScore>=0.75));	#only taggers matching or besting the threshold
		taggers<-subset(taggersSel,taggersSel$regionTaggerId %in% tsdf$regionTaggerId);	#those taggers that shared maps with Michelle and who also have a score higher or equal to threshold
		###REVIEW!!
	}else{
		taggers<-taggersSel
	}
	
	#calculate the quotient ProbS/ProbF for each tagger exceding the threshold crthr - need tags and feature geo
	taggerProbs<-ldply(.data=unique(taggers$taggerId),.fun=function(tt,viewsSel,tagsSel,tagsML,mapsML){
				#calculating ProbF...(number of maps with tags among the maps inspected by this tagger)
				#need the number of maps inspected by this tagger
				taggerMaps<-unique(subset(viewsSel,taggerId==tt)$regionMapId)
				#determine how many of these overlap with MLRs and how many seals and features are we dealing with
				mlMaps<-subset(mapsML,regionMapId %in% taggerMaps); numOverMaps<-nrow(mlMaps)
				mlTags<-subset(tagsML,regionMapId %in% taggerMaps); TS<-nrow(mlTags)
				taggerTags<-subset(tagsSel,taggerId==tt);TF<-nrow(taggerTags)
				# Need TS, TF and TT, where TS is total seals (or SS + SN), TF is total features (or SS + NS), and TT is total tags combined (or SS + SN + NS)
				# SS is seals tagged, NS is features that are not seals, and SN are seals not tagged
				# But since Qval = ProbS/ProbF, and ProbS = TS/TT and ProbF = TF/TT, then Qval=TS/TF
				if(TF>0){
					Qval=TS/TF
				}else{	#no features found in maps inspected by this tagger
					Qval<-NA
				}
				
				ttdf<-data.frame(taggerId=tt,nTags=TF,nSeals=TS,nMaps=numOverMaps,Qval=Qval)
				return(ttdf)
			},viewsSel=viewsSel,tagsSel=tagsSel,tagsML=tagsML,mapsML=mapsML)
	
	taggerProbs<-subset(taggerProbs,!is.na(Qval) & nSeals>0)
	return(taggerProbs)
}


## Estimates by map
## gspm is the output from getTaggerProbabilities_tagsOnlybyTag
## tgvutm is a data.frame with properly merged tags and views, uing UTM coordinates 
## maps is the data file on maps
## overlays is the data file on overlays
## mlCounts is teh data.frame of MLR's counts by mapId and region
## corrMethod is a string naming how the counts were corrected
## dist is the distribution to fit to the values of correction factors
getMapEstimates_byTag<-function(gspm,tgvutm,maps,overlays,mlCounts,dist="gamma"){
	#cutrows<-nrow(gspm)-ceiling(nrow(gspm)*0.05)
	#gspm<-gspm[order(gspm$corrFactor),]
	#gspm<-gspm[1:cutrows,]
	nTaggers<-nrow(gspm)
	
	distvals<-fitDistLimits(dist=dist,qdf=gspm,parn="corrFactor")
	
	meanCorrF<-distvals$meanCorrF
	propUpper<-distvals$vupper/meanCorrF
	propLower<-distvals$vlower/meanCorrF
	
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
	numEst$corrMethod<-"taggerSSCorr"
	numEst$region<-substr(numEst$regionMapId,1,3)
	
	return(numEst)
}

## This function is similar to the one above, but using Qval instead.
getMapEstimates_fromQval<-function(qdf,crthr,taggers,nSealsFilt=10,tgvutm,maps,overlays,dist="gamma"){
	
	tgvoutm<-merge(tgvutm,maps[,c("regionMapId","overlayId")],by="regionMapId",all.x=T)
	tgvoutm<-merge(tgvoutm,overlays[,c("overlayId","satId")],by="overlayId")
	sdf<-aggregate(numTags~region+mapId+year+satId+taggerId+regionTaggerId+regionMapId+overlayId,tgvoutm,sum)
	
	#now apply Qval by tagger
	#seltaggers<-subset(taggers,taggerScore>=crthr)
	seltaggers<-taggers
	qdf<-subset(qdf,nSeals>=nSealsFilt)
	#qdf<-subset(qdf,regionTaggerId %in% seltaggers$regionTaggerId) #Don't need this - the tagger data were already pre-filtered
	#fit the gamma for qvals
	#get the median, lower 2.5 and upper 97.5%
	nTaggers<-nrow(qdf)
	
	distvals<-fitDistLimits(dist=dist,qdf=qdf,parn="Qval")
	
	meanCorrF<-distvals$meanCorrF
	propUpper<-distvals$vupper/meanCorrF
	propLower<-distvals$vlower/meanCorrF
	
	sdf$taggerEstNumSeals<-ceiling(sdf$numTags*meanCorrF)
	sdf$taggerUclNumSeals<-ceiling(sdf$taggerEstNumSeals*propUpper)
	sdf$taggerLclNumSeals<-ceiling(sdf$taggerEstNumSeals*propLower)
	sdf<-merge(sdf,seltaggers[,c("regionTaggerId","taggerScore")],by="regionTaggerId",all.x=TRUE)
	
	#for each map in sdf
		#average the estimated value across all taggers who inspected it
		#also calculate the weighted average using crscores
	sdfp <- sdf %>% 
			group_by(region,mapId,regionMapId,year,overlayId,satId) %>% 
			dplyr::summarise(numTaggers=NROW(unique(regionTaggerId)),estNumSeals=round(mean(taggerEstNumSeals)),uclNumSeals=round(mean(taggerUclNumSeals)),lclNumSeals=round(mean(taggerLclNumSeals)),
					wgtEstNumSeals=round(weighted.mean(taggerEstNumSeals,taggerScore)),wgtUclNumSeals=round(weighted.mean(taggerUclNumSeals,taggerScore)),wgtLclNumSeals=round(weighted.mean(taggerLclNumSeals,taggerScore)))
	sdfp<-as.data.frame(sdfp)
	sdfp$corrMethod<-"taggerQVal"
	sdfp$crThreshold<-crthr
	return(sdfp)
	
}

## The following function returns the number of maps with tags given the crowdrank threshold
## the number of tags in those maps, and the average number of tags per map
# tags is the table with all tags and the views where these were placed
## views is the table of all map views
## maps is the table with information about each map and its link to an image (overlay)
## overlays is the table with data on the satellite images 
## byRegion is a binary that summarizes the countsby region 
getTagStats<-function(tags,views,maps,overlays,byRegion=TRUE){
	#need to know all maps with tags, and the average number of tags per map
	tgv<-merge(tags,views[,c("mapViewId","mapId","region","regionMapId")],by=c("mapViewId","region"),all.x=TRUE)
	tgmp<-merge(tgv,maps[,c("mapId","region","overlayId","mapcoords.x1","mapcoords.x2")],by=c("mapId","region"),all.x=TRUE)
	tgco<-merge(tgmp,overlays[,c("overlayId","acquisition_date")],by="overlayId",all.x=TRUE)
	tgco$satdate<-format(tgco$acquisition_date,"%Y-%m-%d %H:%M:%S")
	tagsPerMapDate<-aggregate(mapId~regionMapId+satdate+region+mapcoords.x1+mapcoords.x2,data=tgco,NROW)
	names(tagsPerMapDate)<-c("RegionMapId","ImgDate","Region","mapLon","mapLat","NumTags")
	if(byRegion){
		res<-aggregate(NumTags~ImgDate+Region,data=tagsPerMapDate,sum)
	}else{
		res<-tagsPerMapDate
	}
	return(tagsPerMapDate)
}


## The following function fits a distribution to the data of correction factors
## dist is the named distribution to fit, either: "gamma", "lognorm"", or "weibull"
## qdf is the data.frame of correction factor values
## parn is the name of the correction factor parameter in qdf
fitDistLimits<-function(dist="gamma",qdf,parn){
	if(dist=="lognorm"){
		## If using log-normal
		ndist<-fitdist(log(qdf[,parn]), "norm")$estimate
		meanCorrF<-exp(ndist[1])
		vupper<-exp(qnorm(0.975,ndist[1],ndist[2]))
		vlower<-exp(qnorm(0.025,ndist[1],ndist[2]))
		
	}else if(dist=="weibull"){
		## If using Weibull
		wdist<-fitdist(qdf[,parn], "weibull")$estimate
		meanCorrF<-wdist[2]*gamma(1+1/wdist[1])
		vupper<-qweibull(0.975,wdist[1],wdist[2])
		vlower<-qweibull(0.025,wdist[1],wdist[2])
		
	}else{	#gamma, default
		## If using Gamma
		gammdist<-fitdist(qdf[,parn], "gamma", method="mle")$estimate
		meanCorrF<-gammdist[1]/gammdist[2]
		vupper<-qgamma(0.975,gammdist[1],gammdist[2])
		vlower<-qgamma(0.025,gammdist[1],gammdist[2])
	}
	res<-list(meanCorrF=meanCorrF, vupper=vupper, vlower=vlower)
	return(res)
}

################################### DATA
mid<-21758509	## This is Michelle LaRue's tagger Id

## load the data, prepare for analyses
load(paste0(pathToLocalGit,"compiledData.RData"))


##############################
#prepare to count tags and convert to UTM
tags<-unique(tags)
tgutm<-convertToUTM(tags,lonfield="tagcoords.x1",latfield="tagcoords.x2")
tgvutm<-merge(tgutm,views[,c("mapViewId","mapId","region")],by=c("mapViewId","region"),all.x=T)
tgvutm$numTags<-1
tgvutm$regionMapId<-paste0(tgvutm$region,tgvutm$mapId)

#prepare the feature data
crt<-merge(crowd[,c("tagId","score","agremnt","sensor","region")],tgutm[,c("tagId","mapViewId","easting","northing","region","regionTagId","regionMapViewId","regionTaggerId")],by=c("tagId","region"),all.x=T)
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
#That is, roughtly we count 1.76 seals for every tag placed
## THIS IS A PROBLEM


############ SHRINK
## ProbF (for each tagger separately) = number of  tags in maps inspected/total number of maps inspected 
## ProbF is the prob that any tag is found on a map.
## corrProbF should be the probability that any tagger will tag a map, which is the ratio of all tags/all maps 
## Version 1A: region-specific estimates
gspm<-getTaggerProbabilities_tagsOnlybyTag(rankdata=rankdata,crthr=0.5,taggersSel=subtaggers,tagsSel=tagsSel,tagsML=tagsML,crtm=crtm,viewsSel=viewsSel,ProbS=estProbS)

## Version 1B: general
gspmG<-getTaggerProbabilities_tagsOnlybyTag_general(rankdata=rankdata,crthr=0.5,taggersSel=subtaggers,tagsSel=tagsSel,tagsML=tagsML,crtm=crtm,viewsSel=viewsSel,ProbS=estProbS)

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
countByQ<-getMapEstimates_fromQval(qdf=gspmQ,crthr=0.5,taggers=subtaggers,nSealsFilt=8,tgvutm=tgvutm,maps=maps,overlays=overlays,dist="weibull")
estByRegionQ<-as.data.frame(countByQ %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)),estNumSeals=ceiling(sum(estNumSeals)),uclNumSeals=ceiling(sum(uclNumSeals))))
estByRegionQ<-rbind(estByRegionQ,data.frame(region="Total",lclNumSeals=round(sum(countByQ$lclNumSeals)),estNumSeals=round(sum(countByQ$estNumSeals)),uclNumSeals=round(sum(countByQ$uclNumSeals))))
print(estByRegionQ)

## If we were to use the general form:
countByQG<-getMapEstimates_fromQval(qdf=gspmQG,crthr=0.5,taggers=subtaggers,nSealsFilt=8,tgvutm=tgvutm,maps=maps,overlays=overlays,dist="weibull")
estByRegionQG<-as.data.frame(countByQG %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)),estNumSeals=ceiling(sum(estNumSeals)),uclNumSeals=ceiling(sum(uclNumSeals))))
estByRegionQG<-rbind(estByRegionQG,data.frame(region="Total",lclNumSeals=round(sum(countByQG$lclNumSeals)),estNumSeals=round(sum(countByQG$estNumSeals)),uclNumSeals=round(sum(countByQG$uclNumSeals))))
print(estByRegionQG)

## And if we used the tagger-specific ss values:
countBySS<-getMapEstimates_byTag(gspm=gspm,tgvutm=tgvutm,maps=maps,overlays=overlays,mlCounts=mlCounts)
estByRegionSS<-as.data.frame(countBySS %>% group_by(region) %>% dplyr::summarize(lclNumSeals=ceiling(sum(lclNumSeals)),estNumSeals=ceiling(sum(estNumSeals)),uclNumSeals=ceiling(sum(uclNumSeals))))
estByRegionSS<-rbind(estByRegionSS,data.frame(region="Total",lclNumSeals=round(sum(countBySS$lclNumSeals)),estNumSeals=round(sum(countBySS$estNumSeals)),uclNumSeals=round(sum(countBySS$uclNumSeals))))
print(estByRegionSS)

## However, these numbers are probably too high - according to this graph, x5 times too high
ggplot(subset(countBySS,!is.na(mlcount)),aes(x=mlcount,y=estNumSeals)) + geom_point() + geom_abline(slope=1,intercept=0, color="blue") + stat_smooth(method="lm",se=F, formula=y~x,color="red")
## CANNOT WORK WITH the ss estimates because the sample of surveyors is too biased!

save(estByRegionQ, estByRegionQG, file=paste0(pathToLocalGit,"estimatesByMap.RData"))

