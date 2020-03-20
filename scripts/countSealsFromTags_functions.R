# TODO: Add comment
# 
# Author: lsalas
###############################################################################


######################  FUNCTIONS FOR THE SCRIPT FILE countSealsFromTags.R

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

## Workhorse function: This function estimates the corection Factor for each tagger and region
# ProbF is defined as the number of  tags in maps inspected/total number of maps inspected			
# ProbF is the prob that any tag is found on a map.
## rankdata is a data.frame listing the taggers and their crowdrank score
## crthr is the threshold to filter the taggers (i.e., only use taggers whose rank is higher than crthr)
## taggersSel is the subset of taggers inspecting maps also inspected by MLR
## tagsSel is the data.frame of tags placed by the taggers in maps also inspected by MLR
## tagsML is the data.frame of tags placed by MLR
## viewsSel is the data.frame of views of maps also inspected by MLR
## ProbS is the probability that a map has a seal: number of MLtags/number of maps inspected by ML
getTaggerProbabilities_tagsOnlybyTag<-function(rankdata,crthr,taggersSel,tagsSel,tagsML,viewsSel,ProbS){
	if(!is.na(crthr)){
		tsdf<-subset(rankdata,(region!="RSS" & taggerScore>=crthr) | (region=="RSS" & taggerScore>=0.75));	#only taggers matching or besting the threshold
		taggers<-subset(taggersSel,taggersSel$regionTaggerId %in% tsdf$regionTaggerId);	#those taggers that shared maps with Michelle and who also have a score higher or equal to threshold
	}else{
		taggers<-taggersSel
	}
	
	#calculate the probabilities for each tagger exceding the threshold crthr - need tags and feature geo
	taggerProbs<-ldply(.data=unique(taggers$regionTaggerId),.fun=function(tt,viewsSel,tagsSel,tagsML){
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
				
				
				return(ttdf)},viewsSel=viewsSel,tagsSel=tagsSel,tagsML=tagsML)
	
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


## Workhorse function - same as above, but not region-specific
getTaggerProbabilities_tagsOnlybyTag_general<-function(rankdata,crthr,taggersSel,tagsSel,tagsML,viewsSel,ProbS){
	if(!is.na(crthr)){
		tsdf<-subset(rankdata,(region!="RSS" & taggerScore>=crthr) | (region=="RSS" & taggerScore>=0.75));	#only taggers matching or besting the threshold
		taggers<-subset(taggersSel,taggersSel$taggerId %in% tsdf$taggerId);	#those taggers that shared maps with Michelle and who also have a score higher or equal to threshold
	}else{
		taggers<-taggersSel
	}
	
	#calculate the probabilities for each tagger exceding the threshold crthr - need tags and feature geo
	taggerProbs<-ldply(.data=unique(taggers$taggerId),.fun=function(tt,viewsSel,tagsSel,tagsML){
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
				
				
				return(ttdf)},viewsSel=viewsSel,tagsSel=tagsSel,tagsML=tagsML)
	
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


## Workhorse function: abbreviated, without the need to estimate SS in the confusion matrix, for each tagger and region 
#   ProbF comes from the confuson matrix: totalFeatures/(totalFeatures + totalNotFeatures), where totalNotFeatures are the seals not tagged
#   ProbS also comes from the confusion matrix: totalSeals/(totalSeals + totalNonSeals), where totalNonSeals are the features that are not seals 
## The arguments are the same as in the previous functions
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


## Workhorse function: abbreviated,  without the need to estimate SS in the confusion matrix, not region-secific
#   ProbF comes from the confuson matrix: totalFeatures/(totalFeatures + totalNotFeatures), where totalNotFeatures are the seals not tagged
#   ProbS also comes from the confusion matrix: totalSeals/(totalSeals + totalNonSeals), where totalNonSeals are the features that are not seals 
## The arguments are the same as in the previous functions
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


## This function provides estimates by map using the output from the getTaggerProbabilities_tagsOnlybyTag function or the getTaggerProbabilities_fromQ function
## cdf is the output from function getTaggerProbabilities_tagsOnlybyTag or getTaggerProbabilities_fromQ
## crthr is the threshold to filter the taggers (i.e., only use taggers whose rank is higher than crthr). See what you used for the above functions.
## taggers is the data.frame of all valid taggers participating in the survey
## nSealsFilt is a value to filter taggers based on how many seals MLR found in maps also inspected by these taggers (for Qvals only). If filtered out, their correction parameter is not used for fitting the distribution.
## tgvutm is a data.frame with properly merged tags and views, uing UTM coordinates 
## maps is the data file on maps
## overlays is the data file on overlays
## corrMethod is a string: either taggerQval or taggerCorrFactor
## dist is the distribution to fit to the values of correction factors
## regional is a boolean that informs if a general or region-specific estimate of the correction is being used
## cival indicates the cofidence limits, defaults to 95%
getMapEstimates<-function(cdf,crthr=0.5,taggers,nSealsFilt=10,tgvutm,maps,overlays,corrMethod,dist="gamma",regional=FALSE,cival=95){
	if(corrMethod=="taggerQval"){
		cdf<-subset(cdf,nSeals>=nSealsFilt)
		parn<-"Qval"
	}else{
		parn<-"corrFactor"
	}
	
	
	tgvoutm<-merge(tgvutm,maps[,c("regionMapId","overlayId")],by="regionMapId",all.x=T)
	tgvoutm<-merge(tgvoutm,overlays[,c("overlayId","satId")],by="overlayId")
	sdf<-aggregate(numTags~region+mapId+year+satId+taggerId+regionTaggerId+regionMapId+overlayId,tgvoutm,sum)
	
	#now apply Qval by tagger
	seltaggers<-taggers
	
	distvals<-fitDistLimits(dist=dist,qdf=cdf,parn=parn,cival=cival)
	
	meanCorrF<-distvals$meanCorrF
	propUpper<-distvals$vupper/meanCorrF
	propLower<-distvals$vlower/meanCorrF
	
	sdf$taggerEstNumSeals<-round(sdf$numTags*meanCorrF)
	
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
	sdfp$corrMethod<-corrMethod
	sdfp$crThreshold<-crthr
	return(sdfp)
	
}

## The following function returns the number of maps with tags given the crowdrank threshold the number of tags in those maps, and the average number of tags per map
## tags is the table with all tags and the views where these were placed
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
## cival indicates the cofidence limits, defaults to 95%
fitDistLimits<-function(dist="gamma",qdf,parn,cival=95){
	ici<-(1-(cival/100))/2
	uci<-1-ici
	if(dist=="lognorm"){
		## If using log-normal
		ndist<-fitdist(log(qdf[,parn]), "norm")$estimate
		meanCorrF<-exp(ndist[1])
		vupper<-exp(qnorm(uci,ndist[1],ndist[2]))
		vlower<-exp(qnorm(ici,ndist[1],ndist[2]))
		
	}else if(dist=="weibull"){
		## If using Weibull
		wdist<-fitdist(qdf[,parn], "weibull")$estimate
		meanCorrF<-wdist[2]*gamma(1+1/wdist[1])
		vupper<-qweibull(uci,wdist[1],wdist[2])
		vlower<-qweibull(ici,wdist[1],wdist[2])
		
	}else{	#gamma, default
		## If using Gamma
		gammdist<-fitdist(qdf[,parn], "gamma", method="mle")$estimate
		meanCorrF<-gammdist[1]/gammdist[2]
		vupper<-qgamma(uci,gammdist[1],gammdist[2])
		vlower<-qgamma(ici,gammdist[1],gammdist[2])
	}
	res<-list(meanCorrF=meanCorrF, vupper=vupper, vlower=vlower)
	return(res)
}


## The following function loads the island and colony models and predicts detection rates for a collection of maps. 
## Returns: keyField, colDetRate, islDetRate
## dat is the data.frame to predict to. Must have: scaledNumTags, scaledNumMaps, sinH
## NOTE: We replace scaledNumTags by scaledAvgTags; we use -1.391183 for scaledNumMaps, which is the equivalent of 1 map; we use acYear to be the year field
## keyFieldName is the unique record identifier. We use: "regionMapId" as default
## islandWeight is the proportion of "island" colonies in the data, defaults to 90% but likely less?
predictDetRates<-function(pathToGit,dat,keyFieldName="regionMapId",islandWeight=0.9){
	## Load the models:
	load(file=paste0(pathToGit,"data/finalModelsAndData.RData"))
	
	## Get the colony and island info df
	colIsldf<-unique(numSealsF[,c("Colony","Island")])
	colIsldf<-subset(colIsldf,Colony!="South Base")
	colIsldf$weightVal<-ifelse(colIsldf$Island==1,islandWeight,1-islandWeight)
	sumWeights<-sum(colIsldf$weightVal)
	
	## Base dataset:
	bdf<-dat[,c(keyFieldName,"scaledNumTags")]
	bdf$sinH<-0	#we already corrected for hour effect
	bdf$scaledNumMaps<-0	#assuming no effect frm number of maps, because we only have 1 image per location
	
	## The colony model is: mdlCol<-lm(detRate~scaledNumTags*scaledNumMaps+sinH+I(sinH^2)+Colony+acYear,data=numSealsF)
	## We predict assuming each map is from each of the colonies, and then average, weighing the mainland colonies at 1-islandWeight, all others at islandWeight
	
	## For 2010
	bdf$acYear<-"2010"
	cdf<-ldply(.data=colIsldf$Colony, .fun=function(x,colIsldf,bdf,mdlCol){
				colWgt<-subset(colIsldf,Colony==x)$weightVal
				tdf<-bdf
				tdf$Colony<-x
				tdf$wgtPredColRate<-predict(mdlCol,newdata=tdf)
				tdf$wgtPredColRate<-tdf$wgtPredColRate*colWgt
				return(tdf)
			},colIsldf=colIsldf,bdf=bdf,mdlCol=mdlCol)
	colMdlEst<-aggregate(as.formula(paste0("wgtPredColRate~",keyFieldName)),data=cdf,sum)
	colMdlEst$wgtPredColRate<-colMdlEst$wgtPredColRate/sumWeights
	
	## The island model is: mdlIsl<-lm(detRate~scaledNumTags+sinH+I(sinH^2)+Island+acYear,data=numSealsF)
	## We predict 100% as islands, 100% as mainland, then add up with 90% weight to island colonies.
	idf<-ldply(.data=c(0,1), .fun=function(x,bdf,mdlIsl,islandWeight){
				tdf<-bdf
				tdf$Island<-x
				islWgt<-ifelse(x==0,1-islandWeight,islandWeight)
				tdf$wgtPredIslRate<-predict(mdlIsl,newdata=tdf)
				tdf$wgtPredIslRate<-tdf$wgtPredIslRate*islWgt
				return(tdf)
			},bdf=bdf,mdlIsl=mdlIsl,islandWeight=islandWeight)
	islMdlEst<-aggregate(as.formula(paste0("wgtPredIslRate~",keyFieldName)),data=idf,sum)
	
	dfRate10<-merge(colMdlEst,islMdlEst,by=keyFieldName)
	dfRate10$Year<-"2010"
	
	## For 2011
	bdf$acYear<-"2011"
	cdf<-ldply(.data=colIsldf$Colony, .fun=function(x,colIsldf,bdf,mdlCol){
				colWgt<-subset(colIsldf,Colony==x)$weightVal
				tdf<-bdf
				tdf$Colony<-x
				tdf$wgtPredColRate<-predict(mdlCol,newdata=tdf)
				tdf$wgtPredColRate<-tdf$wgtPredColRate*colWgt
				return(tdf)
			},colIsldf=colIsldf,bdf=bdf,mdlCol=mdlCol)
	colMdlEst<-aggregate(as.formula(paste0("wgtPredColRate~",keyFieldName)),data=cdf,sum)
	colMdlEst$wgtPredColRate<-colMdlEst$wgtPredColRate/sumWeights
	
	## The island model is: mdlIsl<-lm(detRate~scaledNumTags+sinH+I(sinH^2)+Island+acYear,data=numSealsF)
	## We predict 100% as islands, 100% as mainland, then add up with 90% weight to island colonies.
	idf<-ldply(.data=c(0,1), .fun=function(x,bdf,mdlIsl,islandWeight){
				tdf<-bdf
				tdf$Island<-x
				islWgt<-ifelse(x==0,1-islandWeight,islandWeight)
				tdf$wgtPredIslRate<-predict(mdlIsl,newdata=tdf)
				tdf$wgtPredIslRate<-tdf$wgtPredIslRate*islWgt
				return(tdf)
			},bdf=bdf,mdlIsl=mdlIsl,islandWeight=islandWeight)
	islMdlEst<-aggregate(as.formula(paste0("wgtPredIslRate~",keyFieldName)),data=idf,sum)
	
	dfRate11<-merge(colMdlEst,islMdlEst,by=keyFieldName)
	dfRate11$Year<-"2011"
	dfRate<-rbind(dfRate10,dfRate11)
	return(dfRate)
	
}


## The following function attributes the WESE counts with the cell ID of a 5-km grid cell for geospatial analyses of habitat selection
## It then adds up the counts of seals per 5-km grid cell and returns the resulting table
## gdf is the spatialpoints data.frame
## wesedf is the data.frame of counts of seals by map
getWESEcountsBy5km<-function(gdf,wesedf){
	pdf<-as.data.frame(gdf)
	wcoords<-wesedf[,c("regionMapId","mapcoords.x1","mapcoords.x2")]
	coordinates(wcoords)<-c("mapcoords.x1","mapcoords.x2")
	proj4string(wcoords)<-CRS("+proj=longlat +datum=WGS84")
	wdf<-spTransform(wcoords,CRS(projection(gdf)))
	wdf<-as.data.frame(wdf)
	
	#Find the closest gpoint to each wpoint
	wgdf<-unique(wdf[,c("regionMapId","mapcoords.x1","mapcoords.x2")])
	wgdf$pointid<-NA
	for(m in 1:nrow(wgdf)){
		mc1<-wgdf[m,"mapcoords.x1"];mc2<-wgdf[m,"mapcoords.x2"]
		pdf$dist<-sqrt(((pdf$coords.x1-mc1)^2)+((pdf$coords.x2-mc2)^2))
		pdf<-pdf[order(pdf$dist),]
		ncId<-pdf[1,"pointid"]
		wgdf[m,"pointid"]<-ncId
	}
	return(wgdf)
}


