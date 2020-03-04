# TODO: Add comment
# 
# Author: lsalas
###############################################################################

libs<-c("raster","rgdal","plyr","dplyr","sp","rgeos","fitdistrplus")
lapply(libs, require, character.only = TRUE)

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

## Regions:
## Amundsen (AMU); East Antarctica 1 (EA1) and 2 (EA2), Queen Maud (QMA); Ross (RSS - need to join first with Pilot1 and Pilot2 data) and WAP (WAP); 
## then adding the larger Ross sea as ERE, but first joining with the rest of RSS
regions<-data.frame(region=c("Amundsen","EastAnt1","EastAnt2","QueenMaud","WAPS"),
		name=c("AMU","EA1","EA2","QMA","WAP"),stringsAsFactors=F)
pathToFiles<-"//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/ContinentalEst2/"
pathToGit<-"c:/users/lsalas/git/continentalWESEestimates/data/"

laois<-list();loverlays<-list();lmaps<-list();lviews<-list();
ltaggers<-list();ltags<-list();lcrowd<-list()
for(rr in 1:5){
	rg<-regions[rr,"region"];nm<-regions[rr,"name"]
	pth<-paste0(pathToFiles,rg,"_additional_files/",rg,".RData")
	load(pth)
	aois$region<-nm;overlays$region<-nm;maps$region<-nm;views$region<-nm;taggers$region<-nm;tags$region<-nm;crowd$region<-nm
	tags<-tags[,-6]	#removing tagDate to match Pilot1 & 2
	#must add taggerId to crowd
	crowd<-merge(crowd,tags[,c("tagId","region","taggerId")],by=c("tagId","region"),all.x=T)
	laois[[rr]]<-aois;loverlays[[rr]]<-overlays;lmaps[[rr]]<-maps;lviews[[rr]]<-views
	ltaggers[[rr]]<-taggers;ltags[[rr]]<-tags;lcrowd[[rr]]<-crowd
}

## Load the Ross sea data
pth<-paste0(pathToFiles,"Ross_additional_files/Ross.RData")
load(pth)
aois$region<-"RSS";overlays$region<-"RSS";maps$region<-"RSS";views$region<-"RSS";taggers$region<-"RSS";tags$region<-"RSS";crowd$region<-"RSS"
tags<-tags[,-6]	
crowd<-merge(crowd,tags[,c("tagId","region","taggerId")],by=c("tagId","region"),all.x=T)
aoisrs<-aois;overlaysrs<-overlays;mapsrs<-maps;viewsrs<-views;taggersrs<-taggers;tagsrs<-tags;crowdrs<-crowd

## Need to add Pilot1 and Pilot2 data for 2011
load("//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/Erebus/dataFromGeoJson3.RData")
aois$region<-"ERE";overlays$region<-"ERE";maps$region<-"ERE";views$region<-"ERE";tags$region<-"ERE";crowd$region<-"ERE"
# An important overlay is missing date:
mdt<-as.POSIXlt("2011-11-26 10:21:57.203", tz="Pacific/Auckland", fomat="%Y-%m-%d %H:%M:%S")	#Melissa gave me this date: 2011-11-25T21:21:57.203Z (UTC), which is 13 hrs behind NZDT
odf<-subset(overlays,catalogId=="103001000F5DCB00");odf$acquisition_date<-mdt
tdf<-subset(overlays,overlayId!=15762); overlays<-rbind(tdf,odf)
taggers<-read.csv("//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/Erebus/antarctica_pilot_2_20170714/user_scores.csv")
taggers$region<-"ERE"; names(taggers)<-names(taggersrs)
#must add mapId to crowd
crowd<-merge(crowd,tvmo[,c("tagId","taggerId.x","mapId")],by.x=c("tagId","taggerId"),by.y=c("tagId","taggerId.x"),all.x=T)
crowd<-crowd[,names(lcrowd[[1]])]
## There are two repeated features - different sensors, tho:
crowd$key<-paste0(crowd$mapId,"::",crowd$tagId,"::",crowd$sensor)
crowd<-subset(crowd, !key %in% c("26691::43922::QB02","26905::43952::QB02","1256997::80730::WV03"))
crowd<-crowd[,which(names(crowd)!="key")]

##filter Ross data to exclude what is already in Pilot1 and 2 - based on overlayId, not aoiId, etc.
overlaysrs<-subset(overlaysrs, !overlayId %in% unique(overlays$overlayId))
aoisrs<-subset(aoisrs, overlayId %in% unique(overlaysrs$overlayId)) 
mapsrs<-subset(mapsrs, overlayId %in% unique(overlaysrs$overlayId))
viewsrs<-subset(viewsrs,mapId %in% mapsrs$mapId)
tagsrs<-subset(tagsrs,mapViewId %in% viewsrs$mapViewId)
crowdrs<-subset(crowdrs,mapId %in% mapsrs$mapId)
taggersrs<-subset(taggersrs,!taggerId %in% taggers$taggerId)

## but the aoiIds will match existing in Pilot 1 and 2, yet referring to different data so...
aoisrs$aoiId<-paste0(aoisrs$aoiId,"A")
mapsrs$mapId<-paste0(mapsrs$mapId,"A")
viewsrs$mapViewId<-paste0(viewsrs$mapViewId,"A")
viewsrs$mapId<-paste0(viewsrs$mapId,"A")
tagsrs$tagId<-paste0(tagsrs$tagId,"A")
tagsrs$mapViewId<-paste0(tagsrs$mapViewId,"A")
crowdrs$tagId<-paste0(crowdrs$tagId,"A")
crowdrs$mapId<-paste0(crowdrs$mapId,"A")
taggersrs$taggerId<-paste0(taggersrs$taggerId,"A")

## Now combine pilot 1 & 2 with Ross survey
aois<-rbind(aois,aoisrs);overlays<-rbind(overlays,overlaysrs);maps<-rbind(maps,mapsrs);views<-rbind(views,viewsrs)
tags<-rbind(tags,tagsrs);crowd<-rbind(crowd,crowdrs);taggers<-rbind(taggers,taggersrs)

laois[["RSS"]]<-aois;loverlays[["RSS"]]<-overlays;lmaps[["RSS"]]<-maps;lviews[["RSS"]]<-views
ltaggers[["RSS"]]<-taggers;ltags[["RSS"]]<-tags;lcrowd[["RSS"]]<-crowd

aois<-do.call("rbind",laois);overlays<-do.call("rbind",loverlays);maps<-do.call("rbind",lmaps);views<-do.call("rbind",lviews)
taggers<-do.call("rbind",ltaggers);tags<-do.call("rbind",ltags);crowd<-do.call("rbind",lcrowd)

## For publication purposes, get a sense of the effort:
nrow(maps); nrow(views); nrow(overlays); nrow(tags); nrow(taggers); nrow(crowd)

##############################
#need to filter everything to a collection of years... (all except the taggers table, right?)
#mapViewIds are unique within region, but not across regions, same with maps, tags, taggers...
overlays$year<-as.integer(format(overlays$acquisition_date,"%Y"))
aois<-merge(aois,overlays[,c("overlayId","region","year")], by=c("overlayId","region"), all.x=T)
maps<-merge(maps,overlays[,c("overlayId","region","year")], by=c("overlayId","region"), all.x=T)
views<-merge(views,maps[,c("mapId","region","year")],by=c("mapId","region"), all.x=T)
crowd<-merge(crowd,maps[,c("mapId","region","year")],by=c("mapId","region"), all.x=T)
tags<-merge(tags,views[,c("mapViewId","taggerId","region","year")],by=c("mapViewId","taggerId","region"), all.x=T)

#check - it checks
nrow(maps); nrow(views); nrow(overlays); nrow(tags); nrow(taggers); nrow(crowd)

#filter by year - cuts everything in half - IMPORTANT: 4 overlays have no year, but we don't care for 3 of them (no map views for 2, and the third is from Nov8, 2015)
# Waiting for Melissa and Michelle to help with the fourth
overlays<-subset(overlays,(region=="ERE" & year==2011) | (region!="ERE" & year %in% c(2010,2011)))
aois<-subset(aois,(region=="ERE" & year==2011) | (region!="ERE" & year %in% c(2010,2011)))
maps<-subset(maps,(region=="ERE" & year==2011) | (region!="ERE" & year %in% c(2010,2011)))
views<-subset(views,(region=="ERE" & year==2011) | (region!="ERE" & year %in% c(2010,2011)))
crowd<-subset(crowd,(region=="ERE" & year==2011) | (region!="ERE" & year %in% c(2010,2011)))
tags<-subset(tags,(region=="ERE" & year==2011) | (region!="ERE" & year %in% c(2010,2011)))

#Now convert ERE back to RSS
overlays$region<-ifelse(overlays$region=="ERE","RSS",overlays$region)
aois$region<-ifelse(aois$region=="ERE","RSS",aois$region)
maps$region<-ifelse(maps$region=="ERE","RSS",maps$region)
views$region<-ifelse(views$region=="ERE","RSS",views$region)
tags$region<-ifelse(tags$region=="ERE","RSS",tags$region)
taggers$region<-ifelse(taggers$region=="ERE","RSS",taggers$region)
crowd$region<-ifelse(crowd$region=="ERE","RSS",crowd$region)

#create the unique keys by region
aois$regionAoiId<-paste0(aois$region,aois$aoiId)
maps$regionMapId<-paste0(maps$region,maps$mapId)
views$regionMapViewId<-paste0(views$region,views$mapViewId)
views$regionMapId<-paste0(views$region,views$mapId)
views$regionTaggerId<-paste0(views$region,views$taggerId)
crowd$regionMapId<-paste0(crowd$region,crowd$mapId)
crowd$regionTagId<-paste0(crowd$region,crowd$tagId)
tags$regionTagId<-paste0(tags$region,tags$tagId)
tags$regionMapViewId<-paste0(tags$region,tags$mapViewId)
tags$regionTaggerId<-paste0(tags$region,tags$taggerId)
taggers$regionTaggerId<-paste0(taggers$region,taggers$taggerId)

#####################################
## Filter for only those taggers and only those views for maps that were inspected
## maps will inform which views, which will inform which taggers, which will inform which tags
## maps will filter overlay and then aois
## tags will filter crowdranks
maps<-subset(maps, num_views>0)
views<-subset(views, regionMapId %in% maps$regionMapId)
taggers<-subset(taggers,regionTaggerId %in% unique(views$regionTaggerId))
tags<-subset(tags,(regionMapViewId %in% unique(views$regionMapViewId)) & (regionTaggerId %in% unique(taggers$regionTaggerId)))
overlays<-subset(overlays,overlayId %in% unique(maps$overlayId))
aois<-subset(aois,overlayId %in% unique(maps$overlayId))
#No need to subset crowd or taggers

#####################################
#filtering for crazy counts:
tg<-tags; tg$tagCount<-1
tg<-merge(tg,views[,c("regionMapViewId","regionMapId")],by="regionMapViewId",all.x=T)
tgcounts<-aggregate(tagCount~regionMapId+regionTaggerId+region+year+regionMapViewId,tg,sum)
sum(tgcounts$tagCount>100)  #44 cases in 40 maps
NROW(unique(subset(tgcounts,tagCount>100)$regionMapId))/nrow(tgcounts)  #47 mapviews out of 35,932, or 0.09%
## CAREFUL: there are repeated views of the same map by the same tagger, with different viewIds - here we consider them as independent
tcdf<-ldply(.data=unique(subset(tgcounts,tagCount>100)$regionMapId), .fun=function(z,tgcounts,views){
			ovis<-NROW(unique(subset(tgcounts,regionMapId==z)$regionMapViewId))
			mval<-subset(tgcounts,tagCount>100 & regionMapId==z)
			mval<-mval[,c("regionTaggerId","region","regionMapId","regionMapViewId","tagCount")]
			mval$totalViews<-ovis
			return(mval)
		},tgcounts=tgcounts,views=views)
tcdf<-tcdf[order(tcdf$regionMapId),]
print(tcdf)		#Note the many cases where there are 4+ views and only one is wrong
## We remove the erroneous counts for maps with 4+ views... (these are 36 of the 47 views!)
tcdff<-ldply(.data=unique(tcdf$regionMapId), .fun=function(m,tcdf){
			ttdf<-subset(tcdf,regionMapId==m)
			nr<-nrow(ttdf)
			nh<-unique(ttdf$totalViews)
			if(nh<(nr+3)){
				rdf<-ttdf
			}else{
				rdf<-data.frame(regionTaggerId=unique(ttdf$regionTaggerId), region=unique(ttdf$region), regionMapId=unique(ttdf$regionMapId), regionMapViewId=unique(ttdf$regionMapViewId),tagCount=NA, totalViews=NA)
			}
			
		},tcdf=tcdf)
tcdff$mapTagger<-paste0(tcdff$regionMapId,":",tcdff$regionTaggerId)
tcz<-subset(tcdff,is.na(tagCount))
tgcounts$mapTagger<-paste0(tgcounts$regionMapId,":",tgcounts$regionTaggerId)
tcz$mapTagger<-paste0(tcz$regionMapId,":",tcz$regionTaggerId)
tg$mapTagger<-paste0(tg$regionMapId,":",tg$regionTaggerId)
## now removing from tags - confirm and then remove
tg<-subset(tg,!mapTagger %in% tcz$mapTagger)
nrow(tags)-nrow(subset(tg,!mapTagger %in% tcz$mapTagger))==sum(subset(tgcounts,mapTagger %in% tcz$mapTagger)$tagCount)
tags<-subset(tags,regionTagId %in% tg$regionTagId)
nrow(tags)==nrow(tg)
## need to remove the views too, and if these are the only views of these taggers, well, then the taggers too... Luckily this is not the case, so no need to remove the taggers
views<-subset(views,!regionMapViewId %in% tcz$regionMapViewId)
#####################################
## Finding the 99% percentile value of tag counts...
gd<-fitdist(log(tgcounts$tagCount), "norm", method="mle")
est99<-ceiling(qnorm(0.99,mean=gd$estimate[1],sd=gd$estimate[2])) 	# 4 tags is the result
####################################
## Now for the remaining erroneous 5 views, use the larger value of all other views of the same map, or est99
## All views of the same map, or est99
tcdff<-subset(tcdff,!is.na(tagCount))
tcdfc<-ldply(.data=1:nrow(tcdff), .fun=function(rr,tcdff,tgcounts,est99){
			adf<-tcdff[rr,]
			rmid<-adf$regionMapId
			cdf<-subset(tgcounts,regionMapId==rmid & tagCount<100)
			if(nrow(cdf)==0){
				nval<-est99
			}else{
				nval<-ifelse(max(cdf$tagCount)>est99,max(cdf$tagCount),est99)
			}
			adf$newVal<-nval
			return(adf)
		},tcdff=tcdff,tgcounts=tgcounts,est99=est99)
tga<-subset(tg,regionMapViewId %in% tcdfc$regionMapViewId)
tgb<-subset(tg,!regionMapViewId %in% tcdfc$regionMapViewId)
tgaa<-ldply(.data=tcdfc$regionMapViewId, .fun=function(vv,tga,tcdfc){
			adf<-subset(tga,regionMapViewId==vv)
			nrw<-subset(tcdfc,regionMapViewId==vv)$newVal
			return(adf[1:nrw,])
		},tga=tga,tcdfc=tcdfc)
tgc<-rbind(tgaa,tgb)
## confirm we did it right, then apply
(nrow(tg)-nrow(tgc)) == (sum(tcdfc$tagCount)-sum(tcdfc$newVal))
tags<-subset(tags,regionTagId %in% tgc$regionTagId)
nrow(tags)==nrow(tgc)

## Removing repeated views
## There is abunch of views that are identical except for the assigned mapViewId. Removing the duplicates from the tags and views tables...
tv<-merge(tags,views[,c("regionMapId","regionMapViewId")],by="regionMapViewId",all.x=T)
tv$mapTagger<-paste0(tv$regionMapId,":",tv$regionTaggerId)
rdf<-unique(tv[,c("mapTagger","regionMapViewId")])
rvdf<-data.frame()
for(rr in unique(rdf$mapTagger)){
	tdf<-subset(rdf,mapTagger==rr)
	if(nrow(tdf)>1){
		rvdf<-rbind(rvdf,tdf[-1,])	#this has the unique views we want to delete
	}
}
tags<-subset(tags,!regionMapViewId %in% rvdf$regionMapViewId)
views<-subset(views,!regionMapViewId %in% rvdf$regionMapViewId)


#######################################
## Now converting tag coordinates to UTM.
tags<-convertToUTM(tags,lonfield="tagcoords.x1",latfield="tagcoords.x2")

save(aois,overlays,maps,views,taggers,tags,crowd,file=paste0(pathToGit,"compiledData.RData"))

