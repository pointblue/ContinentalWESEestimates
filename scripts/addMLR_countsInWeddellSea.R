# TODO: Add comment
# 
# Author: lsalas
###############################################################################


libs<-c("raster","rgdal","plyr","dplyr","sp","rgeos","ggplot2")
lapply(libs, require, character.only = TRUE)

pathToGit<-"c:/users/lsalas/git/continentalWESEestimates/data/"

## We'll need this function:
findWithin500<-function(tt,dat){
	tgdf<-subset(dat, regionTagId!=tt,select=c("regionTagId","easting","northing"))
	ct<-subset(dat, regionTagId==tt,select=c("easting","northing"))
	tgdf$tge<-ct$easting;tgdf$tgn<-ct$northing
	tgdf$dist<-sqrt(((tgdf$easting-tgdf$tge)^2) + ((tgdf$northing-tgdf$tgn)^2))
	ttgdf<-subset(tgdf,dist<500)
	rdat<-subset(dat,regionTagId==tt)
	if(nrow(ttgdf)>0){
		rrdat<-subset(dat,regionTagId %in% ttgdf$regionTagId)
		rdat<-rbind(rdat,rrdat)
	}
	return(rdat)
}

## We want to organize the data in "maps"
## The resulting table must have: "estNumSeals","region","satId","scaledTotalTags","numViews","numTaggers","originHour")]
## So I need the mean and sd for scaling the total tags: 
ttmean<-8.91902114734369; ttsd=23.33190853818
## I also need the sensorId, which we'll assume is WV02

## Load the data
wesedat<-readOGR("//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/WeddellSea", layer="WS_MLcheck_seals")
wesedf<-wesedat[,c("Lat","Lon","Image_Aqui")]; wesedf<-as.data.frame(wesedf)
wesedf$acquisition_date<-as.POSIXlt(as.character(wesedf$Image_Aqui), format="%a %b %d %Y %T", tz="MST")
## convert to utm and use the latlon to cluster to generate the maps
## then aggregate based on these "maps" and construct the data.frame we need

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

wesedf<-convertToUTM(data=wesedf,latfield="Lat",lonfield="Lon")
wesedf$tagId<-1:nrow(wesedf)
wesedf$regionTagId<-paste0("WED",wesedf$tagId)

## Let's see how this looks like on the map
ggplot(wesedf,aes(x=easting,y=northing)) + geom_point()
ggplot(subset(wesedf,easting<500000),aes(x=easting,y=northing)) + geom_point()
ggplot(subset(wesedf,easting<250000),aes(x=easting,y=northing)) + geom_point()

## First step: take every seal and assign to its own mapId if it is the only animal within 500m
maps<-data.frame();mm<-0; dat<-wesedf
for(tt in dat$regionTagId){
	if(nrow(subset(dat,regionTagId==tt))>0){
		tdf<-findWithin500(tt,dat)
		#if(nrow(tdf)==1){
		mm<-mm+1
		tdf$regionMapId<-paste0("WED",mm)
		maps<-rbind(maps,tdf)
		dat<-subset(dat,!regionTagId %in% maps$regionTagId)
	#}
	}
	
}

wesedf<-merge(wesedf,maps[,c("regionTagId","regionMapId")],by="regionTagId",all.x=T)
sum(is.na(wesedf$regionMapId))
wesedf$tagCount<-1

wdf<-wesedf[,c("regionMapId","tagCount","coords.x1","coords.x2")]

## Now we aggregate the estimated number of seals per map and calculate the average elat lon for the map too.
weddellMaps<-as.data.frame(wdf %>% group_by(regionMapId) %>% 
				dplyr::summarize(lclNumSeals=sum(tagCount),estNumSeals=sum(tagCount),uclNumSeals=sum(tagCount),mapcoords.x1=mean(coords.x1),mapcoords.x2=mean(coords.x2)))
weddellMaps<-merge(weddellMaps,unique(wesedf[,c("regionMapId","acquisition_date")]),by="regionMapId",all.x=T)

weddellMaps$region<-"WAP"
weddellMaps$satId<-"WV02"
weddellMaps$numViews<-1
weddellMaps$numTaggers<-1
weddellMaps$originHour<-format(weddellMaps$acquisition_date,"%H")
weddellMaps$corrMethod<-"MLRcount"
weddellMaps$crThreshold<-NA
weddellMaps$mlCount<-weddellMaps$estNumSeals
weddellMaps$originTZ<-format(weddellMaps$acquisition_date,"%Z")
weddellMaps$year<-format(weddellMaps$acquisition_date,"%Y")
weddellMaps$totalTags<-weddellMaps$estNumSeals
weddellMaps$scaledTotalTags<-(weddellMaps$totalTags-ttmean)/ttsd

save(weddellMaps,file=paste0(pathToGit,"MLR_WeddellSea_counts.RData"))
##'totalTags''scaledTotalTags'
