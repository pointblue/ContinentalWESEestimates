# TODO: Add comment
# 
# Author: lsalas
###############################################################################

## Regions: (loop)
## Amundsen (AMU); East Antarctica 1 (EA1) and 2 (EA2), Queen Maud (QMA); Ross (RSS) and WAP (WAP)
regions<-data.frame(region=c("Amundsen","EastAnt1","EastAnt2","QueenMaud","Ross","WAPS"),
		name=c("AMU","EA1","EA2","QMA","RSS","WAP"),stringsAsFactors=F)
pathToFiles<-"//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/ContinentalEst2/"
pathToGit<-"c:/users/lsalas/git/continentalWESEestimates/data/"

laois<-list();loverlays<-list();lmaps<-list();lviews<-list();
ltaggers<-list();ltags<-list();lcrowd<-list()
for(rr in 1:6){
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
## Need to add Pilot1 and Pilot2 data for 2011
load("//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/Erebus/dataFromGeoJson3.RData")
aois$region<-"ERE";overlays$region<-"ERE";maps$region<-"ERE";views$region<-"ERE";taggers$region<-"ERE";tags$region<-"ERE";crowd$region<-"ERE"
#must add mapId to crowd
crowd<-merge(crowd,tvmo[,c("tagId","taggerId.x","mapId")],by.x=c("tagId","taggerId"),by.y=c("tagId","taggerId.x"),all.x=T)

crowd<-crowd[,names(lcrowd[[1]])]


laois[["ERE"]]<-aois;loverlays[["ERE"]]<-overlays;lmaps[["ERE"]]<-maps;lviews[["ERE"]]<-views
ltaggers[["ERE"]]<-taggers;ltags[["ERE"]]<-tags;lcrowd[["ERE"]]<-crowd


aois<-do.call("rbind",laois);overlays<-do.call("rbind",loverlays);maps<-do.call("rbind",lmaps);views<-do.call("rbind",lviews)
taggers<-do.call("rbind",ltaggers);tags<-do.call("rbind",ltags);crowd<-do.call("rbind",lcrowd)

##############################
#need to filter everything to a collection of years... (all except the taggers table, right?)
#mapViewIds are unique within region, but not across regions, same with maps, tags, taggers...
overlays$year<-as.integer(format(overlays$acquisition_date,"%Y"))
aois<-merge(aois,overlays[,c("overlayId","region","year")], by=c("overlayId","region"), all.x=T)
maps<-merge(maps,overlays[,c("overlayId","region","year")], by=c("overlayId","region"), all.x=T)
views<-merge(views,maps[,c("mapId","region","year")],by=c("mapId","region"), all.x=T)
crowd<-merge(crowd,maps[,c("mapId","region","year")],by=c("mapId","region"), all.x=T)
tags<-merge(tags,views[,c("mapViewId","taggerId","region","year")],by=c("mapViewId","taggerId","region"), all.x=T)

#filter by year - cuts everything in half
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
taggers<-subset(taggers, regionTaggerId %in% unique(views$regionTaggerId))
tags<-subset(tags,regionTaggerId %in% unique(taggers$regionTaggerId))
crowd<-subset(crowd, regionTagId %in% tags$regionTagId)
overlays<-subset(overlays,overlayId %in% unique(maps$overlayId))
aois<-subset(aois,overlayId %in% unique(maps$overlayId))


save(aois,overlays,maps,views,taggers,tags,crowd,file=paste0(pathToGit,"compiledData.RData"))

