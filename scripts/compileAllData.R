# TODO: Add comment
# 
# Author: lsalas
###############################################################################

## Regions: (loop)
## Amundsen (AMU); East Antarctica 1 (EA1) and 2 (EA2), Queen Maud (QMA); Ross (RSS) and WAP (WAP)
regions<-data.frame(region=c("Amundsen","EastAnt1","EastAnt2","QueenMaud","Ross","WAPS"),
		name=c("AMU","EA1","EA2","QMA","RSS","WAP"),stringsAsFactors=F)
pathToFiles<-"//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/ContinentalEst2/"

laois<-list();loverlays<-list();lmaps<-list();lviews<-list();
ltaggers<-list();ltags<-list();lcrowd<-list()
for(rr in 1:6){
	rg<-regions[rr,"region"];nm<-regions[rr,"name"]
	pth<-paste0(pathToLocalGit,rg,"_additional_files/",rg,".RData")
	load(pth)
	aois$region<-nm;overlays$region<-nm;maps$region<-nm;views$region<-nm;taggers$region<-nm;tags$region<-nm;crowd$region<-nm
	laois[[rr]]<-aois;loverlays[[rr]]<-overlays;lmaps[[rr]]<-maps;lviews[[rr]]<-views
	ltaggers[[rr]]<-taggers;ltags[[rr]]<-tags;lcrowd[[rr]]<-crowd
}
aois<-do.call("rbind",laois);overlays<-do.call("rbind",loverlays);maps<-do.call("rbind",lmaps);views<-do.call("rbind",lviews)
taggers<-do.call("rbind",ltaggers);tags<-do.call("rbind",ltags);crowd<-do.call("rbind",lcrowd)

save(aois,overlays,maps,views,taggers,tags,crowd,file=paste0(pathToFiles,"compiledData.R"))

