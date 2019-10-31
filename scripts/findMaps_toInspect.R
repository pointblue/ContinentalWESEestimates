# TODO: Add comment
# 
# Author: lsalas
###############################################################################


## The goal of this code file is to find maps that ML should inspect to increase the number of taggers that have evaluated against her
## and to ensure we have information for taggers from all regions.
## To keep the repo clean, I am not including the first set of data files used with this code

library(plyr)

## Load the data
pth<-"//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/Data/"
regions<-data.frame(fnam=c("east_ant1_additional_files/","east_ant2_additional_files/","amundsen_additional_files/",
		"queen_maud_additional_files/","ross_additional_files/","waps_additional_files/","weddell_crowdrank/"),
		dnam=c("eastAnt1","eastAnt2","Amundsen","QueenMaud","Ross","WAPS","Weddell")
)

## for each region...
## identify the maps inspected by ML
## For taggers overlapping with ML, determine how many maps and how many tags ML placed in these maps
## Consider taggers not overlapping with at least 20 seals with ML as incomplete
## (need tags and views merged)

## For maps not inspected by ML
## List the number of non-complete taggers in each
## Sort descending
needMaps<-llply(.data=1:(nrow(regions)),.fun=function(rr,regions,pth){
			reg<-regions[rr,"fnam"];datn<-regions[rr,"dnam"]
			load(paste0(pth,reg,datn,".RData"))
			tv<-merge(tags,views,by=c("mapViewId","taggerId"),all.x=T)
			chk<-sum(is.na(tv$mapId))
			if(chk>0){
				print(paste("Region",datn,"has unlinked tags:",chk))
				tv<-subset(tv,!is.na(mapId))
			}
			tv$numTags<-1
			
			#find ML maps and the number of tags she placed in each
			mlTags<-subset(tv,taggerId==21758509)
			mlMaps<-unique(mlTags$mapId)
			mlm<-aggregate(numTags~mapId,mlTags,sum)
			names(mlm)<-gsub("numTags","numMLtags", names(mlm), fixed=T)
			
			#Here's the tag data for maps inspected by ML
			tviml<-subset(tv,mapId %in% mlMaps & taggerId!=21758509)
			
			#Get the list of taggers and see how many mlTags they overlapped with ML
			uniqTaggers<-unique(tags$taggerId);uniqTaggers<-subset(uniqTaggers,uniqTaggers!=21758509)
			overTaggers<-ldply(.data=uniqTaggers,.fun=function(x,tviml,mlm){
						ttger<-subset(tviml,taggerId==x)
						if(nrow(ttger)==0){
							ntt<-0;ntml<-0
						}else{
							stger<-aggregate(numTags~mapId,ttger,sum)
							stml<-merge(stger,mlm,by="mapId",all.x=T)
							ntt<-sum(stml$numTags);ntml<-sum(stml$numMLtags)
						}
						tdf<-data.frame(taggerId=x,numTags=ntt,numMLtags=ntml)
						return(tdf)
					},tviml=tviml,mlm=mlm)
			
			#Those for which we need to find good overlapping maps are the ones with fewer than 20 numMLtags
			incompTaggers<-c(unique(subset(overTaggers,numMLtags<20)$taggerId),unique(subset(tv,!mapId %in% mlMaps)$taggerId))
			compTaggers<-unique(subset(overTaggers,numMLtags>19)$taggerId)
			
			#Find the maps where the compatible taggers overlapped with ML
			compMaps<-sapply(X=compTaggers,FUN=function(q,tviml,mlMaps){
						cmpm<-unique(subset(tviml,taggerId==q & (mapId %in% mlMaps))$mapId)
						return(cmpm)
					},tviml=tviml,mlMaps=mlMaps)
			compMaps<-unlist(compMaps)
			compMaps<-unique(compMaps)
			
			#list the incompatible maps and the number of taggers that visited each
			tvnml<-subset(tv,(!mapId %in% compMaps) & (!mapId %in% mlMaps))
			incompMaps<-unique(tvnml$mapId)
			inspectMaps<-ldply(.data=incompMaps,.fun=function(mm,tvnml){
						timp<-subset(tvnml,mapId==mm)
						nit<-NROW(unique(timp$taggerId))
						aggnit<-aggregate(numTags~taggerId,timp,sum)
						avgtags<-mean(aggnit$numTags);maxtags<-max(aggnit$numTags)
						tdf<-data.frame(mapId=mm,numTaggers=nit,avgtags=avgtags,maxtags=maxtags)
						return(tdf)
					},tvnml=tvnml)
			
			inspectMaps<-inspectMaps[order(inspectMaps$numTaggers,decreasing=TRUE),]
			inspectMaps$region<-datn
			reslst<-list(overTaggers=overTaggers,incompTaggers=incompTaggers,compTaggers=compTaggers,
					compMaps=compMaps,incompMaps=incompMaps,inspectMaps=inspectMaps)
			
			return(reslst)
			
		},regions=regions,pth=pth)

names(needMaps)<-regions$dnam
		

## Inspect results, send to ML
for(nn in regions$dnam){
	tt<-needMaps[[nn]]$inspectMaps
	print(paste("Maps in region",nn,"have a maximum number of overlapping taggers per map:",max(tt$numTaggers)))
}

## bind the results and filter for a reasonable number of overlaps
inspMaps<-ldply(.data=regions$dnam,.fun=function(nn,needMaps){
			tt<-needMaps[[nn]]$inspectMaps
			return(tt)
		},needMaps=needMaps)

setInspMaps<-subset(inspMaps,numTaggers>3 & avgtags>5)

