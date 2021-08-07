# TODO:
# 
# Author: lsalas
###############################################################################


## The following are utility functions to be used with the file getNICice.R


## FUNCTION to download and unzip the target NIC data
# filename is the name of the NIC file to retrieve - obtained with function getNICfilename
# returns the name of the unpacked file, including its path
nicDownload<-function(filename,nicsavedir){
	Filelocation<-paste0("https://www.natice.noaa.gov/pub/weekly/antarctic/",filename[1])
	savename<-paste0(nicsavedir,"/",substr(filename[1],29,46))	
	download.file(Filelocation,destfile=savename,method="libcurl")
	unzip(zipfile=savename,exdir=paste0(nicsavedir,"/",substr(filename[1],29,40)))
	return(savename)
}

## FUNCTION to access the downloaded NIC data
# fileloc is the full path to the downloaded zipped NIC data file 
# dataproj is the projection of the target data, here named "primaryproj"
# nicsavedir is the directory where the tool saves NIC data downloads
getFastIce<-function(fileloc,dataproj,nicsavedir){
	print("processing single fast ice date")
	
	#Read shapefile
	loopdsn <-substr(fileloc,1,nchar(fileloc)-4)
	dirnam<-str_sub(fileloc,-16,-5)
	shplst<-list.files(paste0(nicsavedir,dirnam),pattern=".shp")
	shpname<-substr(shplst[1],1,regexpr(".shp",shplst[1])-1)
	
	#read shapefile
	region <- readOGR(dsn=loopdsn,shpname)
	
	#Project to match full points shapefile
	region<-spTransform(region,dataproj)
	
	# Define fast ice and get region layer for creating land and ocean edges
	region@data$FP <- ifelse(region@data$FP == "08",8,0)
	#get fast ice only regions
	fast<-subset(region,region@data$FP == 8)
	
	functionList <- list("fast" = fast, "region" = region)
	
	return(functionList)
	
}

## FUNCTION to retrieve the NIC data filenames, one of which would be used to generate covariates
# getmonth is the name of a month, first three characters, first inupper case; defaults to "Nov"
# getyear is the 4-number year; defaults to 2011
getNICfilename<-function(getmonth="Nov",getyear=2011){
	# To get NIC data... Remove warns of certificate at NIC
	set_config( config( ssl_verifypeer = 0L ) )
	
	#we send request here
	urlv<-"https://www.natice.noaa.gov/products/weekly_products.html"
	
	#make sure to ignore security cert error
	set_config( config( ssl_verifypeer = 0L ) )
	
	dd<-ifelse(getmonth=="Feb","28",
			ifelse(getmonth %in% c("Apr","Jun","Sep","Nov"),"30","31"))
	#this is what we want
	body <- list(			
			area="Antarctic",
			day0="01",
			day1=dd,
			format="Shapefiles",
			month0=getmonth,
			month1=getmonth,
			oldarea="Antarctic",
			oldformat="Shapefiles",
			subareas="Hemispheric",
			year0=as.character(getyear),
			year1=as.character(getyear))
	
#request the data
	r <- POST(urlv, body = body, encode = "form")
	
#check that it worked
	r$status_code==200
	
#it returns XML, and that can be a pain to parse.
#yet, we know we want something that looks like this: ...shapefiles/hemispheric/antarc999999
#so, using regular expressions to search for it...
	nv<-as.numeric(gregexpr("shapefiles/hemispheric/antarc",content(r,encoding="UTF-8"))[[1]])
	
	filenames<-character()
	#here getting the address including year
	for(nn in nv){
		filenames<-c(filenames,substr(content(r,encoding="UTF-8"),nn-5,nn+38))
	}
	
#the above include the .html and .xml files. Filter for zips only...
	filenames<-subset(filenames,grepl("^[0-9][0-9][0-9][0-9].*hemispheric.*.zip",filenames))  
}

## FUNCTION to get edges from fast ice areas
# areas is a list of shapefiles (spatialPolygons or spatialPolygonsDataFrame) which is the output of the getFastIce function 
# areas has two named elements: region and fast. See getFastIce function for details.
getEdges<-function(areas){
	region<-myareas$region
	#fast<-myareas$fast
	
	# Dissolve regions based on fast ice or not
	unionfp <- gUnaryUnion(region, id = region@data$FP)  
	
	#create line of fast ice edge away from continent
	oceanedge = gDifference(
			as(unionfp,"SpatialLines"),
			as(gUnaryUnion(unionfp),"SpatialLines"),
			byid=TRUE)
	
	#create line of fast ice edge nearest to continent
	landedge<-as(gUnaryUnion(unionfp),"SpatialLines")
	
	return(list=c(oceanedge=oceanedge,landedge=landedge))
}

## FUNCTION to plot grid points, and their nearest land points
# subdf is the data.frame obtained from a spatial_points data.frame with all points on fast ice for the NIC date chosen
# edges is a spatial_lines data.frame outlining the continent's edge
# myareas is a list with the spatial polygons of fast ice for the NIC date chosen
showPoints<-function(subdf,edges,myareas){
	
	edgesldf<-SpatialLinesDataFrame(edges$landedge, data=data.frame(ID=1))
	edgedf<-fortify(edgesldf)
	pc<-ggplot(data=edgedf, aes(x=long, y=lat)) + geom_path(aes(group=group)) + theme_bw()
	pc<-pc + geom_point(data=subdf,aes(x=coords.x1,y=coords.x2),color="blue",size=0.2) + geom_point(data=subdf,aes(x=near_x,y=near_y),color="red",size=0.1)
	
	icedf<-fortify(myareas$fast)
	subedge<-subset(edgedf,long>100000 & long<1000000 & lat< -500000 & lat > -1500000)
	subsubdf<-subset(subdf,coords.x1>100000 & coords.x1<1000000 & coords.x2< -500000 & coords.x2 > -1500000)
	subicedf<-subset(icedf,long>100000 & long<1000000 & lat< -500000 & lat > -1500000)
	pe<-ggplot(data=subedge, aes(x=long, y=lat)) + geom_path(aes(group=group)) + theme_bw()
	pe<-pe + geom_point(data=subsubdf,aes(x=coords.x1,y=coords.x2),color="blue",size=0.2) + geom_point(data=subsubdf,aes(x=near_x,y=near_y),color="red",size=0.1)
	pe<- pe + geom_path(data=subicedf,aes(x=long, y=lat,group=group),color="gray")
	
	print(pc)
	dev.new();print(pe)
	
}

## FUNCTION to make segments in 8 directions from each grid point
# pt is the point to use
# dst is the length in m for the segment
makeSegments<-function(pnt,dst,rosprj){
	px<-as.integer(pnt$coords.x1);py<-as.integer(pnt$coords.x2)
	seqx<-round(c(cos(0),cos(pi/4),0,cos(3*pi/4),cos(pi),cos(3*pi/4),0,cos(pi/4)) *dst)
	seqy<-round(c(0,sin(pi/4),sin(pi/2),sin(3*pi/4),0,sin(5*pi/4),sin(3*pi/2),sin(7*pi/4)) *dst)
	linlst<-list()
	for(ss in c(1:8)){
		ll<-cbind(c(px,px+seqx[ss]), c(py,py+seqy[ss]))
		Sl <- Line(ll)
		Sll <- Lines(list(Sl), ID = letters[ss])
		linlst[[ss]]<-Sll
	}
	ptlines<-SpatialLines(linlst)
	projection(ptlines)<-CRS(rosprj)
	return(ptlines)
}

##FUNCTION to intersect segments with ice polygons and determine shortest distance to edge, fast ice width, total ice distance
# ptcoord are the xy coordinates of the grid point
# gplines is the spatial lines data frame with the 8 line segments
# fip is the fast ice polygons data frame
# neland are the xy coordinates of the nearland point
# edgl is the spatial lines for the land's edge
getNearestIceEdge<-function(ptcoord,gplines,fip,neland,edgl){
	xld<-as.numeric(neland[1]);yld<-as.numeric(neland[2])
	xgp<-as.numeric(ptcoord[1]);ygp<-as.numeric(ptcoord[2])
	distdf<-data.frame()
	for(dd in 1:8){
		tst<-gIntersection(gplines[dd],edgl)
		if(!is.null(tst)){
			nlp<-nrow(coordinates(tst))
			landint<-TRUE
		}else{
			nlp=0;landint<-FALSE
		}
		aa<-gIntersection(gplines[dd],fip)
		outpt<-getIntersectPoints(aa,xgp,ygp)
		xout<-as.numeric(outpt[1]);yout<-as.numeric(outpt[2])
		dld<-round(sqrt(((xld-xout)^2)+((yld-yout)^2)))
		dgp<-round(sqrt(((xgp-xout)^2)+((ygp-yout)^2)))
		tdf<-data.frame(dir=dd,xout=xout,yout=yout,distland=dld,distedge=dgp,nlp=nlp,landint=landint)
		distdf<-rbind(distdf,tdf)
	}
	distdf$diffdist<-distdf$distland-distdf$distedge
	
	#retrieve the edgepoint with smallest positive diffdist that either has no landint, or even-numbered landint
	totID<-sum(distdf$distedge)
	if(FALSE %in% landint){
		distdf<-subset(distdf,landint==FALSE & diffdist>0)
	}else{
		distdf<-subset(distdf, (nlp %% 2 == 0) & diffdist>0)
	}
	if(nrow(distdf)>0){
		distdf<-distdf[order(distdf$distedge,decreasing=FALSE),]
	}else{
		distdf$xout<-NA;distdf$yout<-NA;distdf$distedge<-NA
	}
	
	resdf<-data.frame(edge_x=distdf[1,"xout"],edge_y=distdf[1,"yout"],distEdge=distdf[1,"distedge"],totalIceDist=totID)
	
	return(resdf)
}

## FUNCTION to extract the coordinates from the intersection points of the lines and the fast ice polygon
# aa is the intersection object
# xgp and ygp are the coordinates of the grid point
getIntersectPoints<-function(aa,xgp,ygp){
	intdf<-data.frame()
	for(ii in 1:length(aa)){
		intc<-coordinates(aa)[[ii]][[1]][2,]
		intct<-data.frame(intx=intc[1],inty=intc[2])
		intdf<-rbind(intdf,intct)
	}
	if(nrow(intdf)>1){
		intdf$disttmp<-sqrt(((intdf$intx-xgp)^2)+((intdf$inty-ygp)^2))
		intdf<-intdf[order(intdf$disttmp),]
	}
	return(intdf[1,c("intx","inty")])
}

## FUNCTION to get fast ice width for each grid point
# savename is the name of the downloaded NIC ice data, including its path
# primaryproj is the projection used by the spatial points file of grid points
# myareas is the selected fast ice data (the areas with fast ice for the date chosen) - a list of spatial polygons 
# studyarea_pointswLand is the spatial points table of grid points, already attributed with nearest land point
# buffwidth is the size (in meters) of a circular buffer to count how many 50-m edgepoints are within it - a metric of fast ice stability and abundance
# plotit indicates if plots of the continent-wide and Erebus-only grid points and their nearest land points be made
calcFasIceWidth<-function(savename,primaryproj,myareas,studyarea_pointswLand,buffwidth=20000,plotit=TRUE){
	
	#get the lines for land (landedge) and fast ice (ocenedge) edges: 
	edges<-getEdges(areas=myareas)
	
	#Subset ALL points to only those within fast ice
	fast<-myareas$fast
	subsetpoints <- studyarea_pointswLand[fast,]
	subdf<-as.data.frame(subsetpoints)
	
	if(plotit==TRUE){
		#see what we got:				************************************************ Optional
		showPoints(subdf=as.data.frame(subsetpoints),edges=edges,myareas=myareas)
	}
	
	#get fast ice edge and add points 50m along it
	oedge<-edges$oceanedge
	npts<-round(lineLength(oedge)/50)	
	icedgedf<-spsample(oedge,n=npts,type="regular")
	icedgedf<-data.table(as.data.frame(icedgedf))
	
	#now loop through the subsetpoints to find the nearest point in icedgedf
	subdf$iceedge.x1<-NA;subdf$iceedge.x2<-NA;subdf$distNearestIceEdge<-NA;subdf$fastIceAbund<-NA
	edgpts<-icedgedf[,c("x","y")]
	coordinates(edgpts)<-c("x","y")
	projection(edgpts)<-CRS(projection(subsetpoints))
	
	#time it...
	tm<-Sys.time()
	for(rr in 1:(nrow(subdf))){
		icedgedf$gptx<-subdf[rr,"coords.x1"];icedgedf$gpty<-subdf[rr,"coords.x2"]
		icedgedf[,dist:=sqrt(((x-gptx)^2)+((y-gpty)^2))]
		icedgedf<-setorder(icedgedf,dist)
		
		tdf<-data.frame(gptx=icedgedf[1,"gptx"],gpty=icedgedf[1,"gpty"])
		coordinates(tdf)<-c("gptx","gpty")
		projection(tdf)<-CRS(projection(subsetpoints))
		pntbuff<-gBuffer(tdf,width=20000)
		qq<-over(tdf,pntbuff)
		subdf[rr,"iceedge.x1"]<-icedgedf[1,"x"];subdf[rr,"iceedge.x2"]<-icedgedf[1,"y"];
		subdf[rr,"distNearestIceEdge"]<-icedgedf[1,"dist"];subdf$fastIceAbund<-sum(!is.na(qq))
	}
	proctim<-Sys.time()-tm
	print(paste("Processing time:",proctim))
	
	#ice width is the sum of dist to edge and distance to land
	subdf$fastIceWidth<-subdf$distToShore+subdf$distNearestIceEdge
	
	### Convert back to spatial points...
	FastIcePoints<-subdf
	coordinates(FastIcePoints)<-c("coords.x1","coords.x2")
	proj4string(FastIcePoints)<-primaryproj
	
	return(FastIcePoints)
}

## Function to get the data for fast ice stability, persistence, and predictability
# spdf is the spatial.points data.frame
# iceyear the ice year being sought from NIC
# iceareas the set of fast ice polygons for the year and date requested
getFastIceSPPdata<-function(spdf,iceyear,iceareas,keyfield="pointid",nicsavedir){
	
	primaryproj<-projection(spdf)
	
	if(nrow(iceareas)>0){#Subset ALL points to only those within desired areas
		targetpoints <- spdf[iceareas,]
	}else{
		targetpoints <- spdf
	}
	
	
	## We want some specific metrics of fast ice:
	## Three values: 
	
	###########
	#	stability (is there fast ice at the end of December?)
	basedf<-data.frame(recId=as.character(as.data.frame(targetpoints)[,keyfield]))
	stabilitydf<-basedf
	
	filename<-getNICfilename(getmonth="Dec",getyear=iceyear)	
	lastDec<-NROW(filename);stabilityDec<-nicDownload(filename=filename[lastDec],nicsavedir=nicsavedir)
	lastM_areas<-getFastIce(fileloc=stabilityDec,dataproj=primaryproj,nicsavedir=nicsavedir)
	lastMPoints <- spdf[lastM_areas$fast,]
	lastMdf<-data.frame(recId=as.character(as.data.frame(lastMPoints)[,keyfield]),DecemberIcePresence=1)
	stabilitydf<-merge(stabilitydf,lastMdf,by="recId",all.x=TRUE)
	stabilitydf$DecemberIcePresence<-ifelse(is.na(stabilitydf$DecemberIcePresence),0,1)
	
	###########
	#	persistence (is there fast ice throughout February? Derived metric: is persistent if it’s there the past 2 years)
	persistencedf<-basedf
	thisyear<-as.numeric(iceyear);year1prior<-as.numeric(iceyear)-1;year2prior<-as.numeric(iceyear)-2
	
	filename<-getNICfilename(getmonth="Feb",getyear=thisyear)	
	lastFeb<-NROW(filename);persistFeb0y<-nicDownload(filename=filename[lastFeb],nicsavedir=nicsavedir)	
	lastM_areas<-getFastIce(fileloc=persistFeb0y,dataproj=primaryproj,nicsavedir=nicsavedir)
	lastMPoints <- spdf[lastM_areas$fast,]
	lastMdf<-data.frame(recId=as.character(as.data.frame(lastMPoints)[,keyfield]),February0Prior=1)
	persistencedf<-lastMdf
	persistencedf$February0Prior<-ifelse(is.na(persistencedf$February0Prior),0,1)
	
	filename<-getNICfilename(getmonth="Feb",getyear=year1prior)	
	lastFeb<-NROW(filename);persistFeb1y<-nicDownload(filename=filename[lastFeb],nicsavedir=nicsavedir)	
	lastM_areas<-getFastIce(fileloc=persistFeb1y,dataproj=primaryproj,nicsavedir=nicsavedir)
	lastMPoints <- spdf[lastM_areas$fast,]
	lastMdf<-data.frame(recId=as.character(as.data.frame(lastMPoints)[,keyfield]),February1Prior=1)
	persistencedf<-merge(persistencedf,lastMdf,by="recId",all.x=TRUE)
	persistencedf$February1Prior<-ifelse(is.na(persistencedf$February1Prior),0,1)
	
	filename<-getNICfilename(getmonth="Feb",getyear=year2prior)	
	last2Feb<-NROW(filename);persistFeb2y<-nicDownload(filename=filename[last2Feb],nicsavedir=nicsavedir)	
	lastM_areas<-getFastIce(fileloc=persistFeb2y,dataproj=primaryproj,nicsavedir=nicsavedir)
	lastMPoints <- spdf[lastM_areas$fast,]
	lastMdf<-data.frame(recId=as.character(as.data.frame(lastMPoints)[,keyfield]),February2Prior=1)
	persistencedf<-merge(persistencedf,lastMdf,by="recId",all.x=TRUE)
	persistencedf$February2Prior<-ifelse(is.na(persistencedf$February2Prior),0,1)
	
	persistencedf$Persistence3Years<-apply(persistencedf[,2:4],1,sum)
	persistencedf<-persistencedf[,c("recId","Persistence3Years")]
	
	###########
	# 	predictability (derived from stability: stable in the past 5 years, but use the number of stable ice years).
	#	Must be predictable by earliest NIC date in October
	predictabilitydf<-basedf
	year3prior<-as.numeric(iceyear)-3
	year4prior<-as.numeric(iceyear)-4
	year5prior<-as.numeric(iceyear)-5
	
	#Y1
	filename<-getNICfilename(getmonth="Oct",getyear=year1prior)	
	# This is what we change from here down: lastDec<-NROW(filename);stabilityDec<-nicDownload(filename=filename[lastDec],nicsavedir=nicsavedir)
	predictOct<-nicDownload(filename=filename[1],nicsavedir=nicsavedir)
	lastM_areas<-getFastIce(fileloc=predictOct,dataproj=primaryproj,nicsavedir=nicsavedir)
	lastMPoints <- spdf[lastM_areas$fast,]
	lastMdf<-data.frame(recId=as.character(as.data.frame(lastMPoints)[,keyfield]),Oct1yP=1)
	predictabilitydf<-merge(predictabilitydf,lastMdf,by="recId",all.x=TRUE)
	predictabilitydf$Oct1yP<-ifelse(is.na(predictabilitydf$Oct1yP),0,1)
	
	#Y2
	filename<-getNICfilename(getmonth="Oct",getyear=year2prior)	
	predictOct<-nicDownload(filename=filename[1],nicsavedir=nicsavedir)
	lastM_areas<-getFastIce(fileloc=predictOct,dataproj=primaryproj,nicsavedir=nicsavedir)
	lastMPoints <- spdf[lastM_areas$fast,]
	lastMdf<-data.frame(recId=as.character(as.data.frame(lastMPoints)[,keyfield]),Oct2yP=1)
	predictabilitydf<-merge(predictabilitydf,lastMdf,by="recId",all.x=TRUE)
	predictabilitydf$Oct2yP<-ifelse(is.na(predictabilitydf$Oct2yP),0,1)
	
	#Y3
	filename<-getNICfilename(getmonth="Oct",getyear=year3prior)	
	predictOct<-nicDownload(filename=filename[1],nicsavedir=nicsavedir)
	lastM_areas<-getFastIce(fileloc=predictOct,dataproj=primaryproj,nicsavedir=nicsavedir)
	lastMPoints <- spdf[lastM_areas$fast,]
	lastMdf<-data.frame(recId=as.character(as.data.frame(lastMPoints)[,keyfield]),Oct3yP=1)
	predictabilitydf<-merge(predictabilitydf,lastMdf,by="recId",all.x=TRUE)
	predictabilitydf$Oct3yP<-ifelse(is.na(predictabilitydf$Oct3yP),0,1)
	
	#Y4
	filename<-getNICfilename(getmonth="Oct",getyear=year4prior)	
	predictOct<-nicDownload(filename=filename[1],nicsavedir=nicsavedir)
	lastM_areas<-getFastIce(fileloc=predictOct,dataproj=primaryproj,nicsavedir=nicsavedir)
	lastMPoints <- spdf[lastM_areas$fast,]
	lastMdf<-data.frame(recId=as.character(as.data.frame(lastMPoints)[,keyfield]),Oct4yP=1)
	predictabilitydf<-merge(predictabilitydf,lastMdf,by="recId",all.x=TRUE)
	predictabilitydf$Oct4yP<-ifelse(is.na(predictabilitydf$Oct4yP),0,1)
	
	#Y5
	filename<-getNICfilename(getmonth="Oct",getyear=year5prior)	
	predictOct<-nicDownload(filename=filename[1],nicsavedir=nicsavedir)
	lastM_areas<-getFastIce(fileloc=predictOct,dataproj=primaryproj,nicsavedir=nicsavedir)
	lastMPoints <- spdf[lastM_areas$fast,]
	lastMdf<-data.frame(recId=as.character(as.data.frame(lastMPoints)[,keyfield]),Oct5yP=1)
	predictabilitydf<-merge(predictabilitydf,lastMdf,by="recId",all.x=TRUE)
	predictabilitydf$Oct5yP<-ifelse(is.na(predictabilitydf$Oct5yP),0,1)
	
	predictabilitydf$PredictabilityOct5Years<-apply(predictabilitydf[,2:6],1,sum)
	
	
	## Prepare and return output df
	outdf<-merge(stabilitydf,persistencedf,by="recId")
	outdf<-merge(outdf,predictabilitydf,by="recId")
	
	return(outdf)
}

## FUNCTION to re-attribute with distance to nearest ADPE and EMPE colony and with size of nearby ADPE and EMPE clusters 
# studyarea_pointswLand is the spatial points table of grid points, already attributed with nearest land point
# adpedf is the tada.frame from MLR with the appropriate colonies and sizes to use for 2011 for ADPE
# empedf - see adpedf, mutatis mutandi
getDistanceToPenguins2011<-function(studyarea_pointswLand,adpedf,empedf){
	#subset studyarea_points to only those points with fast ice
	fip<-as.data.frame(subset(studyarea_pointswLand,fastIcePresent==TRUE))
	
	#re-project penguin data to polar stereographic.
	coordinates(adpedf)<-c("Longitude","Latitude")
	proj4string(adpedf)<-CRS("+proj=longlat +datum=WGS84")
	adpe<-spTransform(adpedf,CRS(proj4string(studyarea_pointswLand)))
	adpe<-as.data.frame(adpe)
	
	coordinates(empedf)<-c("Longitude","Latitude")
	proj4string(empedf)<-CRS("+proj=longlat +datum=WGS84")
	empe<-spTransform(empedf,CRS(proj4string(studyarea_pointswLand)))
	empe<-as.data.frame(empe)
	#then calc distance and use the min
	
	#50 km cutoff
	adpedist<-ldply(.data=fip$pointid,.fun=function(ii,fip,df){
				flon<-as.numeric(subset(fip,pointid==ii)$coords.x1)
				flat<-as.numeric(subset(fip,pointid==ii)$coords.x2)
				df$ADPEdist<-sqrt(((df$Longitude-flon)^2)+((df$Latitude-flat)^2))
				seldf<-df[which(df$ADPEdist==min(df$ADPEdist)),c("ADPEname","ADPEdist")]
				seldf$pointid<-ii
				#calculate the abundance: within 100 km, sum of weighted abundances,by log of distance (in km)
				df<-subset(df,ADPEdist<=50000)
				if(nrow(df)>0){
					df$weight<-distanceWeightDecay(coff=50,df$ADPEdist)
					df$weightedAbund<-df$ADPEcount*df$weight
					seldf$ADPEabund<-sum(df$weightedAbund)
				}else{
					seldf$ADPEabund<-0
				}
				
				return(seldf)
			},fip=fip, df=adpe)
	
	#EMPE cutoff is 150 km
	empedist<-ldply(.data=fip$pointid,.fun=function(ii,fip,df, s){
				flon<-as.numeric(subset(fip,pointid==ii)$coords.x1)
				flat<-as.numeric(subset(fip,pointid==ii)$coords.x2)
				df$EMPEdist<-sqrt(((df$Longitude-flon)^2)+((df$Latitude-flat)^2))
				seldf<-df[which(df$EMPEdist==min(df$EMPEdist)),c("EMPEname","EMPEdist")]
				seldf$pointid<-ii
				df<-subset(df,EMPEdist<=150000)
				if(nrow(df)>0){
					df$weight<-distanceWeightDecay(coff=150,df$EMPEdist)
					df$weightedAbund<-df$EMPEcount*df$weight
					seldf$EMPEabund<-sum(df$weightedAbund)
				}else{
					seldf$EMPEabund<-0
				}
				return(seldf)
			},fip=fip, df=empe, s=s)
	
	penguindf<-merge(adpedist,empedist,by="pointid")
	
	return(penguindf)
}

## This function calculates the decay slope and intercept and then estimates the weight based on distance
# coff is the decay cutoff value, 50 km for ADEP and 150 km for EMPE
# cdist is the colony distance in meters
distanceWeightDecay<-function(coff,cdist){
	s<--0.9/(log(coff)-log(2.5))
	b<-1-(s*log(2.5))
	cdkm<-cdist/1000
	weight<-ifelse(cdkm<2.5,1,(log(cdkm)*s)+b)
	return(weight)
}


