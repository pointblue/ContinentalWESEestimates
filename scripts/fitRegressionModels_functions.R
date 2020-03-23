# TODO: Add comment
# 
# Author: lsalas
###############################################################################


## The purpose of this file is to store functions to be used during the linear modeling of WESE distributions in the continent

## This is a utility function to load a shapefile to the environment as a geospatial data.frame (e.g., a spatial polygons data.frame)
## pathToGit is the path to the local Git folder
## folderName is the name of the folder containing the shape file
## shapeName is the name of the shape file (.shp)
readShapeFile<-function(pathToGit,folderName="glaciers",shapeName="glaciers"){
	library("rgdal");library("sp");library("rgeos")
	dsnpth<-paste0(pathToGit,"data/",folderName)
	if(file.exists(paste0(pathToGit,"data/",folderName,"/",shapeName,".shp"))){
		shpobj<-readOGR(dsnpth,shapeName)
	}else{
		shpobj<-paste("ERROR: the shape file was not found at",paste0(pathToGit,"data/",folderName,"/",shapeName,".shp. Check spelling of pathToGit, folderName and shapeName? Is it a shape file?"))
	}
	return(shpobj)
}

## This function overlaps the existing data.frame with a spatial polygons object to get attribution from the polygons
## shpobj is the spatial data object from which we will attribute
## attribName is the name of a field in shpobj that will be used to attribute dat
## dat is the data.frame we want to attribute
## datKey is a unique key for the table dat, usually (and expected to be) gridCellId
## lonfld and latfld are the names of lon and lat fields (or their equivalents if projected) in dat
## datproj is the projection of the dat file, in PROJ format
attributeWithShape<-function(shpobj,attribName,dat,datKey="gridCellId",lonfld,latfld,datproj){
	library("sp");library("raster")
	
	shpNames<-names(shpobj)
	if(!TRUE %in% c(attribName %in% shpNames)){
		datdf<-paste("ERROR: the variable",attribName,"was not found in the shape object. Check spelling?")
	}else{
		datNames<-names(dat)
		if(!TRUE %in% c(lonfld %in% datNames) | !TRUE %in% c(latfld %in% datNames)){
			datdf<-paste("ERROR: either",latfld,"or",lonfld,"or both were not found in the data file. Check spelling?")
		}else{
			if(!TRUE %in% c(datKey %in% datNames)){
				datdf<-paste("ERROR: the key",datKey,"was not found in the data file. Check spelling?")
			}else{
				#make spatial points from dat
				spp<-dat[,c(datKey,lonfld,latfld)]
				coordinates(spp)<-c(londf,latfld)
				proj4string(spp)<-CRS(datproj)
				
				#reproject to shape
				shpprj<-projection(shpobj)
				spp<-spTransform(spp,CRS(shpprj))
				
				#get over from shape to attribute with values of attribName
				sppat<-over(spp,shpobj)
				spp$attribVal<-sppat[,attribName]
				sppdf<-as.data.frame(spp)
				names(sppdf)<-gsub("attribVal",attribName,names(sppdf))
				
				#merge attributed df to original dat
				datdf<-merge(datdf,sppdf[,c(datKey,attribName)],by=datKey,all.x=TRUE)
			}
		}		
	}
	
	#return
	return(datdf)
}
