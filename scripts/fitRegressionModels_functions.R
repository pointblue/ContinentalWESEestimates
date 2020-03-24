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
attributeWithShape<-function(shpobj,attribName,data,datKey="gridCellId",lonfld,latfld,dataproj){
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

## This function creates a bootstrapped list of subsamples of the data so that each subsample has the same number of cells with and without seals
## data is the data.frame to subsample from
## nsamples is the number of bootstrap samples to take, defaults to 100
## hasMapsBehavior determines if the presence of WESE is defined by abundance > 0 (default) or the presence of maps.
## stratifyByClusters 
## countVar is the name of the count variable in the data, either mdlCol or mdlIsl
## setBinomial indicates if the outpout samples are binomial presence/absence or abundance (default)
bootSampleWESEdata<-function(data,nsamples=100,hasMapsBehavior=0,stratifyByClusters=0,countVar="mdlCol",setBinomial=FALSE){
	#limit the clustering to max 5
	kv<-stratifyByClusters;if(kv>5){kv<-5}
	
	#divide into presence and absence sets
	#if hasMapsBehavior==1, then present = hasMap=1, not mdlCol/mdlIsl > 0
	if(setBinomial==FALSE){
		names(data)<-gsub(countVar,"abundance",names(data))
		if(hasMapsBehavior==0){
			pres<-subset(data,abundance>0)
			abst<-subset(data,abundance==0)
		}else{
			pres<-subset(data,hasMap==1)
			abst<-subset(data,hasMap==0)
		}
	}else{
		names(data)<-gsub(countVar,"presence",names(data))
		data$presence<-ifelse(data$presence>0,1,0)
		if(hasMapsBehavior==0){
			pres<-subset(data,presence==1)
			abst<-subset(data,abundance==0)
		}else{
			pres<-subset(data,hasMap==1);pres$presence<-1
			abst<-subset(data,hasMap==0)
		}
		
	}
	
	psize<-nrow(pres)
	
	#if stratifyByClusters > 0, then 
	if(stratifyByClusters==0){
		abst$kval<-1
	}else{
		#create and add the cluster attributions to abst
		library(cluster)
		kmm<-kmeans(abst[,which(!names(abst) %in% c("gridCellId","nearLineId","near_x","near_y","adpecol","empecol","coords.x1","coords.x2","mdlCol","mdlIsl","hasMaps"))], 
				centers=kv) 
		abst$kval<-kmm$cluster
	}
	binvals<-unique(abst$kval)	#this is the set of unique bin names
	asizek<-round(psize/max(binvals))	#this is how much to survey from each cluster
	#For each bootstrap:
	nbdf<-abst[1:(2*psize),]
	#dimensioning...
	res<-llply(.data=1:nsamples, .fun=function(x,nbdf){
				return(nbdf)
			},nbdf=nbdf)
	for(x in 1:nsamples){
		#sample with replacement from each cluster
		if(max(binvals)==1){	#just take a random sample from absents of size psize
			rv<-sample(1:nrow(abst),size=psize,replace=TRUE)
			rdf<-abst[rv,]
		}else{
			#sample at random from each cluster, of size asizek, then combine into df
			rdf<-ldply(.data=binvals, .fun=function(x,abst,asizek){
						tdf<-subset(abst,kval==x)
						rv<-sample(1:nrow(tdf),size=asizek,replce=TRUE)
						ttdf<-tdf[rv,]
						return(ttdf)
					})
		}
		#add the presence data
		rdf<-rbind(pres,rdf[,names(pres)])
		res[[x]]<-rdf
		#take the next sample
	}
	
	#return the list
	return(res)
}

## This function fits a regression model formula to a list of bootstrapped datasets, returning the mean coefficients, etc., and metrics of GOF
## fml is a string specifying the formula to fit
## datalist is a list object with each element being a data.frame on which to fit the model
## fam is a string naming the error link, defaulting to "gaussian", which will result in a normal regression
fitModelToBootstrap<-function(fml="abundance~1",datalist,fam="gaussian"){
	#outputs: a list with (df of coefficients, df of std errors, df of t-values, df of p-values), each df has as many columns as datasets in datalist
	#outputs: a df with adjR2, df, aic, loglik
	#outputs: a df with each column being the residuals from each model fit
	
	mdllst<-llply(.data=datalist, .fun=function(dd,fml,fam){
				if(fam=="gaussian"){
					mdl<-lm(formula=fml,data=dd)
				}else{
					mdl<-glm(formula=fml,data=dd,family=fam)
				}
				return(mdl)
		},fml=fml,fam=fam)
	
	coeflst<-list()
	for(pp in 1:4){
		pardf<-ldply(.data=mdllst, .fun=function(mm){
					svc<-summary(mm)$coefficients
					tdf<-as.data.frame(t(svc[,pp]))
					return(tdf)
				})
		coeflst[[pp]]<-pardf
	}
	names(coeflst)<-c("Coefficients","St.Errors","t_values","p_values")
		
		
	gofdf<-ldply(.data=mdllst, .fun=function(mm,fam){
				dfv<-summary(mm)$df[2]
				aicv<-AIC(mm)
				loglikv<-as.numeric(logLik(mm))
				if(fam=="gaussian"){
					rsqv<-summary(mm)$r.square
					tdf<-data.frame(Df=dfv,AIC=aicv,LogLik=loglikv,Rsq=rsqv)
				}else{
					resdev<-mm$deviance
					tdf<-data.frame(Df=dfv,AIC=aicv,LogLik=loglikv,ResidDeviance=resdev)
				}
				return(tdf)
		},fam=fam)
	
	#better to dimension the residuals data.frame ahead of estimation...
	ndat<-NROW(datalist);nrec<-nrow(datalist[[1]])
	residdf<-as.data.frame(matrix(rep(9.999,(ndat*nrec)),ncol=ndat))
	names(residdf)<-paste("resMdl",1:ndat,sep="_")
	for(rr in 1:ndat){
		resids<-mdllst[[rr]]$residuals
		if(NROW(resids)<nrec){
			diffrec<-nrec-NROW(resids)
			padres<-rep(NA,times=diffrec)
			resids<-c(resids,padres)
		}
		residdf[,rr]<-resids
	}
		
	res<-list(models=mdllst,coefs=coeflst,gofs=gofdf,resids=residdf)
	
	return(res)
		
}

## This function summarizes the results from fitting a model to an ensemble of data with the fitModelToBootstrap function
## fitobj is the object resulting from the fitModelToBootstrap function
## what indicates what we want to summarize: coefs, gof, resids
summarizeResults<-function(fitobj,what="coefs"){
	resdf<-"Nothing summarized. Make sure to use the function for the results of fitModelToBootstrap, with 'what' being either one of: 'coefs' (default), 'gof', or 'resids'."
	if(what=="coefs"){
		obj<-fitobj$coefs
		#this has 4 data.frames, one each for: coefficients, standard errors, t-values and p-values
		#each data frame has as many rows as there are coefficients, and as many columns as there are bootstrap samples
		#so, we calculate the average of each data.frame across columns, and then cbind the results
		coefdf<-obj$Coefficients;avgcoef<-apply(X=coefdf,MARGIN=2,"mean")
		parnams<-names(avgcoef); nbt<-ncol(coefdf)
		stedf<-obj$St.Errors;avgste<-apply(X=stedf,MARGIN=2,"mean")
		tvaldf<-obj$t_values;avgtval<-apply(X=tvaldf,MARGIN=2,"mean")
		pvaldf<-obj$p_values;avgpval<-apply(X=pvaldf,MARGIN=2,"mean")
		resdf<-data.frame(Parameter=parnams,Coefficient=avgcoef,StError=avgste,t_value=avgtval,Prob_t=avgpval,Nboot=nbt)
	}
	if(what=="gof"){
		gofdf<-fitobj$coefs
		gofres<-apply(X=gofdf, MARGIN=1, "mean")
		resdf<-data.frame(Parameter=names(gofdf),Value=gofres)
	}
	if(what=="resids"){
		residsdf<-fitobj$resids
		resdf<-apply(X=residsdf, MARGIN=2, "mean")
	}
	return(resdf)
}

