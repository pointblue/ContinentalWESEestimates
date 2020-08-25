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
## data is the data.frame we want to attribute
## datKey is a unique key for the table dat, usually (and expected to be) gridCellId
## lonfld and latfld are the names of lon and lat fields (or their equivalents if projected) in dat
## datproj is the projection of the dat file, in PROJ format
attributeWithShape<-function(shpobj,attribName,data,datKey="gridCellId",lonfld,latfld,dataproj){
	
	shpNames<-names(shpobj)
	if(!TRUE %in% c(attribName %in% shpNames)){
		datdf<-paste("ERROR: the variable",attribName,"was not found in the shape object. Check spelling?")
	}else{
		datNames<-names(data)
		if(!TRUE %in% c(lonfld %in% datNames) | !TRUE %in% c(latfld %in% datNames)){
			datdf<-paste("ERROR: either",latfld,"or",lonfld,"or both were not found in the data file. Check spelling?")
		}else{
			if(!TRUE %in% c(datKey %in% datNames)){
				datdf<-paste("ERROR: the key",datKey,"was not found in the data file. Check spelling?")
			}else{
				#make spatial points from dat
				spp<-data[,c(datKey,lonfld,latfld)]
				coordinates(spp)<-c(lonfld,latfld)
				proj4string(spp)<-CRS(dataproj)
				
				#reproject to shape
				shpprj<-proj4string(shpobj)
				spp<-spTransform(spp,CRS(shpprj))
				
				#get over from shape to attribute with values of attribName
				sppat<-over(spp,shpobj)
				spp$attribVal<-sppat[,attribName]
				sppdf<-as.data.frame(spp)
				names(sppdf)<-gsub("attribVal",attribName,names(sppdf))
				
				#merge attributed df to original dat
				datdf<-merge(data,sppdf[,c(datKey,attribName)],by=datKey,all.x=TRUE)
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
bootSampleWESEdata<-function(data,nsamples=100,hasMapsBehavior=0,stratifyByClusters=0,clustByVars=NA,countVar="mdlCol",setBinomial=FALSE){
	#must first subset the data to cells <= the largest distance within each region
	#divide into presence and absence sets
	#if hasMapsBehavior==1, then present = hasMap=1, not mdlCol/mdlIsl > 0
	if(setBinomial==FALSE){
		if(hasMapsBehavior==0){ #this is the default and correct behavior
			pres<-subset(data,mdlCol>0)
			abst<-subset(data,mdlCol==0)
		}else{
			pres<-subset(data,hasMaps==1)
			abst<-subset(data,hasMaps==0)
		}
	}else{
		names(data)<-gsub(countVar,"presence",names(data))
		data$presence<-ifelse(data$presence>0,1,0)
		if(hasMapsBehavior==0){
			pres<-subset(data,presence==1)
			abst<-subset(data,presence==0)
		}else{
			pres<-subset(data,hasMap==1);pres$presence<-1
			abst<-subset(data,hasMap==0)
		}
	}
	negdf<-ldply(unique(data$Region),function(rr,abst,pres){
				regneg<-subset(abst,Region==rr)
				regpos<-subset(pres,Region==rr)
				maxdistShore<-max(regpos$distToShore)
				regneg<-subset(regneg,distToShore<=maxdistShore)
				return(regneg)
			},abst=abst,pres=pres)
	abst<-negdf
	cuts<-cut(abst$distToShore, c(-Inf,seq(500, 21500, 500), Inf),labels=c(1:44))
	abst$distbin<-as.integer(as.character(cuts))
	
	#Need to assign probs from normal of positives to the negs
	#fnrm<-fitdist(pres$logdistToShore, distr="norm")  #better fit than a gamma
	#abst$probDTS<-pnorm(abst$logdistToShore,mean=fnrm$estimate[1],sd=fnrm$estimate[2])
	
	#limit the clustering to max 5
	kv<-stratifyByClusters;if(kv>5){kv<-5}
	
	psize<-nrow(pres)
	
	if(stratifyByClusters==0 | is.na(clustByVars)){
		abst$kval<-1
	}else{
		#create and add the cluster attributions to abst
		library(cluster)
		kmm<-kmeans(abst[,which(names(abst) %in% clustByVars)], centers=kv) 
		abst$kval<-kmm$cluster
	}
	binvals<-unique(abst$kval)#this is the set of unique bin names
	asizek<-round(psize/max(binvals))#this is how much to survey from each cluster
	#For each bootstrap:
	nbdf<-abst[1:(2*psize),]
	
		
	#dimensioning...
	res<-llply(.data=1:nsamples, .fun=function(x,nbdf){return(nbdf)},nbdf=nbdf)
		
	for(x in 1:nsamples){
		rdf<-getBootSample(pres=pres,abst=abst,binvals=binvals,asizek=asizek)
		res[[x]]<-rbind(rdf,pres)
	}
	
	#return the list
	return(res)
}

## This function samples from a normal distribution to obtain a probability with which to sample from the distribution of negatives
## The goal is to obtain a distribution of distances to shore from the negative data that resembles the distribution of the positives
getBootSample<-function(pres,abst,binvals,asizek){
	#take onse sample
	#rvals<-rnorm(nrow(pres),0.5,0.2)  #the sd comes from the CV of the above normal fit
	#rvals<-ifelse(rvals<0,-1*rvals,rvals)
	#rvals<-ifelse(rvals>1,0.5,rvals)
	
	#sample with replacement from each cluster
	if(max(binvals)==1){	#just take a random sample from absents of size psize
		rdf<-ldply(1:nrow(pres),function(rr,abst){	#,rvals
					#rv<-rvals[rr]
					#tdf<-subset(abst,probDTS<rv)
					
					#get random bin to sample from
					dbv<-sample(1:44,1)	#we assigned 44 bins to the abst data
					tdf<-subset(abst,distbin==dbv)
					sdf<-tdf[sample(x=1:nrow(tdf),size=1),]
					return(sdf)
				},abst=abst)	#,rvals=rvals
		
	}else{
		#sample at random from each cluster, of size asizek, then combine into df - even if cluster=1, could just sample from this side of the function only...
		rdf<-ldply(.data=binvals, .fun=function(x,abst,asizek){	#,rvals
					ttdf<-subset(abst,kval==x)
					#krvals<-sample(rvals,size=asizek)
					
					bsamp<-ldply(1:asizek,function(rr,ttdf){	#,krvals
								rv<-rvals[rr]
								tdf<-subset(tdf,probDTS<rv)
								sdf<-tdf[sample(x=1:nrow(tdf),size=1),]
								return(sdf)
							},ttdf=ttdf)	#,krvals=krvals
					
					return(bsamp)
				},abst=abst,rvals=rvals,asizek=asizek)
		
	}
	rdf<-rdf[,which(names(rdf) %in% names(pres))]
	return(rdf)
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
	predsdf<-as.data.frame(matrix(rep(9.999,(ndat*nrec)),ncol=ndat))
	obsdf<-as.data.frame(matrix(rep(9.999,(ndat*nrec)),ncol=ndat))
	names(residdf)<-paste("resMdl",1:ndat,sep="_")
	names(predsdf)<-paste("predMdl",1:ndat,sep="_")
	names(obsdf)<-paste("obsMdl",1:ndat,sep="_")
	for(rr in 1:ndat){
		resids<-mdllst[[rr]]$residuals
		preds<-mdllst[[rr]]$fitted.values
		obs<-mdllst[[rr]]$y
		if(NROW(resids)<nrec){
			diffrec<-nrec-NROW(resids)
			padres<-rep(NA,times=diffrec)
			resids<-c(resids,padres)
			preds<-c(preds,padres)
			obs<-c(obs,padres)
		}
		residdf[,rr]<-resids
		predsdf[,rr]<-preds
		obsdf[,rr]<-obs
	}
		
	res<-list(models=mdllst,coefs=coeflst,gofs=gofdf,resids=residdf,preds=predsdf,obs=obsdf,fam=fam)
	
	return(res)
		
}

## This function summarizes the results from fitting a model to an ensemble of data with the fitModelToBootstrap function
## fitobj is the object resulting from the fitModelToBootstrap function
## what indicates what we want to summarize: coefs, gof, residplot
summarizeResults<-function(fitobj,what="coefs"){
	resdf<-"Nothing summarized. Make sure to use the function for the results of fitModelToBootstrap, with 'what' being either one of: 'coefs' (default), 'gof', or 'resids'."
	if(what=="coefs"){
		obj<-fitobj$coefs
		#this has 4 data.frames, one each for: coefficients, standard errors, t-values and p-values
		#each data frame has as many rows as there are coefficients, and as many columns as there are bootstrap samples
		#so, we calculate the average of each data.frame across columns, and then cbind the results
		coefdf<-obj$Coefficients;avgcoef<-apply(X=coefdf,MARGIN=2,"mean")
		parnams<-names(avgcoef); nbt<-nrow(coefdf)
		stedf<-obj$St.Errors;avgste<-apply(X=stedf,MARGIN=2,"mean")
		tvaldf<-obj$t_values;avgtval<-apply(X=tvaldf,MARGIN=2,"mean")
		pvaldf<-obj$p_values;avgpval<-apply(X=pvaldf,MARGIN=2,"mean")
		resdf<-data.frame(Parameter=parnams,Coefficient=round(avgcoef,3),StError=round(avgste,3),t_value=round(avgtval,3),Prob_t=round(avgpval,3),Nboot=nbt)
		if(fitobj$fam=="binomial"){
			names(resdf)<-gsub("t_value","z_value",names(resdf))
		}
	}
	if(what=="gof"){
		gofdf<-fitobj$gofs;nbt<-nrow(gofdf)
		gofres<-apply(X=gofdf, MARGIN=2, "mean")
		resdf<-data.frame(Parameter=names(gofdf),Value=gofres,Nboot=nbt)
	}
	if(what=="residplot"){
		residdf<-fitobj$resids
		predsdf<-fitobj$preds
		obsdf<-fitobj$obs
		resdf<-ldply(1:ncol(residdf),function(cc,residdf,preddf,obsdf){
					tdf<-data.frame(bootstrap=cc,predicted=predsdf[,cc],residual=residdf[,cc],observed=obsdf[,cc])
					return(tdf)
				},residdf=residdf,preddf=preddf,obsdf=obsdf)
	}
	return(resdf)
}

## This function generates a subset of the grid cell data, for only the Weddell Sea
## data is the continent-wide data.frame of WESE grid cell data
getWeddellSeaData<-function(data){
	#get coordinate values in lonlat for easier filtering
	gdf<-data[,c("gridCellId","coords.x1","coords.x2")]
	coordinates(gdf)<-c("coords.x1","coords.x2")
	proj4string(gdf)<-CRS(dataproj)
	ggdf<-spTransform(gdf,CRS("+proj=longlat +datum=WGS84"))
	ggdf<-as.data.frame(ggdf)
	names(ggdf)<-c("gridCellId","lon","lat")
	gdf<-merge(gdf,ggdf,by="gridCellId",all.x=T)
	
	#now filter for the Weddell Sea only
	wdsdf<-subset(gdf,lon > -66 & lon < -10 & lat < -64)	
	
	#there are still some cells on the western side of the peninsula, removing these
	cdf<-as.data.frame(wdsdf)[,c("gridCellId","coords.x1","coords.x2")];
	names(cdf)<-c("gridCellId","cv1","cv2")
	wdsdf<-merge(wdsdf,cdf,by="gridCellId",all.x=T)
	wdsdf<-subset(wdsdf,(cv1 > -2335000) | ((cv1 <= -2335000 & cv1 > -2410000) & cv2 > 1155000) | ((cv1 <= -2410000 & cv1 > -2455000) & cv2 > 1250000))
	
	#now filtering original data and adding the binary
	wesedf<-subset(data,gridCellId %in% wdsdf$gridCellId)
	wesedf<-wesedf[,which(!names(wesedf) %in% c("mdlIsl","hasMaps","RegionName"))]
	names(wesedf)<-gsub("mdlCol","abundance",names(wesedf))
	wesedf$presence<-ifelse(wesedf$abundance>0,1,0)
	
	return(wesedf)
	
}

## This function returns the mean, median, min and max of every numeric variable in the data, and all levels of each factor
## data is the dataset we use to bootstrap
## fml is the formula used in the model
getVarDesc<-function(data,fml){
	fml<-as.formula(fml)
	vars<-all.vars(fml[[3]])
		
	varDesc<-ldply(vars,function(nn,data){
				if(is.numeric(data[,nn])){
					minv<-min(data[,nn], na.rm = FALSE)
					maxv<-max(data[,nn], na.rm = FALSE)
					meanv<-mean(data[,nn], na.rm = FALSE)
					medv<-median(data[,nn], na.rm = FALSE)
					modv<-"NA"
				}else{
					modv<-as.character(toJSON(unique(data[,nn])))
					minv<-NA;maxv<-NA;meanv<-NA;medv<-NA;modv<-modv
				}
				tdf<-data.frame(var=nn,minv=minv,maxv=maxv,meanv=meanv,medv=medv,modv=modv)
				return(tdf)
			},data=data)
	
	return(varDesc)
}

## This function generates a prediction df from var descriptions and 1-2 selected variables to plot
## varDesc is the data frame returned by the function getVarDesc
## pdvars is the character vector naming 1 or 2 variables to vary
## useMedian is aboolean indicating if the median should be used instead of the mean for the fixed values of unvarying variables
getNewData<-function(varDesc,pdvars,useMedian=FALSE){
	if(NROW(pdvars)>2)stop("More than 2 variables to permutate?")
	vd<-subset(varDesc,!var %in% pdvars)
	pd<-subset(varDesc,var %in% pdvars)
	nt<-NROW(pdvars)
	nddf<-as.data.frame(matrix(rep(NA,times=nrow(vd)*100^nt), ncol=nrow(vd), byrow=F))
	names(nddf)<-vd$var
	for(vv in vd$var){
		if(useMedian==FALSE){cv<-subset(vd,var==vv)$meanv}else{cv<-subset(vd,var==vv)$medv}
		if(is.na(cv)){ #it's a factor - get the values and choose one at random
			cv<-subset(vd,var==vv)$modv
			modvals<-fromJSON(as.character(cv))
			cv<-sample(modvals,1)
		}
		nddf[,vv]<-rep(cv,times=nrow(nddf))
	}
	if(nt==1){
		vv<-as.character(pd$var)
		miv<-subset(pd,var==vv)$minv
		mav<-subset(pd,var==vv)$maxv
		nddf[,vv]<-seq(from=miv,to=mav,length.out=100)
	}else{
		vvs<-as.character(pd$var)
		miv<-subset(pd,var==vvs[1])$minv
		mav<-subset(pd,var==vvs[1])$maxv
		v1s<-seq(from=miv,to=mav,length.out=100)
		miv<-subset(pd,var==vvs[2])$minv
		mav<-subset(pd,var==vvs[2])$maxv
		v2s<-seq(from=miv,to=mav,length.out=100)
		vdf<-expand.grid(v1s,v2s)
		names(vdf)<-vvs
		nddf<-cbind(nddf,vdf)
	}
	
	return(nddf)
}

## This function generates a partial dependence dataset from the model bootstrap
## reslst is the list of models
## newdata is the data.frame to predict to
## pdvars is the name of the partial dependence var(s) to plot
## type indicates the type of model: lm or glm
getpdData<-function(reslst,newdata,pdvars,type="glm"){
	if(type=="lm"){
		resdf<-ldply(reslst$models,function(mm,newdata,pdvars){
					mdl<-mm
					preds<-as.data.frame(predict(mdl,newdata=newdata,interval="confidence"))
					ttt<-as.data.frame(newdata[,pdvars]); names(ttt)<-pdvars
					res<-as.data.frame(cbind(ttt,preds))
					return(res)
				}, newdata=newdata,pdvars=pdvars)
	}else if(type=="glm"){
		resdf<-ldply(reslst$models,function(mm,newdata,pdvars){
					mdl<-mm
					preds<-predict(mdl,newdata=newdata,se.fit=TRUE)
					tdf<-data.frame(fit=preds$fit,sefit=preds$se.fit)
					tdf$lwr<-tdf$fit-(1.96*tdf$sefit)
					tdf$upr<-tdf$fit+(1.96*tdf$sefit)
					ttt<-as.data.frame(newdata[,pdvars]); names(ttt)<-pdvars
					res<-as.data.frame(cbind(ttt,tdf))
					return(res)
				}, newdata=newdata,pdvars=pdvars)
	}else{}
	if(NROW(pdvars)==1){
		names(resdf)<-gsub(pdvars,"predictor",names(resdf))
		pddata<-resdf %>% group_by(predictor) %>% dplyr::summarise(predicted=mean(fit),lcl=min(lwr),ucl=max(upr))
		names(pddata)<-gsub("predictor",pdvars,names(pddata))
	}else{
		names(resdf)<-gsub(pdvars[1],"predictor1",names(resdf))
		names(resdf)<-gsub(pdvars[2],"predictor2",names(resdf))
		pddata<-resdf %>% group_by(predictor1,predictor2) %>% dplyr::summarise(predicted=mean(fit),lcl=min(lwr),ucl=max(upr))
		names(pddata)<-gsub("predictor1",pdvars[1],names(pddata))
		names(pddata)<-gsub("predictor2",pdvars[2],names(pddata))
	}
	
	return(as.data.frame(pddata))
}

## This function computes the lielihood ratio test between two models, or between a model and itself minus a named variable
## reslst is the results list from the full model bootstraps
## contr is another list of results from a reduced model's bootstraps
getLRtest<-function(reslst,contr){
	rdf<-ldply(1:NROW(reslst$models),function(mm,reslst,contr){
				mdl2<-reslst$models[[mm]]
				mdl1<-contr$models[[mm]]
				lrtt<-as.data.frame(lrtest(mdl1,mdl2))
				tdf<-data.frame(Ddf=lrtt[2,3],Dloglik=lrtt[1,2]-lrtt[2,2],Chisqv=lrtt[2,4],PrChisq=lrtt[2,5])
				return(tdf)
			},reslst=reslst,contr=contr)
	return(rdf)
}



