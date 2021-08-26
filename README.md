# ContinentalWESEestimates
This repository contains the R markdown notebooks, R scripts, and all datasets used in developing the analyses 
and results presented in the manuscript: "Insights from the first global population estimate of Weddell seals in Antarctica", published in Scicen Advances (2021)
by Michelle LaRue1, Leo Salas, Nadav Nur, David Ainley, Sharon Stammerjohn, Jean Pennycook, Melissa Dozier, Jon Saints, Kostas Stamatiou, Luke Barrington,and Jay Rotella. 

<br/><br/>

### Description of notebook, script, and dataset contents

<br/>

#### Notebooks and dependencies
This project uses the following notebooks in the order shown below, with the datasets indicated therein:
  1 Count Seals From Tags.ipynb, needs:
    * compiledData.RData, which comes from the script compileAllData.R
    * MLR_WeddellSea_counts.RData, which comes from the script addMLR_countsInWeddellSea.R
The Count Seals from Tags Notebook requires the compiled data about “feature tags”, which represent Weddell seals “tagged” on high-resolution images using the crowd-sourcing platform, Tomnod (now called GeoHive). During 2016-2018, in concert with DigitalGlobe, our team selected panchromatic, high-resolution images (QuickBird-02, WorldView-01 and WorldView-02; ~0.4-0.6 m spatial resolution; images acquired across 2009-2017) of the sea ice surrounding the Antarctic coastline and requested assistance from volunteers around the world to identify and then place a marker on top of suspected seals. The datasets in this notebook represent the compilation of all seal tags derived from crowd-sourcing efforts.   CompiledData.RData combines all the crowd-sourcing campaigns, which were necessarily split by geography (e.g., The Ross Sea, the Weddell Sea). The data, originally geoJSON files, are compiled all into one dataset.  
The MLR_WeddellSea_counts dataset was derived from the addMLR_countsInWeddellSea.R script, which simply reads an additional dataset from the Weddell Sea, inspected by MLR. We used these data, too, in the analyses.
               
  2 Correct_countsFromGroundCounts.ipynb, needs:
    * countsCorrectedByHour.Rdata, which originates from notebook Count Seals From Tags.ipynb
To determine a population count for every grid cell around the coastline (grid cells = 5 km x 5 km), we needed to correct the aforementioned tags derived in the CompiledData dataset to account for several things: imperfect detection/counting error; time of day and day of year when the image was taken. This dataset includes the corrected the tags data by hour of the day, due to the fact that Weddell seals experience a well-known diel haulout pattern.  
The notebook uses the reference ground counts to adjust the satellite-based counts. The estimated adjustment function is then used on records for the rest of the continent.
               
  3 Calculate_stableAgeDistribution.ipynb, needs no external data
This notebook determines the likely stable age distribution of all Weddell seals, based off data from the Erebus Bay population (Rotella et al. 2012, Paterson et al. 2018 – see references in the paper for details). We this information is used to correct for the proportions of age classes in the ground truthing colonies, which are different from the stable age distribution, to get an adequate estimate for the local and regional counts.
  4 AttributeGrid_withWESEdata.ipynb, needs:
    * FinalWESEcounts.RData, which originates from Correct_countsFromGroundCounts.ipynb
    * studyarea_points_wNearLandPenguins.RData, which originates from the script changePredictability_DecToOct.R
    * wese5k.RData, which is the data grid with the coordinates of centerpoints in the polar projection we used, provided by Point Blue
    * FastIceGridPoints_weseNoIce_111114.RData, a file of locations with WESE counts but no ice based on the grid's midpoint, provided by Point Blue
FinalWESEcounts is the population estimate per grid cell for the year 2011, that takes into account the tags from crowdsourcing, and then corrects for time of day, imperfect detection, and day of year, and is adjusted based on the ground counts and the bias in age distribution of the ground counts. The wese5k.RData is the grid used to determine the extent of a population around the continent (e.g., the smallest spatial scale for a population of Weddell seals). We attribute it with the WESE abundance estimates (from FinalWESEcounts), biophysical (e.g., fast ice presence, distance to shore – from studyarea_points_wNearLandPenguins) and adjusted values for a few centerpoints with WESE that land on solid ground and a re-positioning within the cell must be done to attribute with biophysical characteristics (from FastIceGridPoints_weseNoIce_111114 – the 111114 is the NIC ice date used to estimate presence and abundance of fast ice).
  
  5 Prepare_WESE_analysisData.ipynb, needs:
    * continentalWESE.RData, which originates from the script changePredictability_DecToOct.R
    * ADPE_colonies_20200416.csv, which was provided by MLR
ADPE Colonies (location and colony size) was derived from Lynch and LaRue (2014)
    * EMPE_colonies_20200327.csv, which was provided by MLR
EMPE colonies (locations and colony size) was derived from Fretwell et al. (2012)
               
  6 ContinentalWESE_regressionAnalyses.ipynb, needs:
    * WESEdata_forContinentalAnalyses.RData, which originates from Prepare_WESE_analysisData.ipynb
    
#### Scripts used to prepare data, or sourced by the notebooks
The notebooks source the following scripts, or these are used to generate preliminary datasets:
  * addMLR_countsInWeddellSea.R is called by MLR_WeddellSea_counts.RData
  * changePredictability_DecToOct.R is needed to update the dataset studyarea_points_wNearLandPenguins.RData used in AttributeGrid_withWESEdata.ipynb. The attribution includes information about the predictability of the fast ice throughout the pupping season and over time in the past 5 years.
  * compileAllData.R produces the starting dataset used in Count Seals From Tags.ipynb
  * countSealsFromTags_functions.R is called by Count Seals From Tags.ipynb, Correct_countsFromGroundCounts.ipynb and AttributeGrid_withWESEdata.ipynb. It provides utility functions called by the notebook.
  * fitRegressionModels_functions.R is called by Prepare_WESE_analysisData.ipynb and ContinentalWESE_regressionAnalyses.ipynb. It provides utility functions called by these notebooks.
  * fastIceCovars_utils.R is called by Prepare_WESE_analysisData.ipynb and provides utility functions used in that notebook.

#### Datasets used
The following datasets are needed to run the notebooks (what is each and what else am I missing?):
  * compiledData.RData - the original data from DG/Maxar
  * MLR_WeddellSea_counts.RData = MLR counts in the Weddell Sea
  * FinalWESEcounts.RData - an intermediate data product
  * countsCorrectedByHour.Rdata - from notebook Count Seals From Tags.ipynb
  * studyarea_points_wNearLandPenguins.RData - the grid points with all ice and geophysical attributes, provided by Point Blue using the FastIce tool
  * wese5k.RData - the grid with coordinates and no additional attribution
  * FastIceGridPoints_weseNoIce_111114.RData - a file of locations with WESE counts but no ice based on the grid's midpoint, provided by Point Blue
  * continentalWESE.RData - an intermediate data product from changePredictability_DecToOct.R
  * ADPE_colonies_20200416.csv - locations and colony sizes of ADPE colonies in 2011, provided by MLR
  * EMPE_colonies_20200327.csv - locations and colony sizes of EMPE colonies in 2011, provided by MLR
  * WESEdata_forContinentalAnalyses.RData - from Prepare_WESE_analysisData.ipynb

