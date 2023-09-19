# .....................................................................
# R script: use gridEZ function to create a sampling frame of gridded
# enumeration zones for Saint Kitts and Nevis
# .....................................................................

# This R script contains the following sections:
# - Instructions for downloading Saint Kitts and Nevis data
# - Load R packages
# - Prepare input data for running gridEZ
# - Run gridEZ function to create a sampling frame of gridded enumeration zones
# - Explore results

# 1. Download data ----

# Download Saint Kitts and Nevis population 2020 from here:
# https://www.worldpop.org/geodata/summary?id=6492

# Download Saint Kitts and Nevis GADM administrative boundaries by selecting the
# country 'Saint Kitts and Nevis' and clicking on 'Shapefile', from here:
# https://gadm.org/download_country_v3.html

# Create a folder called 'KNA'

# Add the population file (kna_ppp_2020.tif) to the KNA folder 

# Extract the files from the GADM shapefile (gadm36_KNA_shp.zip) so that you
# have a folder called gadm36_KNA_shp containing the files. Add this
# gadm36_KNA_shp folder to the KNA folder.

# 2. Load R packages ----

# Packages used in the gridEZ function
library(sp)           # Version 1.6-0
library(raster)       # Version 3.6-20
library(parallel)     # Version 4.2.3 (base package)
library(igraph)       # Version 1.5.1

# Additional packages used in this script (sf used to load a shapefile;
# data.table used to produce summary statistics - neither are needed to run
# gridEZ)
library(sf)           # Version 1.0-12
library(data.table)   # Version 1.14.8

# 3. Prepare input data ----

# Set working directory to the location where your KNA folder and gridEZ code
# are saved 
setwd("C:/Users/clary/Biostat Global Dropbox/Caitlin Clary/CBC Projects/GridEZ")

# Read in population and admin boundaries data
KNApop_raster <- raster("KNA/kna_ppp_2020.tif")
plot(KNApop_raster)

KNAad_poly <- sf::read_sf("KNA/gadm36_KNA_shp/gadm36_KNA_1.shp") |>
  sf::as_Spatial()
plot(KNAad_poly)

# Define a function (gridEZ_rasterize) to convert the GADM boundaries object 
# into a raster

gridEZ_rasterize <- function(strata_shp, population_raster, field){
  r1 <- raster()
  crs(r1) <- crs(population_raster)
  extent(r1) <- extent(population_raster)
  res(r1) <- res(population_raster)
  r2 <- raster::rasterize(strata_shp, r1, field)
  return(r2)
}

# Identify the field that provides unique IDs for the admin units of the
# boundaries data
KNAad_poly@data  # 'GID_1' looks appropriate

# Sometimes the rasterize process doesn't work with character IDs, so convert
# to factor or use a unique ID variable that is already numeric or factor 
KNAad_poly$GID_1 <- as.factor(KNAad_poly$GID_1)
field_ID <- "GID_1"

# Apply the gridEZ_rasterize function to create an admin boundaries raster
KNAad_raster <- gridEZ_rasterize(KNAad_poly, KNApop_raster, field = field_ID)
plot(KNAad_raster)

# Create a settlement layer (this step is simply for this example, there are a
# number of publicly available settlement datasets that can be used)

# Define urban and rural regions
KNAsettagg_raster <- aggregate(KNApop_raster, 10, fun = max, na.rm = TRUE) # this finds the max pop (per 100m x 100m pixel) in each 1km x 1km square
KNAsettagg_raster[KNAsettagg_raster < 8] <- 1              # squares with max pop (per 100m x 100m pixel) smaller than than 8 are classified as '1' representing rural
KNAsettagg_raster[KNAsettagg_raster >= 8] <- 2             # squares with max pop (per 100m x 100m pixel) large than or equal to 8 are classified as '2' representing urban
KNAsettagg_raster[is.nan(KNAsettagg_raster)] <- NA

# Match the extent of the settlement layer to that of the pop layer
KNAsettagg_raster <- projectRaster(
  KNAsettagg_raster, 
  crs = crs(KNApop_raster), 
  res = res(KNApop_raster), method = 'ngb')
new_extent <- extent(
  c(max(extent(KNAsettagg_raster)[1], extent(KNApop_raster)[1]),
    min(extent(KNAsettagg_raster)[2], extent(KNApop_raster)[2]),
    max(extent(KNAsettagg_raster)[3], extent(KNApop_raster)[3]),
    min(extent(KNAsettagg_raster)[4], extent(KNApop_raster)[4])))
KNAsettagg_raster <- crop(KNAsettagg_raster, new_extent)

# Plot the three input rasters
par(mfrow = c(1,3))
plot(KNApop_raster, main = 'population')
plot(KNAad_raster, main = 'admin units')
plot(KNAsettagg_raster, main = 'settlement type (1=rural, 2=urban)')

# 4. Run gridEZ function ----

# Read in gridEZ function
source('gridEZ_fn_public_release_v3.R')

# Specify number of cores to use during parallel processing 
total_ncores <- detectCores()  # This gives the number of cores for your computing system
ncores <- total_ncores - 1     # Use all but 1 core for the parallel processing
if (ncores == 0){ncores <- 1}

# Specify parallel processing type: "PSOCK" for Windows, "FORK" for Linux/Unix/Mac
par_type <- "PSOCK" 

# Run grid EZ code (expect this to take several minutes; run time will depend on
# the number of cores used)

gridEZ(population_raster = KNApop_raster, 
       settlement_raster = KNAsettagg_raster, 
       strata_raster = KNAad_raster, 
       exclude_unsettled = FALSE, 
       using_ghs_smod_pop2015 = FALSE,
       predefined_EZ_size = TRUE, EZ_target_size = "small", 
       output_path = "KNA/", run_ID = "_KNA_smallEZs_2020")

# 5. Explore results ----

# gridEZ saves the output rasters to file. Open them in GIS software or you can
# read them into R:

EZ_pops <- raster("KNA/EZ_pop_raster_master_KNA_smallEZs_2020.tif")
EZ_IDs <- raster("KNA/EZ_raster_master_KNA_smallEZs_2020.tif")
             
# Create data table for EZ results 
EZ_results <- data.table(
  EZ_ID = EZ_IDs[],
  pop = EZ_pops[],
  sett = KNAsettagg_raster[],
  strata = KNAad_raster[])

cell_count_dt <- EZ_results[, .N, by = EZ_ID]

res_dt <- merge(EZ_results, cell_count_dt, by = 'EZ_ID')

# Delete all NA raster cell locations, i.e. outside of country boundary
res_EZ_dt <- res_dt[!is.na(EZ_ID),] |> unique()

# Total number of EZs in sampling frame
nrow(res_EZ_dt)

# Total EZs by sett type
res_EZ_dt[, .N, by = 'sett']     # sett = 1 for rural and 2 for urban

# Total EZs per admin unit
res_EZ_dt[, .N, by = 'strata'] 

# For all EZs, plot histograms for a) number of cells per EZ, and b) pop count per EZ
par(mfrow = c(1,2))

hist(res_EZ_dt[,N], 
     breaks = 100, 
     xlab = "number of cells per EZ", 
     main = paste("All EZs, n =", nrow(res_EZ_dt)))

hist(res_EZ_dt[,pop], 
     breaks = 100, 
     xlab = "population count per EZ", 
     main = paste("All EZs, n =", nrow(res_EZ_dt)))

# Our target population per EZ was 75 so we expect most EZs to have between 56
# and 94 people (+/- 25% of 75)

# Our maximum cells per EZ (i.e. maximum geographic size) was 100 cells so we
# expect most EZs to not exceed 1.5 x the max number of cells (150 in this
# example). The spike at 100 in the number of cells per EZ histogram show this
# restriction being played out in the results.

# For rural EZs, plot histograms for a) number of cells per EZ, and b) pop count per EZ
par(mfrow = c(1,2))

hist(res_EZ_dt[sett==1, N], 
     breaks=100, 
     xlab = "number of cells per EZ", 
     main = paste("Rural EZs, n =", nrow(res_EZ_dt[sett==1])))

hist(res_EZ_dt[sett==1, pop], 
     breaks=100, 
     xlab = "population count per EZ", 
     main = paste("Rural EZs, n =", nrow(res_EZ_dt[sett==1])))

# For urban EZs, plot histograms for a) number of cells per EZ, and b) pop count per EZ
par(mfrow = c(1,2))

hist(res_EZ_dt[sett==2, N], 
     breaks = 100, 
     xlab = "number of cells per EZ", 
     main = paste("Urban EZs, n =", nrow(res_EZ_dt[sett==2])))

hist(res_EZ_dt[sett==2, pop], 
     breaks = 100, 
     xlab = "population count per EZ", 
     main = paste("Urban EZs, n =", nrow(res_EZ_dt[sett==2])))

# From the urban and rural plots you can see that most EZs with population count
# below 45 are found in rural areas. In the very low population density areas
# the EZs have the maximum geographic size per EZ and therefore have low EZ
# population counts.

par(mfrow = c(1,2))
plot(KNApop_raster, main = 'population per 100m x 100m pixel', cex.main=0.7)
plot(EZ_pops, main = 'EZ population count', cex.main=0.7)  

