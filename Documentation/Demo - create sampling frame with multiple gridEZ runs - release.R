# Prevent the sp package from calling rgdal
options("sp_evolution_status" = 2)

# Load packages 
library(sf)
library(sp)
library(raster)
library(parallel)
library(igraph)
library(tidyverse)

# Section 1: prep and run gridEZ ----

## Setup ----

# GridEZ function 
source("Q:/GridPopSurvey Frames/2.CODE-gridEZ-master/gridEZ_fn_public_release_v3.R")

# Output folder
output_dir <- "Q:/GridPopSurvey Frames/frame-Demo_ESW/Output/"

# Buffered shapefile folder
buffer_dir <- "Q:/GridPopSurvey Frames/frame-Demo_ESW/Dissolved_Layers/"

# Gridded population raster file path 
gridded_pop_path <- "Q:/GridPopSurvey Frames/3.DATA-global-GHSL/GHS_POP_E2030_GLOBE_R2023A_54009_100_V1_0/GHS_POP_E2030_GLOBE_R2023A_54009_100_V1_0.tif"

# Gridded settlement raster file path 
gridded_sett_path <- "Q:/GridPopSurvey Frames/3.DATA-global-GHSL/GHS_SMOD_E2030_GLOBE_R2023A_54009_1000_V1_0/GHS_SMOD_E2030_GLOBE_R2023A_54009_1000_V1_0.tif"

if (!dir.exists(output_dir)){dir.create(output_dir)}

# Specify number of cores to use during parallel processing and type of parallel processing
ncores <- detectCores() - 1
par_type <- "PSOCK"     # parallel processing type - "PSOCK" for windows

# Projections
proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
proj.moll <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

## Data ----

# Load population data
pop_full <- raster(gridded_pop_path)

# Load settlement layer
settlement_full <- raster(gridded_sett_path)

# Data frame to define gridEZ loops 
loop_var <- "GID_1"
loop_df <- data.frame(
  loop_id = c("SWZ.1_1", "SWZ.2_1", "SWZ.3_1", "SWZ.4_1")
)

errs <- NULL # track any failed loops
for(i in 1:nrow(loop_df)){

  # Set run ID
  loop_runid <- paste0("_", loop_var, "_", loop_df$loop_id[i])
  
  if (
    !file.exists(paste0(output_dir, "EZ_raster_master", loop_runid, ".tif")) & 
    !dir.exists(paste0(output_dir, "temp_folder", loop_runid))
  ){
    
    boundaries <- sf::st_read(
      paste0(
        buffer_dir, "Dissolved_Merged_GID_1_", loop_df$loop_id[i], ".gpkg")) %>% 
      sf::st_transform(., proj) %>% 
      filter(!sf::st_is_empty(.)) %>% 
      sf::as_Spatial()
    
    boundaries@data$fid <- 1:nrow(boundaries@data)
    
    boundaries_moll <- sf::st_read(
      paste0(
        buffer_dir, "Dissolved_Merged_GID_1_", loop_df$loop_id[i], ".gpkg")) %>% 
      sf::st_transform(., proj.moll) %>% 
      filter(!sf::st_is_empty(.)) %>% 
      sf::as_Spatial()

    # Crop population data to boundaries extent
    pop <- crop(pop_full, extent(boundaries_moll))
    pop <- mask(pop, boundaries_moll)
    
    # Crop settlement layer to boundaries extent, transform to same grid as pop
    settlement <- crop(settlement_full, extent(boundaries_moll))
    settlement <- mask(settlement, boundaries_moll)
    settlement <- resample(settlement, pop, method = 'ngb')
    
    pop_gridez <- projectRaster(pop, crs = proj)
    settlement_gridez <- projectRaster(settlement, crs = proj, method = "ngb")

    # Rasterize boundaries to the same grid as population
    boundaries_gridez <- rasterize(boundaries, pop_gridez, field = "fid")
    
    tryCatch({
      gridEZ(
        population_raster = pop_gridez,
        settlement_raster = settlement_gridez,
        strata_raster = boundaries_gridez,
        exclude_unsettled = FALSE,
        using_ghs_smod_pop2015 = TRUE,
        predefined_EZ_size = TRUE,
        EZ_target_size = "medium",
        output_path = output_dir,
        run_ID = loop_runid)
    },
    error = function(x){
      errs <- c(errs, loop_df$loop_id[i])
      assign("errs", errs, envir = globalenv())
    })
    
  } else {next}
}

# Section 2: Clip output rasters ----
# Don't run if there were errors in section 1
if (is.null(errs)){
  
  # Folder with the original split shapefile layers
  boundaries_dir <- "Q:/GridPopSurvey Frames/frame-Demo_ESW/Split_Layers/"
  
  # Data frame to define gridEZ loops 
  loop_var <- "GID_1"
  loop_df <- data.frame(
    loop_id = c("SWZ.1_1", "SWZ.2_1", "SWZ.3_1", "SWZ.4_1")
  )
  
  my_crs <- 2054
  
  for(i in 1:nrow(loop_df)){
    
    loop_runid <- paste0(loop_var, "_", loop_df$loop_id[i])
    
    # Read raster from gridEZ 
    loop_ez <- terra::rast(paste0(output_dir, "EZ_raster_master_", loop_runid, ".tif")) %>% 
      terra::as.polygons() %>% 
      sf::st_as_sf() %>%
      sf::st_transform(., my_crs) 
    
    names(loop_ez) <- c("gridEZ_ID_old", "geometry")
    
    loop_boundary <- sf::st_read(paste0(boundaries_dir, loop_runid, ".gpkg")) %>% 
      summarize()
    
    loop_ez_crop <- sf::st_intersection(loop_ez, loop_boundary)
    orig_nrow <- nrow(loop_ez_crop)
    # If any point/line artifacts hiding in a geometrycollection, extract just
    # the polygons and reattach to main data frame
    if (any(st_geometry_type(loop_ez_crop) %in% "GEOMETRYCOLLECTION")){
      extract_collection <- loop_ez_crop[sf::st_geometry_type(loop_ez_crop) %in% "GEOMETRYCOLLECTION",]
      loop_ez_crop <- loop_ez_crop[!sf::st_geometry_type(loop_ez_crop) %in% "GEOMETRYCOLLECTION",]
      
      extract_poly <- sf::st_collection_extract(extract_collection, "POLYGON")
      
      loop_ez_crop <- rbind(loop_ez_crop, extract_poly)
    }
    
    loop_ez_crop <- sf::st_cast(loop_ez_crop, "MULTIPOLYGON") %>% 
      mutate(gridEZ_ID = row_number(),
             gridEZ_ID = paste0(loop_df$loop_id[i], "_", gridEZ_ID))
    
    sf::write_sf(loop_ez_crop, paste0(output_dir, "R_cropped_", loop_runid, ".gpkg"),
                 append = FALSE)
  }
  
}
