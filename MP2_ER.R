
# library(usethis) #you may need to install this using install.packages('usethis')
# use_git_config(user.name = "edenravecca", user.email = "edenravecca@u.boisestate.edu") #your info here
# library(gitcreds) #may need to install this too
# gitcreds_set() #should prompt you for your pat - paste it here # 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

rm ( list = ls() )
gc()

library(tidyverse)
library(sf)
library(terra)
library(tidycensus)
library(tmap)
library(viridis)
library(RColorBrewer)
library(leaflet)

################################.... Load raster data ....################################
# Land cover data for 2016, NLCD, using Terra

nlcd_data_terr <- terra::rast("/opt/data/MP/Ravecca_Eden_MP/NLCD/nlcd_2016_land_cover_l48_20210604.img") 
nlcd_data_terr2 <- terra::rast("/opt/data/MP/Ravecca_Eden_MP/NLCD/nlcd_2016_land_cover_l48_20210604.img")

#################################.... Load PRFA data ....#################################
# Occurrences of PRFA in 2016
# filter occurrences in states of interest, keep relevant columns, rename columns

falcon_csv <- read.csv("/opt/data/MP/Ravecca_Eden_MP/FalconData/PRFA_Data.csv")
states <- c( "Idaho", "Montana", "Wyoming", "Utah", "Nevada", "Arizona", "New Mexico", 
             "Colorado", "California", "Oregon", "Washington" ) # to filter states of interest at once
west_prfa <- falcon_csv %>% 
  dplyr::filter(., stateProvince %in% states) %>% # keep western states
  dplyr::select(., c(occurrenceID, locality, stateProvince, decimalLatitude, decimalLongitude, eventDate,
                     day, month, year)) %>% # keep relevant columns
  dplyr::filter(., year == "2016") %>% # keep only 2016 observations
  dplyr::rename(., c(state = stateProvince, lat = decimalLatitude, lon = decimalLongitude, date = eventDate, ID
                     = occurrenceID)) # change names to make sense
prfa.sf <- st_as_sf(west_prfa, coords = c("lon", "lat"), crs = 4326) # convert prfa gps points to sf object
prfa <- st_transform(prfa.sf, crs= crs(nlcd_data_terr))

################################.... Load states data ....################################
# Poly for all Western States of interest (Rocky Mountain states and west to pacific)

states_data <- st_read("/opt/data/MP/Ravecca_Eden_MP/StatesPoly/tl_2012_us_state/tl_2012_us_state.shp")
w_states_shp <- states_data %>%
  dplyr::select(., c(STUSPS, NAME, geometry)) %>% # keep relevant columns
  dplyr::filter(., NAME %in% states) %>% # keep western states of interest
  dplyr::rename(., c(States = NAME, State_Abbreviation = STUSPS))
all(st_is_valid(w_states_shp)) # make sure geometries are valid; TRUE
w_states_reproj <-  st_transform(w_states_shp, crs = crs(nlcd_data_terr)) # reproject states polys to match NLCD data

###############################.... Load climate data ....#############################
# https://www.fs.fed.us/rm/boise/AWAE/projects/NFS-regional-climate-change-maps/downloads/NationalForestClimateChangeMapsMetadata.pdf

precip_hist_terr <- terra::rast("/opt/data/MP/Ravecca_Eden_MP/NFCC_data/Ppt_annual_historical/Ppt_annual_historical.tif") 
x <- rast(extent= ext(nlcd_data_terr), crs= crs(nlcd_data_terr), resolution= 771)
ppt_hist_reproj <- terra::project(precip_hist_terr, x) # project to NLCD crs

#################################.... Crop Rasters ....################################

w.states.vect <- as(w_states_reproj, "SpatVector") # turn states shapefiles into SpatVector for cropping (in Terra)
nlcd.fullres.crop <- terra::crop(nlcd_data_terr2, w.states.vect) # crop NLCD full res SpatRaster to States poly
ppt.crop <- terra::crop(ppt_hist_reproj, w.states.vect) # crop Precip. SpatRaster to States poly

#########################.... Aggregate & Resample Rasters ....#########################

cls <- c("Unclassified", "Open Water", "Perennial Ice/Snow", "Developed, Open Space", 
         "Developed, Low Intensity", "Developed, Medium Intensity", 
         "Developed High Intensity", "Barren Land", "Deciduous Forest", 
         "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous", 
         "Hay/Pasture", "Cultivated Crops", "Woody Wetlands", "Emergent Herbaceous Wetlands")

nlcd.cats <- terra::cats(nlcd_data_terr2, layer=1) # define nlcd categories
nlcd_modeval <- terra::aggregate(nlcd.fullres.crop, fact=26, fun="modal") # aggregate to 780 m resolution and assign modal value from NLCD classes
m <- rast(extent= ext(w.states.vect), crs= crs(nlcd.fullres.crop), resolution= 771) # blank raster for template
nlcd_771_mode <- terra::resample(nlcd_modeval, m, method= "near", factors=TRUE) # resample to 771 m resolution and assign NLCD value based on "nearest neighbor"
levels(nlcd_771_mode) <- nlcd.cats # nlcd Land Cover classes as levels
names(nlcd_771_mode) <- "NLCD.Land.Cover.Class"
is.factor(nlcd_771_mode)
head(nlcd_771_mode)

##########################.... Create Population Summary ....##########################

census_api_key("ccd5de8eb03adae16092bc5e1c6c40cb104f1e97", install = TRUE, overwrite = TRUE)
readRenviron("~/.Renviron")
options(tigris_use_cache = TRUE)
w.census <- tidycensus::get_estimates(geography = "county",
                                      product = "population",
                                      state = c( "ID", "MT", "WY", "UT", "NV", "AZ", "NM", "CO", "OR", "WA", "CA" ), 
                                      year = 2016,
                                      key = key,
                                      geometry = TRUE,
                                      time_series = T) %>%
  st_transform(., st_crs(w_states_reproj)) %>%
  filter(., variable == "POP") %>%
  spread(., variable, value)

cp.summary <- w.census %>% 
  group_by(NAME) %>% 
  summarize(., avg_pop = mean(POP)) %>%
  separate(NAME,c("county","state"),sep=",")

################################.... Spatial Joins ....################################

cp_prfa <- st_join(prfa, cp.summary, by = "state", left=TRUE)
cp_prfa.vect <- as(cp_prfa, "SpatVector") # turn county population & prfa data into SpatVector for rasterize (in Terra)

################################.... Extract Values ....###############################

ppt_values <- terra::extract(ppt.crop, cp_prfa.vect, fun= mean, method= "bilinear", na.rm=T)
lc_values <- terra::extract(nlcd_771_mode, cp_prfa.vect, fun= max, method= "simple", na.rm=T)
ppt_lc_val <- left_join(ppt_values, lc_values, by= "ID")


tmaptools::palette_explorer()

tm_shape(cp.summary) + 
  tm_polygons(col = "avg_pop",  border.col = "gray", legend.show=TRUE, n=8, style= "cont",
              title = "Population by County (million)") + 
  tm_add_legend(type= "line", lwd = 1, lty= 1, col= "dark gray", title="County Boundaries") +
  tm_legend(outside = TRUE) +
  tm_shape(prfa) + 
  tm_dots(col= "blue", size= 0.01, title= "Prairie Falcons", legend.show=TRUE) +
  tm_add_legend(type= "symbol", size = 0.2, shape = 19, col= "blue", title="Prairie Falcons") +
  tm_legend(outside = TRUE) +
  tm_shape(w_states_reproj) +
  tm_polygons(col = NA, border.col = "black", alpha= 0, title= "Western US States") +
  tm_add_legend(type= "line", lwd = 1, lty= 1, col= "black", title="Western US States") +
  tm_legend(outside = TRUE) +
  tm_layout(main.title = "Prairie Falcon Occurrences in Western US", main.title.size = 1) +
  tm_compass(north = 0, text.size = 0.8, show.labels = 1, 
             cardinal.directions = c("N", "E", "S", "W"), lwd = 1,
             position=c("left", "top"))

brew.ygb <- brewer.pal(8, "YlGnBu")
tmap_mode("plot")
tm_shape(ppt.crop) +
  tm_raster("Ppt_annual_historical", palette= brew.ygb, n=8, legend.show=TRUE, 
            title = "Annual Historical Precipitation (mm)") +
  tm_legend(outside = TRUE) +
  tm_shape(prfa) + 
  tm_dots(col= "black", size= 0.01, title= "Prairie Falcons", legend.show=TRUE) +
  tm_add_legend(type= "symbol", size = 0.2, shape = 19, col= "black", title="Prairie Falcons") +
  tm_legend(outside = TRUE) +
  tm_shape(w_states_reproj) +
  tm_polygons(col = NA, border.col = "dark gray", alpha= 0, title= "Western US States") +
  tm_add_legend(type= "line", lwd = 1, lty= 1, col= "dark gray", title="Western US States") +
  tm_legend(outside = TRUE) +
  tm_layout(main.title = "Prairie Falcon Occurrences in Western US", main.title.size = 1) +
  tm_compass(north = 0, text.size = 0.8, show.labels = 1, 
             cardinal.directions = c("N", "E", "S", "W"), lwd = 1,
             position=c("left", "top"))

brew.set3 <- brewer.pal(7, "Set3")
tm_shape(nlcd_771_mode) +
  tm_raster("NLCD.Land.Cover.Class", palette= brew.set3, n=17, style= "cat", legend.show=TRUE, 
            title = "NLCD Land Cover Types") +
  tm_legend(outside = TRUE) +
  tm_shape(prfa) + 
  tm_dots(col= "navy", size= 0.01, title= "Prairie Falcons", legend.show=TRUE) +
  tm_add_legend(type= "symbol", size = 0.2, shape = 19, col= "navy", title="Prairie Falcons") +
  tm_legend(outside = TRUE) +
  tm_shape(w_states_reproj) +
  tm_polygons(col = NA, border.col = "black", alpha= 0, title= "Western US States") +
  tm_add_legend(type= "line", lwd = 1, lty= 1, col= "black", title="Western US States") +
  tm_legend(outside = TRUE) +
  tm_layout(main.title = "Prairie Falcon Occurrences in Western US", main.title.size = 1) +
  tm_compass(north = 0, text.size = 0.8, show.labels = 1, 
             cardinal.directions = c("N", "E", "S", "W"), lwd = 1,
             position=c("left", "top"))

































































































