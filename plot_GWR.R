
library(raster)
library(leaflet)
library(htmlwidgets)

setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/PhD_DG/to federico")

stations_PM25 <- read.csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/WRF_Chem/validation/PM25_stations_UAE.csv")

# load raster for in situ data
PM25_in_situ_January <- raster("in_situ_kriging_UAE_1KM.tif")
plot(PM25_in_situ_January)
res(PM25_in_situ_January)

# load raster for satellite data converted to PM2.5
PM25_Sat_January <- raster("AOD_mean_jan_1km.tif")
plot(PM25_Sat_January)
res(PM25_Sat_January)

# load raster for satellite data converted to PM2.5
AOD_2016_08_26 <- raster("2016-08-26.tif")
plot(AOD_2016_08_26)
res(AOD_2016_08_26)


# load raster for urban fraction 
URBAN_UAE <- raster("urban_fraction.tif")
plot(URBAN_UAE)


# load raster for desert fraction
DESERT_UAE <- raster("desert_fraction.tif")
plot(DESERT_UAE)

# load Elevation difference ED
ED_UAE <- raster("ED_raster_1km.tif")
plot(ED_UAE)


# load Organic Matter raster for AOD
OM_AOD <- raster("OM_raster_1km.tif")
plot(OM_AOD)

# load Black Carbon raster for AOD
OM_AOD_ECMWF <- raster("OM_raster_1km.tif")
plot(OM_AOD_ECMWF)


# load SO4
SO4_AOD_ECMWF <- raster("SO4_ECMWF_1km.tif")
plot(SO4_AOD_ECMWF)

# load DUST 
DUST_AOD_ECMWF <- raster("DUST_raster_1km.tif")
plot(DUST_AOD_ECMWF)


# load raster for adjusted PM2.5
ADJ_PM25_SATELLITE <- raster("corrected_sat_1KM_LU_ED_DU_SO4.tif")
plot(ADJ_PM25_SATELLITE)

# load raster for adjusted PM2.5
BIAS <- raster("BIAS_predicted.tif")
plot(BIAS)

# load standard error SE
SE_regression <- raster("SE_regression.tif")
plot(SE_regression)

# load R2 regression\
R2_regression <- raster("R2_regression.tif")
plot(R2_regression)


########################################################################
### plots ##############################################################

MIN_PAL <- 15
MAX_PAL <- 70

pal <- colorNumeric(c("#9999FF", "#FFFF00", "#FF0000", "#ff8000"),
                    c(MIN_PAL, MAX_PAL),na.color = "transparent")

pal_August <- colorNumeric(c("#9999FF", "#FFFF00", "#FF0000", "#ff8000"),
                    getValues(AOD_2016_08_26),na.color = "transparent")

pal_LU_URB <- colorNumeric(c("#0000ff", "#ffff00", "#ff0000"),
                           getValues(URBAN_UAE),na.color = "transparent")

pal_LU_DESERT <- colorNumeric(c("#0000ff", "#ffff00", "#ff0000"),
                           getValues(DESERT_UAE),na.color = "transparent")

pal_DUST <- colorNumeric(c("#0000ff", "#ffff00", "#ff0000"),
                              getValues(DUST_AOD_ECMWF),na.color = "transparent")

pal_ED <- colorNumeric(c("#0000ff", "#ffff00", "#ff0000"),
                         getValues(ED_UAE),na.color = "transparent")

pal_OM <- colorNumeric(c("#0000ff", "#ffff00", "#ff0000"),
                       getValues(OM_AOD_ECMWF),na.color = "transparent")

pal_SO4 <- colorNumeric(c("#0000ff", "#ffff00", "#ff0000"),
                       getValues(SO4_AOD_ECMWF),na.color = "transparent")

pal_SE <- colorNumeric(c("#0000ff", "#ffff00", "#ff0000"),
                        getValues(SE_regression),na.color = "transparent")

pal_R2 <- colorNumeric(c("#0000ff", "#ffff00", "#ff0000"),
                        getValues(R2_regression),na.color = "transparent")

pal_BIAS <- colorNumeric(c("#0000ff", "#ffff00", "#ff0000"),
                       getValues(BIAS),na.color = "transparent")


map <- leaflet() %>%
  addTiles(group = "OSM (default)") %>%
  addProviderTiles("OpenStreetMap.Mapnik", group = "Road map") %>%
  addProviderTiles("Thunderforest.Landscape", group = "Topographical") %>%
  addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
  addProviderTiles("Stamen.TonerLite", group = "Toner Lite") %>%
  
  addCircleMarkers(data = stations_PM25,
  lng = ~ Longitude, lat = ~ Latitude,
  radius = 7, stroke = FALSE, fillOpacity = 0.5, popup = ~ Site,
  group = "sites") %>%
  
  # addMarkers(data = stations_PM25, lng = ~ Longitude, lat = ~ Latitude,
  #            popup = ~ Site, group = "Sites") %>%
  addRasterImage(AOD_2016_08_26, colors = pal_August, opacity = 0.4,
                 group = "raster_August") %>%
  addRasterImage(PM25_in_situ_January, colors = pal, opacity = 0.5,
                 group = "in situ PM25") %>%
  addRasterImage(PM25_Sat_January, colors = pal, opacity = 0.5,
                 group = "PM25 AOD 1km") %>%
  addRasterImage(URBAN_UAE, colors = pal_LU_URB, opacity = 0.5,
                 group = "URBAN_UAE") %>%
  addRasterImage(DESERT_UAE, colors = pal_LU_DESERT, opacity = 0.5,
                 group = "DESERT_UAE") %>%
  addRasterImage(ADJ_PM25_SATELLITE, colors = pal, opacity = 0.5,
                 group = "PM25 AOD ADJ 1km") %>%
  addRasterImage(DUST_AOD_ECMWF, colors = pal_DUST, opacity = 0.5,
                 group = "DUST") %>%
  addRasterImage(OM_AOD_ECMWF, colors = pal_BC, opacity = 0.5,
                 group = "OM") %>%
  addRasterImage(SO4_AOD_ECMWF, colors = pal_SO4, opacity = 0.5,
                 group = "SO4") %>%
  addRasterImage(ED_UAE, colors = pal_ED, opacity = 0.5,
                 group = "ED_UAE") %>%
  addRasterImage(SE_regression, colors = pal_SE, opacity = 0.5,
                 group = "SE_UAE") %>%
  addRasterImage(R2_regression, colors = pal_R2, opacity = 0.5,
                 group = "R2_UAE") %>%
  addRasterImage(BIAS, colors = pal_BIAS, opacity = 0.5,
                 group = "BIAS") %>%
  addLegend("bottomright", pal = pal, values = c(MIN_PAL, MAX_PAL), # values = getValues(BIAS), , #  values = getValues(R2_regression), # ,  , 
            title = "<br><strong>PM2.5</strong>",
       #  title = "<br><strong>AOD</strong>",
            labFormat = labelFormat(prefix = ""),
            opacity = 0.6) %>%
  
  addLayersControl(
    baseGroups = c("Road map", "Topographical", "Satellite", "Toner Lite"),
    overlayGroups = c("BIAS",  "SE_UAE","R2_UAE", "SO4", "OM",  "ED_UAE", "DUST", "sites", "raster_August", "in situ PM25", "PM25 AOD 1km", "PM25 AOD ADJ 1km", "URBAN_UAE", "DESERT_UAE"),
    options = layersControlOptions(collapsed = TRUE)) %>%
  hideGroup(c("BIAS",  "SE_UAE", "R2_UAE",  "SO4", "OM", "ED_UAE", "sites",  "DUST", "raster_August",  "URBAN_UAE", "in situ PM25","PM25 AOD 1km", "DESERT_UAE")) 

map



# save map
saveWidget(map, paste0("ADJ_AOD_PM25_1KM_LU_ED_DU_SO4.html"), selfcontained = FALSE)
# saveWidget(map, paste0("Example_Stations_PM25_Sat.html"), selfcontained = FALSE)
