#!/usr/bin/env Rscript

# FedData provides functions for getting GHCN data, 
# and the `pkg_test` function for installing/loading other packages
#install.packages("devtools")
# devtools::install_cran("FedData")
# Python-style argument parsing
library("optparse")

# Set the number of parallel cores
# We'll be using the snow-like functionality with Rmpi
slurm_cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))

# Use optparse for Python-style argument parsing
# Define a set of options; only output file and force_redo for now
option_list = list(
  make_option(c("-o", "--output_dir"),
              type = "character",
              default = "./OUTPUT/", 
              help = "Output directory [default = %default]",
              metavar = "character"),
  make_option(c("-n", "--cores"),
              type = "numeric",
              default = ifelse(is.na(slurm_cores), 2, slurm_cores),
              help = "Number of cores for multicore run [default = %default]"),
  make_option(c("-c", "--clean"),
              action = "store_true",
              default = FALSE,
              help = "Delete output directory and re-run all analyses? [default= %default]")
); 

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

### BEGIN SCRIPT ###

start_time <- Sys.time()
message("asian_niche.R started at ", start_time)

## This is the code for the pan-Asian niche reconstructions

## Load all packages
all_packages <- c("foreach", "doParallel", # Packages for parallel processeing
                  "FedData", # Package for data aquisition
                  "R.utils", "Hmisc", "zoo", "abind", "mgcv", "rgbif", "fields", # Packages offering general utilities
                  "sf", "rgdal", "ncdf4", "raster", "geomapdata", "maptools", "mapproj", # Packages for spatial processing
                  "Bchron", "mclust", # Packages for chronometric analysis
                  "magrittr", "tidyverse", "ggthemes", "purrrlyr", # Packages for tidy code
                  "RColorBrewer", "htmlwidgets", "plotly", "bibtex", "knitcitations") # Plotting and rmarkdown

purrr::walk(all_packages, library, character.only = TRUE)

# Load all functions
list.files("./src",full.names=T) %>%
  purrr::walk(source)

# Force Raster to load large rasters into memory
rasterOptions(chunksize=2e+08,maxmemory=2e+09)

##### SET PARAMETERS #####

# Delete the output directory if requested, then create it
opt$output_dir <- ifelse(stringr::str_sub(opt$output_dir,-1) != "/",
                         stringr::str_c(opt$output_dir,"/"),
                         opt$output_dir) 
if(opt$clean) unlink(opt$output_dir,
                     recursive = TRUE,
                     force = TRUE)
dir.create(opt$output_dir,
           showWarnings = FALSE,
           recursive = TRUE)

# a function that builds output paths
out <- function(...){
  stringr::str_c(opt$output_dir,...)
}

# Create a directory for writing tables
c("DATA","MODELS","RECONS","FIGURES","TABLES","SITE_DENSITIES") %>%
  purrr::walk(~ dir.create(out(.),
                           showWarnings = FALSE,
                           recursive = TRUE))

# A function to set an objects class
set_class <- function(x, classes){
  class(x) <- classes
  return(x)
}

crop_custom <- function(x, y) {
  x.sp <- as(x, "Spatial")
  y.sp <- as(y, "Spatial")
  x.sp.crop <- raster::crop(x.sp, y.sp)
  st_as_sf(x.sp.crop)
}

#Define the study region
ASIA_poly <-
  extent(55,125,5,60) %>% ## THE REAL EURASIAN POLYGON
  # extent(66,72,37,43) %>% ## FOR TESTING PURPOSES ONLY
  FedData::polygon_from_extent("+proj=longlat +ellps=GRS80")

# Set the calibration period for paleoclimate reconstructions
calibration.years <- 1961:1990

# Set this to your google maps elevation api key
# https://developers.google.com/maps/documentation/elevation/start
google_maps_elevation_api_key = "AIzaSyDi4YVDZPt6uH1C1vF8YRpbp1mxqsWbi5M"

# Report the number of cores for parallel processing
message("Number of cores set to SLURM_NTASKS_PER_NODE = ",opt$cores)

# A ggplot2 theme for Nature Publishing Group
nature_theme <- ggplot2::theme(panel.grid.major = element_line(size = 0.5, color = "grey"),
                               axis.line = element_line(size = 0.7, color = "black"),
                               legend.position = c(0.85, 0.7),
                               text = element_text(size = 14))

##### END SET PARAMETERS #####



##### LOAD ETOPO5 GRID #####

message("Preparing the ETOPO5 grid-aligned dataset")
time_check <-  Sys.time()
if(!opt$clean & file.exists(out("DATA/ASIA_rast_etopo5.tif"))){
  ASIA_rast_etopo5 <- raster(out("DATA/ASIA_rast_etopo5.tif"))
}else{

  # Get the Natural Earth country lakes data
  dir.create(out("DATA/NaturalEarth/"), showWarnings = FALSE, recursive = TRUE)
  FedData::download_data(url="http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_lakes.zip",destdir=out("DATA/NaturalEarth/"))
  unzip(out("DATA/NaturalEarth/ne_10m_lakes.zip"), exdir=out("DATA/NaturalEarth/ne_10m_lakes/"))
  ne_10m_lakes <- out("DATA/NaturalEarth/ne_10m_lakes/ne_10m_lakes.shp") %>%
    sf::st_read() %>%
    crop_custom(ASIA_poly)

  # Get the Natural Earth country boundaries data
  dir.create(out("DATA/NaturalEarth/ne_10m_admin_0_countries_lakes/"), showWarnings = FALSE, recursive = TRUE)
  FedData::download_data(url="http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries_lakes.zip",destdir=out("DATA/NaturalEarth/"))
  unzip(out("DATA/NaturalEarth/ne_10m_admin_0_countries_lakes.zip"), exdir=out("DATA/NaturalEarth/ne_10m_admin_0_countries_lakes/"))
  ne_10m_admin_0_countries_lakes <- out("DATA/NaturalEarth/ne_10m_admin_0_countries_lakes/ne_10m_admin_0_countries_lakes.shp") %>%
    sf::st_read() %>%
    crop_custom(ASIA_poly) %>%
    dplyr::filter(!(NAME %in% c("Indonesia",
                                "Scarborough Reef",
                                "Malaysia",
                                "Philippines",
                                "Maldives",
                                "Spratly Is.",
                                "Oman",
                                "United Arab Emirates",
                                "Saudi Arabia")))

  data("ETOPO5")
  ASIA_rast_etopo5 <- ETOPO5 %>%
    t() %>%
    raster(xmn = 0,
           xmx = 360,
           ymn = -90,
           ymx = 90,
           crs = CRS("+proj=longlat +ellps=clrk66 +no_defs")) %>%
    raster::crop(y = sp::spTransform(ASIA_poly, CRSobj = raster::projection(.))) %>%
    raster::mask(ne_10m_admin_0_countries_lakes %>%
                   sf::st_transform("+proj=longlat +ellps=clrk66 +no_defs") %>%
                   as("Spatial")) %>%
    raster::mask(ne_10m_lakes %>%
                   sf::st_transform("+proj=longlat +ellps=clrk66 +no_defs") %>%
                   as("Spatial"),
                 inverse = TRUE) %T>%
    writeRaster(filename = out("DATA/ASIA_rast_etopo5.tif"),
                datatype="INT2S",
                options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND", "PHOTOMETRIC=MINISWHITE"),
                overwrite=T,
                setStatistics=FALSE)
  rm(ETOPO5)
  # unlink(out("DATA/NaturalEarth/ne_10m_lakes/"), recursive = T, force = T)
  # unlink(out("DATA/NaturalEarth/ne_10m_admin_0_countries_lakes/"), recursive = T, force = T)
}

countries <- out("DATA/NaturalEarth/ne_10m_admin_0_countries_lakes/ne_10m_admin_0_countries_lakes.shp") %>%
  sf::st_read()  %>%
  crop_custom(ASIA_poly) %>%
  dplyr::filter(!(NAME %in% c("Indonesia",
                              "Scarborough Reef",
                              "Malaysia",
                              "Philippines",
                              "Maldives",
                              "Spratly Is.",
                              "Oman",
                              "United Arab Emirates",
                              "Saudi Arabia"))) %>%
  dplyr::select(NAME, geometry)

message("ETOPO5 grid-aligned dataset preparation complete: ", capture.output(Sys.time() - time_check))

##### END LOAD ETOPO5 GRID #####



##### PREPARE THE GHCN DATA #####

## Downloads and cleans daily climate records from the Global Historical Climate Database.
message("Preparing the daily climate records from the Global Historical Climate Database")
time_check <-  Sys.time()
GHCN.data.final <- prepare_ghcn(region = ASIA_poly,
                                label = "ASIA_poly",
                                calibration.years = calibration.years,
                                google_maps_elevation_api_key = google_maps_elevation_api_key,
                                force.redo = opt$clean)
message("GHCN preparation complete: ", capture.output(Sys.time() - time_check))
## An example of plotting the GHCN data
# climate_plotter(data = GHCN.data.final, station = "CHM00051334", element = "TMIN")

##### END PREPARE THE GHCN DATA #####



##### PREPARE THE MARCOTT DATA #####

# Run the script that transforms the Marcott et al. 2013 data into standard scores.
message("Preparing the Marcott et al. 2013 data")
time_check <-  Sys.time()
marcott2013 <- prepare_marcott(calibration.years = calibration.years)
message("Marcott et al. 2013 data preparation complete: ", capture.output(Sys.time() - time_check))

##### END PREPARE THE MARCOTT DATA #####



##### MODULATING CLIMATOLOGY BY MARCOTT SD #####

#### Calculating growing degree days ####
message("Modulating local climatology by Marcott et al. 2013 data")
time_check <-  Sys.time()
# How often to sample GDD, in z-space, for model tuning
# Here, we sample from -20 to 20 SD, at 1 SD interval
sample.points <- -20:20

# Read in data on different crop GDD needs
crop_GDD <- readr::read_csv("./DATA/crops.csv",
                            col_types = readr::cols(
                              cultivar_long = col_character(),
                              cultivar = col_character(),
                              crop_long = col_character(),
                              crop = col_character(),
                              t_base = col_double(),
                              min_gdd = col_integer()
                            )) %>%
  dplyr::filter(crop %in% c("foxtail_millet",
                            "broomcorn_millet",
                            "wheat",
                            "barley",
                            "buckwheat"))

# create the cluster for parallel computation
cl <- makeCluster(opt$cores, type = "PSOCK")
registerDoParallel(cl)

# Transform GHCN data to GDDs of each base, and modulate to Marcott
GDDs <- sort(unique(crop_GDD$t_base))
GHCN.GDD.incremented.sd <- foreach::foreach(base = GDDs,
                                            .packages = c("foreach","magrittr")) %dopar% {

                                              out.list <- foreach::foreach(change = sample.points,
                                                                           .packages = c("foreach","magrittr")) %do% {

                                                                             GHCN.GDDs <- foreach::foreach(station = GHCN.data.final$climatology, .combine = c) %do% {

                                                                               return(sdModulator(data.df = station,
                                                                                                  temp.change.sd = change,
                                                                                                  t.base = base,
                                                                                                  t.cap = 30))

                                                                             }

                                                                             names(GHCN.GDDs) <- names(GHCN.data.final$climatology)

                                                                             return(tibble::tibble(SD_change = change,
                                                                                                   ID = names(GHCN.GDDs),
                                                                                                   GDD = GHCN.GDDs))
                                                                           }
                                              names(out.list) <- sample.points

                                              return(out.list %>%
                                                       dplyr::bind_rows() %>%
                                                       dplyr::left_join(GHCN.data.final$spatial %>%
                                                                          dplyr::as_data_frame() %>%
                                                                          dplyr::rename(x = LONGITUDE, y = LATITUDE), by = "ID"))
                                            }
names(GHCN.GDD.incremented.sd) <- GDDs

# stop the cluster (will free memory)
stopCluster(cl)
message("Modulation of local climatology by Marcott et al. 2013 data complete: ", capture.output(Sys.time() - time_check))

##### END MODULATING CLIMATOLOGY BY MARCOTT SD #####



##### GENERATE CROP NICHE MODELS #####

message("Calculating indicator Krige models")
time_check <-  Sys.time()
# Create a spatialPointsDataFrame of the etopo5 data, and convert to WGS84 ellipsoid
ASIA_rast_etopo5.sp <-
  ASIA_rast_etopo5 %>%
  magrittr::set_names("elevation") %>%
  raster::rasterToPoints(spatial = T) %>%
  sp::spTransform(sp::CRS(raster::projection(GHCN.data.final$spatial)))

# A function that generates the kriging model, then predicts
krige_and_predict <- function(dt){
  model <- fields::mKrig(x = dt[,c("x","y")],
                         y = dt$GDD_thresh,
                         Z = dt$elevation,
                         Covariance = "Exponential",
                         Distance = "rdist.earth")

  prediction <- ASIA_rast_etopo5.sp %>%
    tibble::as_tibble() %>%
    dplyr::mutate(chunk = rep(1:ceiling(nrow(.)/chunk_size),length.out = nrow(.)) %>%
                    sort()
    ) %>%
    dplyr::group_by(chunk) %>%
    dplyr::do(prediction = fields::predict.mKrig(model,
                                                 xnew = .[,c("x","y")],
                                                 Z = .$elevation) %>%
                as.vector()) %$%
    prediction %>%
    unlist()

  return(list(model = model, prediction = prediction))

}

# Calculate gdd kriging models for each crop
gdd_model_files <- out("MODELS/",crop_GDD$cultivar,"_models.rds")

if(opt$clean){
  unlink(out("MODELS"), recursive = TRUE, force = TRUE)
  dir.create(out("MODELS/"), showWarnings = F)
  crop_GDD_run <- crop_GDD
}else{
  dir.create(out("MODELS/"), showWarnings = F)
  crop_GDD_run <- crop_GDD[!file.exists(gdd_model_files),]
}

if(nrow(crop_GDD_run) == 0){
  message("All indicator Krige models have already been calculated. Continuing.")
}else message("Calculating indicator Krige models for ",
              nrow(crop_GDD_run),
              " cultivars:\n",
              paste0(capture.output(crop_GDD_run), collapse = "\n"))

# A function to reduce the size of a loess model by half
skinny.loess <- function(x){
  x[c("fitted",
      "residuals",
      "enp",
      "one.delta",
      "two.delta",
      "trace.hat",
      "call",
      "terms",
      "xnames")] <- NULL
  return(x)
}

# A function of correct the indication predictions and estimate a smooth
# monotonic function
# This first uses isotonic regression, then loess smoothing with a degree of 1
smooth.preds <- function(y){
  y[y<0] <- 0
  y[y>1] <- 1
  y <- loess(isoreg(y~sample.points)$yf~sample.points, span=0.1, degree=1) %>%
    skinny.loess()
  return(y)
}

if(nrow(crop_GDD_run) > 0){
  # create the cluster for parallel computation
  cl <- makeCluster(min(opt$cores,nrow(crop_GDD_run)), type = "PSOCK")
  registerDoParallel(cl)

  options(dplyr.show_progress = FALSE)
  chunk_size <- 10000

  gdd.models <- foreach::foreach(crop = 1:nrow(crop_GDD_run),
                                 .packages = c("fields",
                                               "dplyr",
                                               "magrittr",
                                               "foreach",
                                               "doParallel",
                                               "readr"),
                                 .export = c("sample.points")) %dopar% {

                                   # Threshold for indicator kriging
                                   GHCN.GDD.incremented.sd[[as.character(crop_GDD_run[crop,"t_base"])]] %>%
                                     dplyr::mutate(GDD_thresh = {GDD >= as.numeric(crop_GDD_run[crop,"min_gdd"])}) %>%
                                     dplyr::group_by(SD_change) %>%
                                     dplyr::do(out_preds = krige_and_predict(.)) %$%
                                     out_preds %>%
                                     sapply("[[","prediction") %>%
                                     apply(1,smooth.preds) %>%
                                     tibble::tibble(model = .) %>%
                                     maptools::spCbind(ASIA_rast_etopo5.sp,.) %>%
                                     readr::write_rds(out("MODELS/",crop_GDD_run[crop,"cultivar"],"_models.rds"), compress = "xz")

                                   return(out("MODELS/",crop_GDD_run[crop,"cultivar"],"_models.rds"))
                                 }

  # stop the cluster (will free memory)
  stopCluster(cl)

}

rm(ASIA_rast_etopo5,
   ASIA_rast_etopo5.sp,
   GHCN.data.final,
   GHCN.GDD.incremented.sd)

message("Calculation indicator Krige models complete: ", capture.output(Sys.time() - time_check))

##### END GENERATE CROP NICHE MODELS #####



##### PREDICT CROP NICHE THROUGH TIME #####

## Predicting crop niche from smoothed Krige models
# Calculate niches for each crop using the Marcott et al. 2013.
# create the cluster for parallel computation
message("Generating niche reconstructions")
time_check <-  Sys.time()
if(opt$clean){
  unlink(out("RECONS"), recursive = TRUE, force = TRUE)
}
dir.create(out("RECONS/"), showWarnings = F)

purrr::walk(crop_GDD$cultivar, function(crop){

  if(!file.exists(out("MODELS/",crop,"_models.rds"))) stop("Models for ",
                                                           crop,
                                                           " are missing! Aborting.")
  Zs <- c("Z_Lower","Z","Z_Upper")
  Zs %<>%
    magrittr::extract({
      Zs %>%
        out("RECONS/",crop,"_",.,".nc") %>%
        file.exists() %>%
        magrittr::not()
    })

  if(length(Zs)==0) return()

  crop.models <- readr::read_rds(out("MODELS/",crop,"_models.rds"))

  purrr::walk(Zs, .f = function(z){
    suppressWarnings(
      crop.models@data %$%
        model %>%
        purrr::map(.f = function(x){
          x %>%
            predict(newdata = marcott2013[[z]]) %>%
            magrittr::multiply_by(100) %>%
            round()
        }) %>%
        do.call(rbind, .) %>%
        tibble::as_tibble() %>%
        magrittr::set_colnames(marcott2013$YearBP) %>%
        new("SpatialPointsDataFrame",
            data = .,
            coords.nrs = crop.models@coords.nrs,
            coords = crop.models@coords,
            bbox = crop.models@bbox,
            proj4string = crop.models@proj4string) %>%
        as("SpatialPixelsDataFrame") %>%
        raster::brick() %>%
        raster::setZ(marcott2013$YearBP, name="Years BP") %>%
        raster::writeRaster(out("RECONS/",crop,"_",z,".nc"),
                            format = "CDF",
                            datatype = "INT2S",
                            varname = "niche_probability",
                            varunit = "unitless",
                            longname = "Probability of being in the crop niche x 100",
                            xname = "Longitude",
                            yname = "Latitude",
                            zname = "Years BP",
                            zunit = "Years BP",
                            compression = 9,
                            overwrite = TRUE)
    )
  })
})

message("Generation of niche reconstructions complete: ", capture.output(Sys.time() - time_check))

##### END PREDICT CROP NICHE THROUGH TIME #####



##### BEGIN COMBINE LIKE CROP NICHES #####
## Combining crop niches from similar crops by taking an arithmatic mean
message("Combining like crop niches")
time_check <-  Sys.time()

combine_varieties <- function(x, file_tail){
  if(file.exists(out("RECONS/All_",x$crop[[1]],"_",file_tail,".nc"))) return(NULL)
  n_crops <- length(x$crop)
  purrr::map(x$cultivar, function(cultivar){
    out("RECONS/",cultivar,"_",file_tail,".nc") %>%
      raster::brick() %>%
      raster::getValues()
  }) %>%
    Reduce(f = "+", x = .) %>%
    magrittr::divide_by(n_crops) %>%
    round() %>%
    raster::setValues(raster::brick(out("RECONS/",x$cultivar[[1]],"_",file_tail,".nc")), .) %>%
    raster::writeRaster(out("RECONS/All_",x$crop[[1]],"_",file_tail,".nc"),
                        format = "CDF",
                        datatype = "INT2S",
                        varname = "niche_probability",
                        varunit = "unitless",
                        longname = "Probability of being in the crop niche x 100",
                        xname = "Longitude",
                        yname = "Latitude",
                        zname = "Years BP",
                        zunit = "Years BP",
                        compression = 9,
                        overwrite = TRUE)
  return(NULL)
}

# create the cluster for parallel computation
# cl <- makeCluster(min(opt$cores,
#                       crop_GDD %$%
#                         crop %>%
#                         unique() %>%
#                         length()),
#                   type = "PSOCK")
# registerDoParallel(cl)
#
# # Get mean niche for Marcott predictions
# foreach::foreach(crop = crop_GDD %>%
#                    split(as.factor(crop_GDD$crop)),
#                  .packages = c("magrittr",
#                                "foreach"),
#                  .combine = c) %dopar% {
purrr::walk(crop_GDD %>%
              split(as.factor(crop_GDD$crop)),
            function(crop){
              combine_varieties(crop, file_tail = "Z") # Get mean niche for Marcott
              combine_varieties(crop, file_tail = "Z_Upper") # Get mean niche for Marcott upper CI
              combine_varieties(crop, file_tail = "Z_Lower") # Get mean niche for Marcott lower CI
            })

# stop the cluster (will free memory)
# stopCluster(cl)

message("Combining like crop niches complete: ", capture.output(Sys.time() - time_check))

##### END COMBINE LIKE CROP NICHES #####



##### BEGIN PLOT CULTIVAR NICHE THROUGH TIME #####

## Plotting cultivar niche
# create the cluster for parallel computation
message("Plotting cultivar niche reconstructions")
time_check <-  Sys.time()

gdd.recons <- foreach::foreach(n = 1:nrow(crop_GDD),
                               .combine = c) %do% {

                                 cultivar <- crop_GDD[n,]$cultivar

                                 title <- stringr::str_c(crop_GDD[n,]$cultivar_long,
                                                         " — T_base: ",
                                                         crop_GDD[n,]$t_base,
                                                         "°C, Required GDD: ",
                                                         crop_GDD[n,]$min_gdd)

                                 if(file.exists(out("FIGURES/",cultivar,".pdf")) & file.exists(out("FIGURES/",cultivar,".mov")))
                                   return(out("FIGURES/",cultivar,".pdf"))

                                 rast <- raster::brick(out("RECONS/",cultivar,"_Z.nc")) %>%
                                   magrittr::extract2(which(.@z$`Years BP` > 1000)) %>%
                                   raster:::readAll() %>%
                                   magrittr::divide_by(100) %>%
                                   magrittr::extract2(nlayers(.):1)

                                 rast.lower <- raster::brick(out("RECONS/",cultivar,"_Z_Lower.nc")) %>%
                                   magrittr::extract2(which(.@z$`Years BP` > 1000)) %>%
                                   raster:::readAll() %>%
                                   magrittr::divide_by(100) %>%
                                   magrittr::extract2(nlayers(.):1)

                                 rast.upper <- raster::brick(out("RECONS/",cultivar,"_Z_Upper.nc")) %>%
                                   magrittr::extract2(which(.@z$`Years BP` > 1000)) %>%
                                   raster:::readAll() %>%
                                   magrittr::divide_by(100) %>%
                                   magrittr::extract2(nlayers(.):1)

                                 years <- rast %>%
                                   names() %>%
                                   gsub(pattern = "X",
                                        replacement = "",
                                        x = .) %>%
                                   as.numeric()

                                 pal <- c(rev(colorRampPalette(brewer.pal(9, "Blues")[2:9],
                                                               bias = 2,
                                                               space = "Lab")(76)),
                                          colorRampPalette(brewer.pal(9, "Reds")[2:9],
                                                           bias = 1.5,
                                                           space = "Lab")(26))

                                 if(!file.exists(out("FIGURES/",cultivar,".pdf")))
                                   space_time_plot(the_brick = rast,
                                                   the_brick_lower = rast.lower,
                                                   the_brick_upper = rast.upper,
                                                   out_file = out("FIGURES/",cultivar,".pdf"),
                                                   title = title,
                                                   time = years,
                                                   timelim = c(max(years),min(years)),
                                                   timeaxis =  seq(from = max(years)-500,
                                                                   to = min(years),
                                                                   by = -500),
                                                   timelab = "Years BP",
                                                   zbreaks = seq(0,1,0.01),
                                                   zlab = "Probability of being in niche",
                                                   zaxis = seq(0,1,0.1),
                                                   zcolors = pal
                                   )

                                 if(!file.exists(out("FIGURES/",cultivar,".mov")))
                                   space_time_video(the_brick = rast,
                                                    the_brick_lower = rast.lower,
                                                    the_brick_upper = rast.upper,
                                                    out_file = out("FIGURES/",cultivar,".mov"),
                                                    title = title,
                                                    time = years,
                                                    timelim = c(max(years),min(years)),
                                                    timeaxis =  seq(from = max(years)-500,
                                                                    to = min(years),
                                                                    by = -500),
                                                    timelab = "Years BP",
                                                    zbreaks = seq(0,1,0.01),
                                                    zlab = "Probability of being in niche",
                                                    zaxis = seq(0,1,0.1),
                                                    zcolors = pal
                                   )

                                 return(out("FIGURES/",cultivar,".pdf"))
                               }

message("Plotting of cultivar niche reconstructions complete: ", capture.output(Sys.time() - time_check))

##### END PREDICT CULTIVAR NICHE THROUGH TIME #####



##### BEGIN CHRONOMETRIC ANALYSIS #####
message("Beginning chronometric co-analysis")

# Read in the site/chronometric data
chronometric_data <- readxl::read_excel("DATA/DALPOIMGUEDES_BOCINSKY_2017.xlsx",
                                        sheet = "chronometric_data",
                                        col_types = c("text",
                                                      "numeric",
                                                      "numeric",
                                                      "numeric",
                                                      "logical",
                                                      "text",
                                                      "text",
                                                      "text",
                                                      "logical",
                                                      "numeric",
                                                      "numeric",
                                                      "text",
                                                      "numeric",
                                                      "numeric",
                                                      "text",
                                                      rep("logical",18))
) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(`Exclude?` = ifelse(is.na(`Exclude?`), FALSE, `Exclude?`)) %>%
  dplyr::filter(!is.na(Site),
                !is.na(Longitude),
                !is.na(Latitude),
                !(is.na(`14C age BP`) & is.na(`Age range lower (BP)`)),
                !is.na(`14C date on cereal?`),
                !`Exclude?`) %>%
  dplyr::select(-`Exclude?`, -Notes) %>%
  group_by(Site,Period)


# Filter out sites not in study area
chronometric_names <- names(chronometric_data)
chronometric_data %<>%
  sf::st_as_sf(coords = c("Longitude", "Latitude"), remove = FALSE) %>%
  sf::st_set_crs("+proj=longlat") %>%
  sf::st_intersection(ASIA_poly %>%
                        sf::st_as_sf() %>%
                        sf::st_transform("+proj=longlat") %>%
                        sf::st_geometry()) %>%
  sf::st_intersection(countries %>%
                        sf::st_union() %>%
                        sf::st_transform("+proj=longlat")) %>%
  dplyr::as_data_frame() %>%
  dplyr::select(-geometry) %>%
  magrittr::set_names(chronometric_names)

# Write table of dates
chronometric_data %>%
  dplyr::mutate(`Estimated age range (BP)` = stringr::str_c(`Age range lower (BP)`,"–",`Age range upper (BP)`)) %>%
  dplyr::select(Site,
                Period,
                `Lab sample identifier`,
                Material,
                `14C age BP`,
                `1-sigma uncertainty`,
                `Estimated age range (BP)`,
                Reference) %>%
  # dplyr::filter(!is.na(`14C age BP`)) %>%
  dplyr::arrange(Site, Period) %>%
  readr::write_csv(out("TABLES/sites_dates_raw.csv"))

# Write some stats about the sites data used in this analysis
write_lines(c("Sites data (after spatial and crop filtering)",
              paste0("Number of entries: ",
                     chronometric_data %>%
                       nrow()),
              paste0("Number of sites: ",
                     chronometric_data %$%
                       Site %>%
                       unique() %>%
                       length()),
              paste0("Number of distinct occupations: ",
                     chronometric_data %>%
                       dplyr::select(Site,Period) %>%
                       distinct() %>%
                       nrow()),
              paste0("Number of radiocarbon dates: ",
                     chronometric_data %>%
                       dplyr::filter(!is.na(`14C age BP`)) %>%
                       nrow()),
              paste0("Number of radiocarbon sites: ",
                     chronometric_data %>%
                       dplyr::filter(!is.na(`14C age BP`)) %$%
                       Site %>%
                       unique() %>%
                       length())
              ),
            out("sites.txt"))


# A function to calibrate either 14C dates (using BchronCalibrate), or
# site age estimates (using a flat density estimation)
calibrate_density <- function(x){
  if(!is.na(x$`14C age BP`)){
    out_calib <- Bchron::BchronCalibrate(ages = x$`14C age BP`,
                                         ageSds = x$`1-sigma uncertainty`,
                                         calCurves = 'intcal13') %>%
      magrittr::extract2(1)
    out_calib <- out_calib[c("ageGrid","densities")]
  } else {
    out_calib <- list()
    out_calib$ageGrid <- x$`Age range lower (BP)` : x$`Age range upper (BP)`
    out_calib$densities <- rep(1 / (x$`Age range lower (BP)` - x$`Age range upper (BP)`), length(out_calib$ageGrid))
  }
  return(out_calib)
}

# A function to extract density estimates over a vector of ages (years)
extract_density <- function(dens, vect){
  out_extract <- rep(0,length(vect))
  out_extract[which(vect %in% dens$ageGrid)] <- dens$densities[which(dens$ageGrid %in% vect)]
  return(out_extract)
}

densities <- chronometric_data %>%
  dplyr::select(Site,
                Period,
                `14C age BP`,
                `1-sigma uncertainty`,
                `Age range lower (BP)`,
                `Age range upper (BP)`) %>%
  purrrlyr::by_row(calibrate_density, # Generate probability densities for 14C dates and age ranges
                   .to = "Density") %>%
  dplyr::mutate(Prediction = purrr::map(Density, # Extract probabilities for Marcott years
                                        extract_density,
                                        vect = marcott2013$YearBP)) %>%
  dplyr::select(Site, Period, Prediction) %>%
  dplyr::group_by(Site, Period) %>%
  purrrlyr::by_slice(map,~ Reduce(x = .x, f = "+"), .to = "Density") %>% # Sum probabilities over sites and periods (as in Oxcal SUM command)
  dplyr::mutate(Density = map(Density, "Prediction"),
                Density = map(Density, ~ .x / sum(.x, na.rm = T))) # re-normalize to 1


# Plot each site's probablility distribution
junk <- foreach(site = unique(densities$Site)) %do% {
  if(!file.exists(out("SITE_DENSITIES/",site,".pdf"))){
    cairo_pdf(filename = out("SITE_DENSITIES/",site,".pdf"))
    g <- densities %>%
      dplyr::filter(`Site` == site) %$%
      Density %>%
      magrittr::set_names(1:length(.)) %>%
      as_tibble() %>%
      mutate(`Years BP` = marcott2013$YearBP) %>%
      tidyr::gather(Period, Density, num_range("",1:(ncol(.)-1))) %>%
      ggplot2::ggplot(aes(x = `Years BP`,
                          y = Density,
                          colour = Period)) +
      ggplot2::geom_line() +
      ggplot2::xlim(6000,0) +
      ggplot2::ggtitle(site)
    print(g)
    dev.off()
  }
}

# # or, plot them all together with
# plot(1,type = "n",ylim = c(0,0.05), xlim = c(6000,1))
# for(site in unique(densities$`Site`)){
#   densities %>%
#     dplyr::filter(`Site` == site) %$%
#     Density %>%
#     magrittr::extract2(1) %>%
#     lines(x = marcott2013$YearBP,
#           y = .,
#           xlab = "Year BP",
#           ylab = "Density")
# }

# # Summing across all distributions yields the net distribution
# densities %$%
#   Density %>%
#   do.call(cbind,.) %>%
#   {rowSums(.)/ncol(.)} %>%
#   magrittr::multiply_by(20) %>%
#   plot(x = marcott2013$YearBP,
#        y = .,
#        xlim = c(6000,1),
#        type = "l",
#        xlab = "Year BP",
#        ylab = "Density")


# Create a spatial object of the sites
sites <- chronometric_data %>%
  mutate(foxtail_millet = ifelse(is.na(`Foxtail millet`), FALSE, TRUE),
         broomcorn_millet = ifelse(is.na(`Broomcorn millet`), FALSE, TRUE),
         wheat = ifelse(is.na(Wheat), FALSE, TRUE),
         barley = ifelse(is.na(Barley), FALSE, TRUE),
         buckwheat = ifelse(is.na(Buckwheat), FALSE, TRUE)) %>%
  dplyr::select(Site,
                Period,
                Longitude,
                Latitude,
                foxtail_millet,
                broomcorn_millet,
                wheat,
                barley,
                buckwheat) %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  tidyr::gather(Crop, Present, -Site, -Period, -Longitude, -Latitude) %>%
  dplyr::filter(Present) %>%
  dplyr::select(-Present) %>%
  dplyr::distinct() %>%
  dplyr::arrange(Site, Period) %>%
  sf::st_as_sf(coords = c("Longitude", "Latitude")) %>%
  sf::st_set_crs("+proj=longlat")

# Map the sites
cairo_pdf(out("FIGURES/crop_map.pdf"), width = 7.2, height = 7)
print(sites %>%
        dplyr::filter(Crop != "buckwheat") %>%
        dplyr::mutate(Crop = dplyr::recode_factor(Crop,
                                                  foxtail_millet = "Foxtail millet",
                                                  broomcorn_millet = "Broomcorn millet",
                                                  wheat = "Wheat",
                                                  barley = "Barley")) %>%
        ggplot() +
        geom_sf(data = countries) +
        geom_sf(aes(colour = Crop),
                size = 0.1) +
        theme_minimal(base_size = 7) +
        theme(legend.position="none") +
        ggplot2::facet_wrap( ~ Crop, ncol = 2))
dev.off()

# A function to extract niches for each crop
extract_niches <- function(x){
  out_niche <- raster::brick(out("RECONS/All_",x$Crop[[1]],"_Z.nc")) %>%
    raster:::readAll() %>%
    raster::extract(x %>%
                      as("Spatial")) %>%
    split(row(.))
  x %<>%
    set_class(c("tbl_df", "tbl", "data.frame")) %>%
    dplyr::mutate(Niche = out_niche) %>%
    sf::st_as_sf()
  return(x)
}

# Extract niches for each crop
niches <- sites %>%
  split(as.factor(sites$Crop)) %>%
  purrr::map(extract_niches) %>%
  purrr::map(~ set_class(., c("tbl_df", "tbl", "data.frame"))) %>%
  do.call(what = rbind, args = .) %>%
  sf::st_as_sf()

# A function to extract the 95% confidence interval of a distribution
get_CI <- function(dens){
  dens[is.na(dens)] <- 0
  cum_dens <- dens %>%
    cumsum() %>%
    magrittr::divide_by(sum(dens, na.rm = T))
  lower <- which(cum_dens >= 0.025)[1]
  mid <- which(cum_dens >= 0.5)[1]
  upper <- which(cum_dens >= 0.975)[1]
  return(list(Density = dens, Median = mid, CI = c(lower = lower, upper = upper)))
}

# A function to calculate local niches
local_niche <- function(niche, dens){
  if(is.na(dens$CI[["lower"]]) | is.na(dens$CI[["upper"]])) return(NULL)
  out_niche <- list()
  out_niche$Niche <- rep(NA, length(niche))
  out_niche$Niche[dens$CI[["lower"]]:dens$CI[["upper"]]] <- niche[dens$CI[["lower"]]:dens$CI[["upper"]]]
  out_niche$Niche <- out_niche$Niche / 100 # onvert to probability
  out_niche$Median <- quantile(out_niche$Niche, probs = 0.5, na.rm = TRUE)
  out_niche$CI <- c(quantile(out_niche$Niche, probs = 0.025, na.rm = TRUE),
                    quantile(out_niche$Niche, probs = 0.975, na.rm = TRUE))
  names(out_niche$CI) <- c("lower","upper")
  return(out_niche)
}

# Get the "Local" niche densities, by multiplying the site occupation probability densities
# by the crop niche probabilities (summing will get the average)
# Join the niche and density tables
niche_densities <- niches %>%
  set_class(c("tbl_df", "tbl", "data.frame")) %>%
  dplyr::left_join(densities, by = c("Site", "Period")) %>%
  dplyr::select(Site:Niche, Density) %>%
  # Get local densities
  dplyr::mutate(Density = purrr::map(Density, get_CI),
                Niche = purrr::map2(Niche, Density, local_niche)
  )

# Create biplots of each crop/site
junk <- foreach::foreach(crop = unique(niche_densities$Crop)) %do% {
  pdf(out("FIGURES/",crop,"_crossplot.pdf"))
  p <- niche_densities %>%
    dplyr::filter(Crop == crop) %>%
    dplyr::filter(!purrr::map_lgl(Niche, is.null)) %>%
    dplyr::mutate(Density_Median = marcott2013$YearBP[purrr::map_int(Density, "Median")],
                  Density_Lower = marcott2013$YearBP[purrr::map(Density, "CI") %>%
                                                       purrr::map_dbl("lower")],
                  Density_Upper = marcott2013$YearBP[purrr::map(Density, "CI") %>%
                                                       purrr::map_dbl("upper")],
                  Crop_Median = purrr::map_dbl(Niche, "Median"),
                  Crop_Lower = purrr::map(Niche, "CI") %>%
                    purrr::map_dbl("lower"),
                  Crop_Upper = purrr::map(Niche, "CI") %>%
                    purrr::map_dbl("upper")) %>%
    ggplot2::ggplot(aes(x = Density_Median, y = Crop_Median, label = Site)) +
    geom_point(na.rm = TRUE) +
    geom_errorbarh(aes(xmin = Density_Lower,
                       xmax = Density_Upper),
                   na.rm = TRUE) +
    geom_errorbar(aes(ymin = Crop_Lower,
                      ymax = Crop_Upper),
                  na.rm = TRUE) +
    xlim(6000,0) +
    ylim(0,1) +
    xlab("Years BP") +
    ylab(stringr::str_c("Probability of Being in the ",crop," Niche"))
  print(p)
  dev.off()
}

p <- niche_densities %>%
  dplyr::mutate(`Site, Period` = stringr::str_c(Site,ifelse(is.na(Period),"",stringr::str_c(", ",Period)))) %>%
  dplyr::filter(!purrr::map_lgl(Niche, is.null),
                Crop != "buckwheat") %>%
  dplyr::mutate(Crop = dplyr::recode_factor(Crop,
                                            foxtail_millet = "Foxtail millet",
                                            broomcorn_millet = "Broomcorn millet",
                                            wheat = "Wheat",
                                            barley = "Barley")) %>%
  dplyr::mutate(Density_Median = marcott2013$YearBP[purrr::map_int(Density, "Median")],
                Density_Lower = marcott2013$YearBP[purrr::map(Density, "CI") %>%
                                                     purrr::map_dbl("lower")],
                Density_Upper = marcott2013$YearBP[purrr::map(Density, "CI") %>%
                                                     purrr::map_dbl("upper")],
                Crop_Median = purrr::map_dbl(Niche, "Median"),
                Crop_Lower = purrr::map(Niche, "CI") %>%
                  purrr::map_dbl("lower"),
                Crop_Upper = purrr::map(Niche, "CI") %>%
                  purrr::map_dbl("upper")) %>%
  ggplot2::ggplot(aes(x = Density_Median,
                      y = Crop_Median,
                      colour = Crop,
                      label = `Site, Period`)) +
  geom_point(na.rm = TRUE) +
  geom_errorbarh(aes(xmin = Density_Lower,
                     xmax = Density_Upper),
                 na.rm = TRUE) +
  geom_errorbar(aes(ymin = Crop_Lower,
                    ymax = Crop_Upper),
                na.rm = TRUE) +
  xlim(6000,0) +
  ylim(0,1) +
  xlab("Years BP") +
  ylab(stringr::str_c("Probability of Being in the Niche"))

cairo_pdf(out("FIGURES/All_crossplot.pdf"), width = 7.2, height = 5.91)
print(p +
        theme_minimal(base_size = 7) +
        theme(legend.position="none") +
        ggplot2::facet_wrap( ~ Crop, ncol=2))
dev.off()

pp <- plotly::ggplotly(tooltip = c("label"))
htmlwidgets::saveWidget(widget = as_widget(pp),
                        file = "All_crossplot.html")
file.copy("All_crossplot.html",
          to = out("FIGURES/All_crossplot.html"),
          overwrite = TRUE)
unlink("All_crossplot.html")

# Write out a table
niche_densities_extract <- niche_densities %>%
  dplyr::filter(!purrr::map_lgl(Niche, is.null)) %>%
  dplyr::mutate(`Age median (BP)` = marcott2013$YearBP[purrr::map_int(Density, "Median")],
                `Age lower CI (BP)` = marcott2013$YearBP[purrr::map(Density, "CI") %>%
                                                           purrr::map_dbl("lower")],
                `Age upper CI (BP)` = marcott2013$YearBP[purrr::map(Density, "CI") %>%
                                                           purrr::map_dbl("upper")],
                `Niche probability median` = purrr::map_dbl(Niche, "Median"),
                `Niche probability lower CI` = purrr::map(Niche, "CI") %>%
                  purrr::map_dbl("lower"),
                `Niche probability upper CI` = purrr::map(Niche, "CI") %>%
                  purrr::map_dbl("upper")) %>%
  dplyr::select(-Niche:-Density) %>%
  dplyr::distinct() %>%
  dplyr::mutate(Crop = dplyr::recode_factor(Crop,
                                            foxtail_millet = "Foxtail millet",
                                            broomcorn_millet = "Broomcorn millet",
                                            wheat = "Wheat",
                                            barley = "Barley",
                                            buckwheat = "Buckwheat")) %>%
  dplyr::arrange(Site, Period, Crop) %T>%
  readr::write_csv(out("TABLES/age_niche_estimates.csv"))

# Plot relationship between number of crops and age
cairo_pdf(out("FIGURES/Age_counts.pdf"),
          width = 3.5,
          height = 2)
print({
  niche_densities_extract %>%
    dplyr::filter(Crop != "Buckwheat") %>%
    dplyr::select(-Crop) %>%
    dplyr::group_by(Site, Period) %>%
    dplyr::summarise_all(.funs = mean) %>%
    dplyr::left_join(niche_densities_extract %>%
                       dplyr::group_by(Site, Period) %>%
                       dplyr::count(),
                     by = c("Site", "Period")) %>%
    dplyr::rename(`Number of crops` = n) %>%
    ggplot() +
    geom_point(aes(x = `Number of crops`,
                   y = `Niche probability median`),
               size = 1) +
    geom_smooth(aes(x = `Number of crops`,
                    y = `Niche probability median`),
                method = "lm")+
    theme_minimal(base_size = 7)
})
dev.off()

##### END CHRONOMETRIC ANALYSIS #####



##### BEGIN PLOT CROP NICHE THROUGH TIME #####

## Plotting crop niche with associated sites
# create the cluster for parallel computation
message("Plotting crop niche reconstructions")
time_check <-  Sys.time()

# Combine the sites with the density data
site_densities <- niche_densities %>%
  dplyr::filter(!purrr::map_lgl(Niche, is.null)) %>%
  dplyr::mutate(Density_Lower = marcott2013$YearBP[purrr::map(Density, "CI") %>%
                                                     purrr::map_dbl("lower")],
                Density_Upper = marcott2013$YearBP[purrr::map(Density, "CI") %>%
                                                     purrr::map_dbl("upper")]) %>%
  dplyr::select(-Density) %>%
  dplyr::filter(!is.na(Density_Lower)) %>%
  dplyr::left_join(sites, by = c("Site", "Period", "Crop")) %>%
  dplyr::mutate(Niche = purrr::map(Niche, function(x){
    x$Niche <-
      tibble::tibble(`Year BP` = marcott2013$YearBP, Niche = x$Niche)
    return(x)
  }))

crops <- crop_GDD %>%
  dplyr::filter(crop %in% site_densities$Crop) %>%
  dplyr::select(crop_long, crop) %>%
  dplyr::distinct()

message("Generating niche videos.")
for(n in 1:nrow(crops)){

  crop <- crops[n,]$crop

  title <- stringr::str_c(crops[n,]$crop_long)

  if(file.exists(out("FIGURES/All_",crop,".pdf")) & file.exists(out("FIGURES/All_",crop,".mov")))
    next

  rast <- raster::brick(out("RECONS/All_",crop,"_Z.nc")) %>%
    magrittr::extract2(which(.@z$`Years BP` > 1000)) %>%
    raster:::readAll() %>%
    magrittr::divide_by(100) %>%
    magrittr::extract2(nlayers(.):1)

  rast.lower <- raster::brick(out("RECONS/All_",crop,"_Z_Lower.nc")) %>%
    magrittr::extract2(which(.@z$`Years BP` > 1000)) %>%
    raster:::readAll() %>%
    magrittr::divide_by(100) %>%
    magrittr::extract2(nlayers(.):1)

  rast.upper <- raster::brick(out("RECONS/All_",crop,"_Z_Upper.nc")) %>%
    magrittr::extract2(which(.@z$`Years BP` > 1000)) %>%
    raster:::readAll() %>%
    magrittr::divide_by(100) %>%
    magrittr::extract2(nlayers(.):1)

  years <- rast %>%
    names() %>%
    gsub(pattern = "X",
         replacement = "",
         x = .) %>%
    as.numeric()

  pal <- c(rev(colorRampPalette(brewer.pal(9, "Blues")[2:9],
                                bias = 2,
                                space = "Lab")(76)),
           colorRampPalette(brewer.pal(9, "Reds")[2:9],
                            bias = 1.5,
                            space = "Lab")(26))

  if(!file.exists(out("FIGURES/All_",crop,".pdf")))
    space_time_plot(the_brick = rast,
                    the_brick_lower = rast.lower,
                    the_brick_upper = rast.upper,
                    out_file = out("FIGURES/All_",crop,".pdf"),
                    title = title,
                    time = years,
                    timelim = c(max(years),min(years)),
                    timeaxis =  seq(from = max(years)-500,
                                    to = min(years),
                                    by = -500),
                    timelab = "Years BP",
                    zbreaks = seq(0,1,0.01),
                    zlab = "Probability of being in niche",
                    zaxis = seq(0,1,0.1),
                    zcolors = pal
    )


  if(!file.exists(out("FIGURES/All_",crop,".mov"))){
    # A function to extract site locations for a given year and crop, and plot them
    sites_plot <- function(years = NULL, crops = NULL, scale = scales::area_pal(c(0,1))){
      x_plot <- site_densities
      year_idx <- which(marcott2013$YearBP == years)

      if(!is.null(crops)){
        x_plot %<>%
          dplyr::filter(Crop %in% crops)
      }

      cexs <- purrr::map(x_plot$Niche,c("Niche","Niche")) %>%
        purrr::map(year_idx) %>%
        unlist() %>%
        scale()

      if(cexs %>%
         is.na() %>%
         all()) return() #Don't plot if no sites during this period

      x_plot %<>%
        dplyr::filter(!is.na(cexs))

      cexs <- cexs[!is.na(cexs)]

      x_plot %$%
        plot(geometry,
             cex = cexs,
             pch = 1,
             lwd = 3,
             col = "white",
             add = T)
      x_plot %$%
        plot(geometry,
             cex = cexs,
             pch = 1,
             lwd = 1.5,
             col = "black",
             add = T)
    }

    sites_legend <- function(){
      legend("bottomright",
             title = "Site probability \nof being in niche",
             legend = seq(0,1,0.2),
             pch = 1,
             pt.lwd = 1.5,
             col = "black",
             pt.cex = scales::area_pal(c(0.5,3))(seq(0,1,0.2)),
             bty = "n",
             y.intersp = 1.5,
             cex = 0.8,
             text.font = 2
      )
    }

    space_time_video(the_brick = rast,
                     the_brick_lower = rast.lower,
                     the_brick_upper = rast.upper,
                     out_file = out("FIGURES/All_",crop,".mov"),
                     title = title,
                     time = years,
                     timelim = c(max(years),min(years)),
                     timeaxis =  seq(from = max(years)-500,
                                     to = min(years),
                                     by = -500),
                     timelab = "Years BP",
                     zbreaks = seq(0,1,0.01),
                     zlab = "Probability of being in niche",
                     zaxis = seq(0,1,0.1),
                     zcolors = pal,
                     extra_plot_fun = purrr::partial(sites_plot,
                                                     crops = crop,
                                                     scale = scales::area_pal(c(0.5,3))),
                     extra_legend_fun = sites_legend
    )

  }
}

# Create a static plot for paper publication
# A function to extract data from a crop raster brick
tidyCrop <- function(x, years){
  x %>%
    out("RECONS/All_",.,"_Z.nc") %>%
    raster::brick() %>%
    magrittr::extract2(stringr::str_c("X",years)) %>%
    raster:::readAll() %>%
    magrittr::divide_by(100) %>%
    as("SpatialPixelsDataFrame") %>%
    as_tibble() %>%
    tidyr::gather(Year, Niche, -x:-y) %>%
    dplyr::rename(Longitude = x,
                  Latitude = y) %>%
    dplyr::mutate(Year = gsub(pattern = "X",
                              replacement = "",
                              x = Year),
                  Year = factor(Year, levels = years),
                  Year = forcats::fct_relabel(Year, function(x){stringr::str_c(x," cal. BP")}),
                  Crop = x)
}

# A color palette for filling
pal <- c(rev(colorRampPalette(brewer.pal(9, "Blues")[2:9],
                              bias = 2,
                              space = "Lab")(76)),
         colorRampPalette(brewer.pal(9, "Reds")[2:9],
                          bias = 1.5,
                          space = "Lab")(26))

# A function to create a plot of years and crops
facet_niche <- function(crops, years){
  these_sites <- years %>%
    map(function(year){
      site_densities %>%
        dplyr::mutate(Year = year) %>%
        dplyr::filter(Density_Lower <= Year & Density_Upper >= Year,
                      Crop %in% crops) %>%
        dplyr::mutate(Longitude = sf::st_coordinates(geometry)[,"X"],
                      Latitude = sf::st_coordinates(geometry)[,"Y"]) %>%
        dplyr::select(Site, Longitude, Latitude, Crop, Year)
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(Year = factor(Year, levels = years),
                  Year = forcats::fct_relabel(Year, function(x){stringr::str_c(x," cal. BP")}),
                  Crop = dplyr::recode_factor(Crop,
                                              wheat = "Wheat",
                                              barley = "Barley",
                                              broomcorn_millet = "Broomcorn millet",
                                              foxtail_millet = "Foxtail millet"))

  these_rasters <- crops %>%
    purrr::map(tidyCrop, years = years) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(Crop = dplyr::recode_factor(Crop,
                                              wheat = "Wheat",
                                              barley = "Barley",
                                              broomcorn_millet = "Broomcorn millet",
                                              foxtail_millet = "Foxtail millet"))

  p <- these_rasters %>%
    ggplot() +
    geom_raster(mapping = aes(x = Longitude,
                              y = Latitude,
                              fill = Niche)) +
    scale_fill_gradientn(colours = pal,
                         limits=c(0, 1),
                         name = "Probability of\nbeing in niche") +
    geom_point(data = these_sites,
               mapping = aes(x = Longitude,
                             y = Latitude),
               size = 1) +
    coord_quickmap() +
    facet_grid(Crop ~ Year,
               switch = "y")

  return(p)
}

cairo_pdf(out("FIGURES/facet_niche.pdf"), width = 7.2, height = 8.25)
print({
  facet_niche(crops = c("wheat", "barley", "broomcorn_millet", "foxtail_millet"),
              years = c(4030, 3550, 1690)) +
    scale_y_continuous(position = "right") +
    # xlab("") +
    # ylab("") +
    theme_minimal(base_size = 7) +
    theme(strip.text.x = element_text(size = 8),
          strip.text.y = element_text(size = 8))
})
dev.off()

message("Plotting of crop niche reconstructions complete: ", capture.output(Sys.time() - time_check))

##### END PREDICT CROP NICHE THROUGH TIME #####

# Output sessionInfo()
unlink(out("session_info.txt"))
devtools::session_info() %>%
  capture.output() %>%
  readr::write_lines(out("session_info.txt"))

# Output all package citations
unlink(out("packages.bib"))
bibtex::write.bib(all_packages, out("packages.bib"))

# Render the README.Rmd file
rmarkdown::render("README.Rmd")

message("asian_niche.R complete! Total run time: ", capture.output(Sys.time() - start_time))

### END SCRIPT ###