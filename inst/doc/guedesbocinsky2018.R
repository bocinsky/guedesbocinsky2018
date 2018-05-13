params <-
list(cores = 2L, clean = FALSE, google_maps_elevation_api_key = "", 
    tdar_un = "rbocinsk", tdar_pw = "Dreams5683")

## ---- setup, echo = FALSE, cache = FALSE---------------------------------
# params <- list(cores = 2,
#                clean = FALSE)

# Set the knitting behavior
knitr::opts_chunk$set(echo = TRUE,
                      warning = FALSE,
                      message = FALSE,
                      collapse = TRUE, 
                      # cache = FALSE,
                      cache = !params$clean,
                      results = 'hold',
                      out.width = "100%",
                      fig.height = 5,  
                      comment = "#>",
                      fig.path = "./figures/"
)

library(guedesbocinsky2018)
library(sf)
library(magrittr)
library(fields)

# Set the behavior for printing scientific notation
# options(digits = 3, scipen = -2)
options(scipen = 999,
        knitr.table.format = "html")

# Force Raster to load large rasters into memory
raster::rasterOptions(chunksize = 2e+08,
                      maxmemory = 2e+09)

# Create the plan for (possibly parallel) execution
future::plan(future::multisession, 
             workers = min(future::availableCores(), 
                           params$cores),
             gc = TRUE)

# Create output directories
dir.create("./data/raw_data",
           recursive = TRUE,
           showWarnings = FALSE)
dir.create("./data/derived_data",
           recursive = TRUE,
           showWarnings = FALSE)
dir.create("./figures",
           recursive = TRUE,
           showWarnings = FALSE)


## ----plot-niches-through-time, cache = FALSE-----------------------------
## Plotting cultivar niche
# create the cluster for parallel computation
message("Plotting cultivar niche reconstructions")
time_check <-  Sys.time()

dir.create("./figures/cultivar_niches/",
           recursive = TRUE,
           showWarnings = FALSE)

gdd.recons <- 
  purrr::map_chr(1:nrow(crop_GDD),
                 function(n){
                   
                   cultivar <- crop_GDD[n,]$cultivar
                   
                   title <- stringr::str_c(crop_GDD[n,]$cultivar_long,
                                           " \u2014 T_base: ",
                                           crop_GDD[n,]$t_base,
                                           "°C, Required GDD: ",
                                           crop_GDD[n,]$min_gdd)
                   
                   if(file.exists(paste0("./figures/cultivar_niches/",cultivar,".pdf")) & file.exists(paste0("./figures/cultivar_niches/",cultivar,".mp4")))
                     return(paste0("./figures/cultivar_niches/",cultivar,".pdf"))
                   
                   rasts <- readr::read_rds(paste0("./data/derived_data/recons/",cultivar,"_recons.rds")) %>%
                     purrr::map(function(x){
                       x %>%
                         magrittr::extract2(which(.@z$`Years BP` > 1000)) %>%
                         # raster:::readAll() %>%
                         magrittr::extract2(raster::nlayers(.):1) 
                     })
                   
                   
                   years <- rasts[[1]] %>%
                     names() %>%
                     gsub(pattern = "X",
                          replacement = "",
                          x = .) %>%
                     as.numeric()
                   
                   breaks <- seq(0, 100, 10)
                   
                   pal <- (RColorBrewer::brewer.pal(10, "Spectral") %>%
                             rev() %>%
                             colorRampPalette(.))(length(breaks) - 1)
                   
                   if(!file.exists(paste0("./figures/cultivar_niches/",cultivar,".pdf")))
                     guedesbocinsky2018:::space_time_plot(
                       the_brick = rasts$Z,
                       the_brick_lower = rasts$Z_Lower,
                       the_brick_upper = rasts$Z_Upper,
                       out_file = paste0("./figures/cultivar_niches/",cultivar,".pdf"),
                       title = title,
                       time = years,
                       timelim = c(max(years),min(years)),
                       timeaxis =  seq(from = max(years)-500,
                                       to = min(years),
                                       by = -500),
                       timelab = "Years BP",
                       zbreaks = breaks,
                       zlab = "Probability of being in niche",
                       zaxis = seq(0,100,10),
                       zcolors = pal
                     )
                   
                   if(!file.exists(paste0("./figures/cultivar_niches/",cultivar,".mp4")))
                     guedesbocinsky2018:::space_time_video(                       the_brick = rasts$Z,
                                                                                  the_brick_lower = rasts$Z_Lower,
                                                                                  the_brick_upper = rasts$Z_Upper,
                                                                                  out_file = paste0("./figures/cultivar_niches/",cultivar,".mp4"),
                                                                                  title = title,
                                                                                  time = years,
                                                                                  timelim = c(max(years),min(years)),
                                                                                  timeaxis =  seq(from = max(years)-500,
                                                                                                  to = min(years),
                                                                                                  by = -500),
                                                                                  timelab = "Years BP",
                                                                                  zbreaks = breaks,
                                                                                  zlab = "Probability of being in niche",
                                                                                  zaxis = seq(0,100,10),
                                                                                  zcolors = pal
                     )
                   
                   return(paste0("./figures/cultivar_niches/",cultivar,".pdf"))
                   
                 })

message("Plotting of cultivar niche reconstructions complete: ", capture.output(Sys.time() - time_check))


## ----colophon, cache = FALSE---------------------------------------------
# which R packages and versions?
devtools::session_info()

