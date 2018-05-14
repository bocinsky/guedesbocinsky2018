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
                      cache.lazy = FALSE,
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


## ----submission, cache = FALSE-------------------------------------------
unlink("./submission",
       recursive = TRUE, 
       force = TRUE)

dir.create("./submission",
           recursive = TRUE,
           showWarnings = FALSE)

# Copy figures
file.copy(from = "./figures/Figure_1_crop_map.pdf",
          to = "./submission/Figure_1_crop_map.pdf",
          overwrite = TRUE)
file.copy(from = "./figures/Figure_2_facet_niche.pdf",
          to = "./submission/Figure_2_facet_niche.pdf",
          overwrite = TRUE)
file.copy(from = "./figures/Figure_3_all_crossplot.pdf",
          to = "./submission/Figure_3_all_crossplot.pdf",
          overwrite = TRUE)
file.copy(from = "./figures/Figure_4_storage.pdf",
          to = "./submission/Figure_4_storage.pdf",
          overwrite = TRUE)
file.copy(from = "./figures/Figure_5_modern_crops.pdf",
          to = "./submission/Figure_5_modern_crops.pdf",
          overwrite = TRUE)
file.copy(from = "./figures/Figure_6_modern_crops_density.pdf",
          to = "./submission/Figure_6_modern_crops_density.pdf",
          overwrite = TRUE)

# Copy Tables
file.copy(from = "./data/derived_data/Table_S1_crops.csv",
          to = "./submission/Table_S1_crops.csv",
          overwrite = TRUE)
file.copy(from = "./data/derived_data/Table_S2_age_niche_estimates.csv",
          to = "./submission/Table_S2_age_niche_estimates.csv",
          overwrite = TRUE)
file.copy(from = "./data/derived_data/Table_S3_storage.csv",
          to = "./submission/Table_S3_storage.csv",
          overwrite = TRUE)

# Copy Data Files
file.copy(from = "./data/derived_data/Data_file_S1_sites_dates.xlsx",
          to = "./submission/Data_file_S1_sites_dates.xlsx",
          overwrite = TRUE)
file.copy(from = "./figures/Data_file_S2_all_crossplot.html",
          to = "./submission/Data_file_S2_all_crossplot.html",
          overwrite = TRUE)

# Copy Movies
file.copy(from = "./figures/crop_niches/all_wheat.mp4",
          to = "./submission/Movie_S1_wheat.mp4",
          overwrite = TRUE)
file.copy(from = "./figures/crop_niches/all_barley.mp4",
          to = "./submission/Movie_S2_barley.mp4",
          overwrite = TRUE)
file.copy(from = "./figures/crop_niches/all_broomcorn_millet.mp4",
          to = "./submission/Movie_S3_broomcorn_millet.mp4",
          overwrite = TRUE)
file.copy(from = "./figures/crop_niches/all_foxtail_millet.mp4",
          to = "./submission/Movie_S4_foxtail_millet.mp4",
          overwrite = TRUE)
file.copy(from = "./figures/crop_niches/all_buckwheat.mp4",
          to = "./submission/Movie_S5_buckwheat.mp4",
          overwrite = TRUE)
file.copy(from = "./figures/crop_niches/all_rice.mp4",
          to = "./submission/Movie_S6_rice.mp4",
          overwrite = TRUE)

# Write the compendium
system(paste0("cp -r ../ ",tempdir(),"/Data_file_S3_research_compendium/"))

system(paste("tar -zcf submission/Data_file_S3_research_compendium.tar.gz",
             "-C ",tempdir(),
             "--exclude='vignettes/data'",
             "--exclude='vignettes/figures'",
             "--exclude='vignettes/submission'",
             "--exclude='vignettes/guedesbocinsky2018_cache'",
             "--exclude='vignettes/zenodo'",
             "--exclude='.git'",
             "Data_file_S3_research_compendium"))


## ----zenodo, cache = FALSE-----------------------------------------------
unlink("./zenodo",
       recursive = TRUE, 
       force = TRUE)

dir.create("./zenodo",
           recursive = TRUE,
           showWarnings = FALSE)

system("tar -zcf ./zenodo/guedesbocinsky2018-1.0.0-output.tar.gz submission figures data")

## ----colophon, cache = FALSE---------------------------------------------
# which R packages and versions?
devtools::session_info()

## ----git details, cache = FALSE------------------------------------------
# what commit is this file at? You may need to change the path value
# if your Rmd is not in analysis/paper/
git2r::repository("..")

