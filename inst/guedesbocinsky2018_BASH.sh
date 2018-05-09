#!/bin/bash
## Set the version as a local environment variable
VERSION="1.0.0"

## Set these parameters if not pre-loaded as environment variables
# google_maps_elevation_api_key=""
# tdar_un=""
# tdar_pw=""

### Build the vignette on the local machine
## Change 'clean' to TRUE to remove all previous output
Rscript -e "install.packages('devtools')"
Rscript -e "devtools::install_github('r-lib/devtools')"
Rscript -e "devtools::install(dependencies = TRUE, upgrade_dependencies = FALSE)"
Rscript -e "rmarkdown::render('./vignettes/guedesbocinsky2018.Rmd', \
                                  params = list(cores = 2, \
                                  clean = FALSE, \
                                  google_maps_elevation_api_key = '$google_maps_elevation_api_key', \
                                  tdar_un = '$tdar_un',\
                                  tdar_pw = '$tdar_pw'))"
