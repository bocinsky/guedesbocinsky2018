#!/bin/bash
## Set the version as a local environment variable
VERSION="1.0.0"

## Set these parameters if not pre-loaded as environment variables
# google_maps_elevation_api_key=""
# tdar_un=""
# tdar_pw=""

# ### Build the vignette on the local machine
# ## Change 'clean' to TRUE to remove all previous output
# Rscript -e "rmarkdown::render('./vignettes/guedesbocinsky2018.Rmd', \
#                                   params = list(cores = 2, \
#                                   clean = FALSE, \
#                                   google_maps_elevation_api_key = '$google_maps_elevation_api_key', \
#                                   tdar_un = '$tdar_un',\
#                                   tdar_pw = '$tdar_pw'))"

### Or, build the vignette in a docker container
## Remove any previous containers
docker rm guedesbocinsky2018
docker rmi guedesbocinsky2018

## Build the Docker image from the github repo
docker build -t guedesbocinsky2018 .

## Create the Docker container
docker create --name guedesbocinsky2018 guedesbocinsky2018

## Start the Docker container
docker start guedesbocinsky2018

## Build the vignette
docker exec guedesbocinsky2018 r -e "rmarkdown::render('/guedesbocinsky2018/vignettes/guedesbocinsky2018.Rmd', \
                                                        params = list(cores = 1, \
                                                                      clean = FALSE, \
                                                                      google_maps_elevation_api_key = '$google_maps_elevation_api_key', \
                                                                      tdar_un = '$tdar_un',\
                                                                      tdar_pw = '$tdar_pw'))"

## Make a Zenodo directory
rm -r docker_out; mkdir docker_out

## Copy the output from the container to the host
docker cp guedesbocinsky2018:/guedesbocinsky2018 ./docker_out/

## Stop the Docker container
docker stop guedesbocinsky2018

## Create a compressed tar archive of the output
tar -zcf ./docker_out/guedesbocinsky2018-$VERSION-output.tar.gz data ./docker_out/

## Make the Submission directory
rm -r Submission; mkdir Submission

## Copy and rename the figures, tables, and supplementary data sets for submission
cp ./figures/crop_map.pdf ./Submission/Figure_1.pdf
cp ./OUTPUT/FIGURES/facet_niche.pdf ./Submission/Figure_2.pdf
cp ./OUTPUT/FIGURES/All_crossplot.pdf ./Submission/Figure_3.pdf
cp ./OUTPUT/FIGURES/All_crossplot.html ./Submission/Supplementary_Data_1.html
cp ./OUTPUT/FIGURES/All_wheat.mov ./Submission/Supplementary_Video_1.mov
cp ./OUTPUT/FIGURES/All_barley.mov ./Submission/Supplementary_Video_2.mov
cp ./OUTPUT/FIGURES/All_broomcorn_millet.mov ./Submission/Supplementary_Video_3.mov
cp ./OUTPUT/FIGURES/All_foxtail_millet.mov ./Submission/Supplementary_Video_4.mov
cp ./OUTPUT/FIGURES/All_buckwheat.mov ./Submission/Supplementary_Video_5.mov
cp ./DATA/crops.csv ./Submission/Supplementary_Table_1.csv
cp ./OUTPUT/TABLES/sites_dates_raw.csv ./Submission/Supplementary_Table_2.csv
cp ./OUTPUT/TABLES/age_niche_estimates.csv ./Submission/Supplementary_Table_3.csv
curl -L https://codeload.github.com/bocinsky/guedesbocinsky2018/tar.gz/$VERSION
