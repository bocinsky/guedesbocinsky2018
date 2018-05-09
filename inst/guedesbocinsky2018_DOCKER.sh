#!/bin/bash
## Set the version as a local environment variable
VERSION="1.0.0"

## Set these parameters if not pre-loaded as environment variables
# google_maps_elevation_api_key=""
# tdar_un=""
# tdar_pw=""

### Build the vignette in a docker container
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
