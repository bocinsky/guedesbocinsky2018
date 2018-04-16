#!/bin/bash
## Set the version as an environment variable
VERSION="0.9.0"
ARCH_SITES="../../DALPOIMGUEDES_BOCINSKY_2017.xlsx"

## Build the Docker image from the github repo
docker build -t bocinsky/guedesbocinsky2018 https://github.com/bocinsky/guedesbocinsky2018.git#$VERSION

## Remove any previous containers
docker rm guedesbocinsky2018

## Create the Docker container
docker create -w /guedesbocinsky2018 --name guedesbocinsky2018 bocinsky/guedesbocinsky2018

## Start the Docker container
docker start guedesbocinsky2018

## Copy the archaeological site data to the Docker container
docker cp $ARCH_SITES guedesbocinsky2018:/guedesbocinsky2018/DATA/guedesbocinsky2018.xlsx

## Download and copy pre-run output into the container
#docker cp ~/Desktop/OUTPUT guedesbocinsky2018:/guedesbocinsky2018/

## Run the analysis in the docker container
docker exec guedesbocinsky2018 Rscript guedesbocinsky2018.R

## Copy the output from the container to the host
docker cp guedesbocinsky2018:/guedesbocinsky2018/OUTPUT ./

## Copy README.md from the container to the host
docker cp guedesbocinsky2018:/guedesbocinsky2018/README.md ./README.md

## Stop the Docker container
docker stop guedesbocinsky2018

## Make a Zenodo directory
rm -r Zenodo; mkdir Zenodo

## Create a compressed tar archive of the output
tar -zcf ./Zenodo/guedesbocinsky2018-$VERSION-OUTPUT.tar.gz OUTPUT

## Make the Submission directory
rm -r Submission; mkdir Submission

## Copy and rename the figures, tables, and supplementary data sets for submission
cp ./OUTPUT/FIGURES/crop_map.pdf ./Submission/Figure_1.pdf
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
