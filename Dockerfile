# get the base image, the rocker/verse has R, RStudio and pandoc
FROM bocinsky/bocin_base:latest

# required
MAINTAINER Kyle Bocinsky <bocinsky@gmail.com>

COPY . /guedesbocinsky2018

# go into the repo directory
RUN . /etc/environment \

  # build this compendium package
  && R -e "devtools::install('/guedesbocinsky2018', dep=TRUE)"

 # render the manuscript into a docx, you'll need to edit this if you've
 # customised the location and name of your main Rmd file
 # && R -e "rmarkdown::render('/guedesbocinsky2018/analysis/paper/paper.Rmd')"
