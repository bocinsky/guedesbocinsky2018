# get the base image
FROM bocinsky/bocin_base:latest

# required
MAINTAINER Kyle Bocinsky <bocinsky@gmail.com>

COPY . /guedesbocinsky2018

# go into the repo directory
RUN . /etc/environment \
  # Install dev version of devtools
  && R -e 'devtools::install_github("r-lib/devtools")' \
  # build this compendium package
  && R -e "devtools::install('/guedesbocinsky2018', dep = TRUE, upgrade_dependencies = FALSE)" \
  # render the analysis
  && R -e "rmarkdown::render('/guedesbocinsky2018/analysis/analysis.Rmd')"
