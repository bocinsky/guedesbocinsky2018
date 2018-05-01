# get the base image
FROM bocinsky/bocin_base:latest

# required
MAINTAINER Kyle Bocinsky <bocinsky@gmail.com>

COPY . /guedesbocinsky2018

# Install dev version of devtools to facilitate installing from "remotes" field in DESCRIPTION
RUN r -e 'devtools::install_github("r-lib/devtools")'

# Check the package
RUN r -e 'devtools::check(vignettes = FALSE, args = "--no-vignettes")'

# build this compendium package
RUN r -e 'devtools::install("/guedesbocinsky2018", dep = TRUE, upgrade_dependencies = FALSE)'

# render the analysis
# && r -e "rmarkdown::render('/guedesbocinsky2018/vignettes/guedesbocinsky2018.Rmd')"
