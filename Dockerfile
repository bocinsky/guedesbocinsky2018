FROM bocinsky/bocin_base:latest

MAINTAINER Kyle Bocinsky <bocinsky@gmail.com>

## Install R package dependencies from stable MRAN repo
RUN install2.r --error \
    ## Packages for parallel processing
    foreach \
    doParallel \
    ## Packages for chronometric analysis
    Bchron \
    mclust \

## Add the build context to the root
ADD . /guedesbocinsky2018/
